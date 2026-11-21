"""Closed-loop campaign runner (Python API).   from yield_surface_bo import run_campaign;  res = run_campaign(cfg)

cfg (dict, or YAML/JSON file via load_config):
  method: v2_mes | v3_weighted | v3_signed | angular_maximin | random           (default v2_mes)
  n_initial: 8          budget: 33          start: 0   (matched-start rule of the frozen studies: start 0 -> candidate 0, s>0 -> rng(7000+s))
  library: {builtin: tl1_signed_8192 | analytic_axes_4096 | b1b2_axes_207} | {csv: path}
  specimen:
    kind: analytic | replay_release | replay_store | moose_truss | external_command
    symmetry: none | antipodal          symmetry_validated: false     symmetry_evidence: "..."
    analytic:        surface: yld91_h1 | von_mises | cb2004 | python      realisation: 0     function: "module:fn"
    replay_release:  trials: path/to/trials.jsonl                         (B1/B2 release format)
    replay_store:    stores: [path, ...]                                  (TL1 trial-store format)
    moose_truss:     lattice: path.json (omit = built-in TL1)  moose_exe: path  env_prefix: path  threshold: 1e-4  solver: {...}
    external_command: command: "python my_sim.py {direction_json} {workdir}"  timeout_s: 3600
  evaluation: {references_csv: path (z1,z2,z3,yield_radius[,set]) , checkpoints: [8, 16, 33, 64]}   (evaluator only; never seen by a controller)
  limits: {cpu_hours: 2.0, stop_file: STOP}
  outdir: path
Resume: run again with the same cfg and outdir; stored observations are replayed exactly and the deterministic controller continues."""
import json, os, time, warnings
from pathlib import Path
import numpy as np
from . import _frozen as F
from . import core

SQ2 = np.sqrt(2.0)


def load_config(path):
    p = Path(path); txt = p.read_text()
    if p.suffix in (".yaml", ".yml"):
        import yaml; cfg = yaml.safe_load(txt)
    else:
        cfg = json.loads(txt)
    base = p.resolve().parent                                             # relative paths in a config file are relative to that file
    def fix(d, keys):
        for k in keys:
            if isinstance(d.get(k), str) and not d[k].startswith("@") and not Path(os.path.expanduser(d[k])).is_absolute():
                d[k] = str(base / d[k])
    fix(cfg, ["outdir"]); fix(cfg.get("library", {}), ["csv"]); sp = cfg.get("specimen", {})
    for sub in ("replay_release", "moose_truss", "analytic", "external_command"):
        fix(sp.get(sub, {}), ["trials", "lattice", "calibration", "moose_exe", "env_prefix"])
    if "stores" in sp.get("replay_store", {}):
        sp["replay_store"]["stores"] = [s if (s.startswith("@") or Path(s).is_absolute()) else str(base / s) for s in sp["replay_store"]["stores"]]
    fix(cfg.get("evaluation", {}), ["references_csv"])
    return cfg


def _resolve_store_alias(s):
    """'@campaign3' / '@v3transfer' -> the archived TL1 trial stores inside the release."""
    return {"@campaign3": str(F.C3D), "@v3transfer": str(F.V3X)}.get(s, s)


def build_oracle(cfg):
    sp = cfg["specimen"]; k = sp["kind"]
    if k == "analytic":
        a = sp.get("analytic", {}); return core.AnalyticOracle(a.get("surface", "yld91_h1"), a.get("realisation", 0), a.get("scale", 1.0), a.get("function"))
    if k == "replay_release":
        t = sp["replay_release"]["trials"]; t = t.replace("@release", str(F.PS / "release_B1_B2_BO_benchmark_20261119")); return core.ReleaseTrialsReplay(t)
    if k == "replay_store":
        return core.StoreReplay([_resolve_store_alias(s) for s in sp["replay_store"]["stores"]])
    if k == "moose_truss":
        from .moose_truss import TrussMooseOracle
        m = sp.get("moose_truss", {})
        return TrussMooseOracle(Path(cfg["outdir"]) / "moose_workspace", m.get("lattice"), m.get("moose_exe"), m.get("env_prefix"), m.get("solver"),
                                m.get("threshold", 1e-4), m.get("calibration"))
    if k == "external_command":
        e = sp["external_command"]; env = dict(YSBO_ROOT=str(F.ROOT), **(e.get("env") or {}))
        return core.ExternalCommandOracle(e["command"], Path(cfg["outdir"]) / "external_runs", e.get("timeout_s"), env)
    raise ValueError(f"unknown specimen kind {k!r}")


def _controller_and_design(M, lib, cfg, oracle):
    B = M["B"]; sym = cfg["specimen"].get("symmetry", "none"); meth = cfg.get("method", core.DEFAULT_LIVE_METHOD)
    st = B.V2Settings(symmetry=sym, schedule="A20", n_initial=int(cfg.get("n_initial", 8)), budget=int(cfg.get("budget", 33)), record_scores=True)
    cls = {"v2_mes": B.V2Controller, "v3_weighted": M["V3"].V3Controller, "v3_signed": M["V3S"].V3SController}.get(meth, B.V2Controller)
    return cls(lib, oracle, st), st


def _run_design(ctl, st, cfg, oracle, start):
    """angular maximin / random on the controller's own candidate set (axis representatives in validated antipodal mode)."""
    B = F.load()["B"]; meth = cfg["method"]; init = ctl.initial_design(start); n = st.budget; spares = max(16, n // 4)
    if meth == "angular_maximin":
        if st.symmetry == "antipodal":
            order = [int(c) for c in B.axial_maximin(ctl.C, ctl.first(start), min(n + spares, len(ctl.C)), st.maximin_tol_rad)]
        else:
            from yield_surface_discovery.initialization import angular_maximin
            order = [int(c) for c in angular_maximin(ctl.C, ctl.first(start), min(n + spares, len(ctl.C)), "deterministic", st.maximin_tol_rad)]
        assert order[:st.n_initial] == init, "maximin prefix must equal the shared initial design"
    else:
        rest = [int(c) for c in np.random.default_rng(9000 + start).permutation(len(ctl.C)) if int(c) not in set(init)]; order = list(init) + rest
    seq, obs, failed = [], [], []
    for c in order:
        if len(seq) == n:
            break
        o = oracle.evaluate(ctl.C[c])
        if o.success:
            seq.append(int(c)); obs.append(o)
        else:
            failed.append(int(c))
    if len(seq) < n:
        raise core.StopCampaign(f"{meth}: only {len(seq)} valid observations from {len(order)} candidates tried")
    return dict(sequence=seq, observations=obs, steps=[dict(test_number=i + 1, candidate=c) for i, c in enumerate(seq)], scores=[], failed=failed)


def fit_gp(directions, radii, symmetry):
    M = F.load(); B = M["B"]; sym = B.ANTIPODAL if symmetry == "antipodal" else B.NO_SYMMETRY
    with warnings.catch_warnings():
        warnings.simplefilter("ignore"); return B.GPv2(sym, B.SCHEDULES["A20"]).fit(np.asarray(directions, float), np.asarray(radii, float))


def fibonacci_sphere(n):
    i = np.arange(n) + 0.5; z = 1 - 2 * i / n; rr = np.sqrt(1 - z ** 2); ph = np.pi * (1 + 5 ** 0.5) * i
    return np.c_[rr * np.cos(ph), rr * np.sin(ph), z]


def predict(gp, D):
    """GP prediction of the yield radius (median exp(mu)) and the total predictive sd of log r at directions D."""
    mu, vl = gp.predict_log(np.asarray(D, float)); sd = np.sqrt(vl + gp.theta()[2])
    return np.exp(mu), sd


def run_campaign(cfg, quiet=False):
    """run (or resume) one campaign.  returns a result dict; writes outputs to cfg['outdir']."""
    t0 = time.process_time(); M = F.load(); out = Path(cfg["outdir"]); out.mkdir(parents=True, exist_ok=True)
    cfg = dict(cfg); cfg.setdefault("method", core.DEFAULT_LIVE_METHOD); cfg.setdefault("start", 0)
    lib = core.load_library(cfg.get("library", {"builtin": "tl1_signed_8192"})); rep = core.library_report(lib)
    comp = core.check_compat(cfg, lib, rep); (out / "compat_report.json").write_text(json.dumps(dict(ok=comp.ok, errors=comp.errors, warnings=comp.warnings, info=comp.info), indent=1, default=str))
    if not comp.ok:
        raise ValueError("configuration rejected:\n  " + "\n  ".join(comp.errors))
    prev = out / "config.json"
    if prev.exists():
        old = json.loads(prev.read_text()); keys = ("method", "n_initial", "budget", "start", "library", "specimen")
        if any(old.get(k) != cfg.get(k) for k in keys):
            raise ValueError(f"{out} holds a different campaign; use a new outdir (resume requires an identical configuration)")
    prev.write_text(json.dumps(cfg, indent=1, default=str))
    lim = cfg.get("limits", {}); inner = build_oracle(cfg)
    oracle = core.CheckpointOracle(inner, out / "observations.jsonl", lim.get("cpu_hours"), (out / lim["stop_file"]) if lim.get("stop_file") else None)
    ctl, st = _controller_and_design(M, lib, cfg, oracle); meth = cfg["method"]; start = int(cfg["start"]); status = "complete"; err = None
    try:
        res = ctl.run(start) if core.METHODS[meth]["kind"] == "bo" else _run_design(ctl, st, cfg, oracle, start)
    except core.StopCampaign as e:
        status, err, res = "stopped", str(e), None
    info = dict(status=status, stop_reason=err, method=meth, maturity=core.METHODS[meth]["maturity"], method_note=core.METHODS[meth]["note"],
                warnings=comp.warnings, n_new_observations=oracle.n_new, n_resumed_observations=oracle.n_reused, oracle_cpu_h=oracle.cpu_s / 3600,
                controller_cpu_s=time.process_time() - t0, frozen_sources=F.verify(strict=False), version=__import__("yield_surface_bo").__version__)
    if res is None:
        (out / "run_info.json").write_text(json.dumps(info, indent=1, default=str))
        if not quiet:
            print(f"campaign {status}: {err}  (resume by re-running the same configuration)")
        return dict(info=info)
    seq = res["sequence"]; D = np.asarray(ctl.C)[seq]; r = np.array([o.yield_radius for o in res["observations"]])
    import pandas as pd
    pd.DataFrame(dict(test_number=np.arange(1, len(seq) + 1), phase=["initial"] * st.n_initial + ["selected"] * (len(seq) - st.n_initial), candidate=seq,
                      z1=D[:, 0], z2=D[:, 1], z3=D[:, 2], sigma_xx=D[:, 0] * r, sigma_yy=D[:, 1] * r, sigma_xy=D[:, 2] * r / SQ2, yield_radius=r)).to_csv(out / "selections.csv", index=False)
    pd.DataFrame([{k: v for k, v in x.items() if k != "tie_set"} for x in res["steps"]]).to_csv(out / "steps.csv", index=False)
    if res.get("scores"):
        np.save(out / "acquisition_scores_f32.npy", np.array([np.where(np.isfinite(x), x, np.nan) for x in res["scores"]], np.float32))
    sym = cfg["specimen"].get("symmetry", "none"); gp = fit_gp(D, r, sym); G = fibonacci_sphere(4000); rg, sd = predict(gp, G)
    pd.DataFrame(dict(z1=G[:, 0], z2=G[:, 1], z3=G[:, 2], r_pred=rg, sd_log=sd, r_lo95=rg * np.exp(-1.96 * sd), r_hi95=rg * np.exp(1.96 * sd))).to_csv(out / "reconstruction_grid.csv", index=False)
    amp, ls, wn = gp.theta(); info.update(n_observations=len(seq), gp_hyperparameters=dict(amplitude=amp, length_scale=ls, white_noise=wn))
    err_rows = evaluate_errors(cfg, D, r, oracle.truth)
    if err_rows:
        pd.DataFrame(err_rows).to_csv(out / "error_curve.csv", index=False); info["error_curve"] = err_rows
    try:
        from .plots import plot_campaign
        plot_campaign(out, D, r, gp, oracle.truth, err_rows, cfg, st.n_initial)
    except Exception as e:                                                                   # plotting never invalidates a campaign
        info["plot_error"] = repr(e)
    (out / "run_info.json").write_text(json.dumps(info, indent=1, default=str))
    if not quiet:
        print(f"campaign complete: {meth} ({info['maturity']}), {len(seq)} observations ({oracle.n_new} new, {oracle.n_reused} resumed) -> {out}")
        for w in comp.warnings:
            print("  WARNING:", w)
    return dict(info=info, sequence=seq, directions=D, radii=r, gp=gp, steps=res["steps"], scores=res.get("scores"), controller=ctl, settings=st)


def evaluate_errors(cfg, D, r, truth=None):
    """reconstruction error vs n (refit of the frozen GPv2 on the first n observations).  Truth: analytic surface or a references CSV."""
    ev = cfg.get("evaluation", {}); sym = cfg["specimen"].get("symmetry", "none"); P = None
    if ev.get("references_csv"):
        import pandas as pd
        ref = pd.read_csv(ev["references_csv"]); ref = ref[ref.get("accepted", True) == True] if "accepted" in ref else ref   # noqa: E712
        if "set" in ref and ev.get("reference_set"):
            ref = ref[ref.set == ev["reference_set"]]
        P = ref[["z1", "z2", "z3"]].to_numpy(); rt = ref["yield_radius"].to_numpy()
    elif truth is not None:
        P = fibonacci_sphere(int(ev.get("n_references", 4000)))
        if sym == "antipodal":
            P = P[P[:, 2] >= 0]
        rt = truth(P)
    if P is None:
        return []
    cps = [n for n in ev.get("checkpoints", [8, 16, 33, 64, 128]) if n <= len(r)]
    rows = []
    for n in cps:
        g = fit_gp(D[:n], r[:n], sym); rp, sd = predict(g, P); e = np.abs(rp / rt - 1); z = np.log(rt / rp) / sd
        rows.append(dict(n=n, mean_rel_error=float(e.mean()), max_rel_error=float(e.max()), p95_rel_error=float(np.percentile(e, 95)), coverage95=float(np.mean(np.abs(z) <= 1.96))))
    return rows
