#!/usr/bin/env python3
"""Closed-loop Bayesian experimental design on the MOOSE plane-stress composite benchmark.

BO side (this process, `ysd` env): the UNMODIFIED yield_surface_discovery package (frozen method: 8 angular-maximin tests ->
direct Matern-3/2 GP on log yield radius -> maximum predictive-variance acquisition -> repeat), driven by its own
AdaptiveController.run().  This file only supplies a YieldOracle subclass and records outputs.

Simulator side (separate processes, `moose` env): closed_loop/moose_worker.py -> psbench.observe.observe(direction, specimen_id).
This process never imports psbench; the controller receives only the specimen id it names, the directions it requested, and the
public measured observation.

    python run_bo_moose.py --config configs/S2_blind.yaml --out results/S2_blind_run1
    python run_bo_moose.py --config configs/S2_blind.yaml --out results/S2_blind_replay --replay results/S2_blind_run1/trials.jsonl
"""
import os, sys
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")                     # the package's bit-reproducibility setting
import argparse, hashlib, json, subprocess, threading, time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np, pandas as pd

HERE = Path(__file__).resolve().parent
PLANE_STRESS = HERE.parent
BO_PKG = Path.home() / "projects/yield_surface_discovery_release/software/yield_surface_discovery"
sys.path.insert(0, str(BO_PKG / "src"))

import yield_surface_discovery as ysd
from yield_surface_discovery import YieldOracle, YieldObservation, CachedOracle, mandel_to_cauchy
from yield_surface_discovery.io import load_config, build_controller, environment_info, write_json
from yield_surface_discovery.reconstruction import Reconstruction
from yield_surface_discovery.geometry import icosphere, angle_rad

MICROMAMBA = os.environ.get("MAMBA_EXE", str(Path.home() / ".micromamba/bin/micromamba"))
MAMBA_ROOT = os.environ.get("MAMBA_ROOT_PREFIX", str(Path.home() / "micromamba"))
MARK = "@@OBSERVATION@@"
BO_CORE = ["gp.py", "acquisition.py", "initialization.py", "controller.py", "symmetry.py", "geometry.py", "candidates.py", "oracle.py",
           "reconstruction.py", "io.py", "__init__.py"]


def _key(d):
    return tuple(np.round(np.asarray(d, float), 12))


def to_yield_observation(d, obs):
    """public psbench observation -> package YieldObservation.  success only if MOOSE accepted the trial and r is finite > 0."""
    d = np.asarray(d, float)
    r = obs.get("yield_radius"); r = float("nan") if r is None else float(r)
    sm = obs.get("stress_mandel")
    sm = np.asarray(sm, float) if sm is not None and all(v is not None for v in sm) else np.full(3, np.nan)
    ok = bool(obs.get("accepted")) and np.isfinite(r) and r > 0 and np.all(np.isfinite(sm))
    ach = sm / np.linalg.norm(sm) if ok else d.copy()
    status = "ok" if ok else ("error: " + obs["error"] if obs.get("error") else "rejected_by_moose_checks")
    meta = {k: obs.get(k) for k in ("trial_id", "log_yield_radius", "converged", "accepted", "direction_error", "sigma_zz",
                                     "mean_matrix_eqps", "runtime_seconds", "macro_strain", "error")}
    return YieldObservation(d, ach, mandel_to_cauchy(sm), r if ok else float("nan"), ok, status, meta)


class PsbenchMooseOracle(YieldOracle):
    """YieldOracle backed by psbench.observe, one fresh worker process per trial (initial batch run concurrently)."""

    def __init__(self, specimen_id, run_name, log_path, max_parallel=8):
        self.specimen_id, self.run_name, self.log_path, self.max_parallel = specimen_id, run_name, Path(log_path), max_parallel
        self.n_calls, self.lock, self.records = 0, threading.Lock(), []

    def _run_one(self, d, call_no, batch):
        tag = f"{self.run_name}_t{call_no:02d}"
        cmd = [MICROMAMBA, "run", "-r", MAMBA_ROOT, "-n", "moose", "python", str(HERE / "moose_worker.py"), "--specimen", self.specimen_id,
               "--direction", *[repr(float(x)) for x in d], "--tag", tag]
        t0, t_wall0 = time.perf_counter(), time.time()
        p = subprocess.run(cmd, cwd=PLANE_STRESS, capture_output=True, text=True)
        wall = time.perf_counter() - t0
        line = next((l for l in p.stdout.splitlines()[::-1] if l.startswith(MARK)), None)
        obs = json.loads(line[len(MARK):]) if line else dict(trial_id=f"{self.specimen_id}/{tag}", yield_radius=None, accepted=False,
                                                                converged=False, error=f"worker exit {p.returncode}: {p.stderr[-2000:]}")
        rec = dict(call_no=call_no, batch=batch, tag=tag, specimen_id=self.specimen_id, direction=[float(x) for x in d],
                   started_unix=t_wall0, worker_wall_s=wall, worker_returncode=p.returncode, observation=obs)
        with self.lock:
            self.records.append(rec)
            with open(self.log_path, "a") as fh:
                fh.write(json.dumps(rec, default=float) + "\n")
        print(f"  [{time.strftime('%H:%M:%S')}] trial {call_no:2d} d={np.round(d, 4)}  r={obs.get('yield_radius')}  "
              f"accepted={obs.get('accepted')}  {wall:.1f}s", flush=True)
        return to_yield_observation(d, obs)

    def _next(self):
        with self.lock:
            self.n_calls += 1; return self.n_calls

    def evaluate(self, direction):
        return self._run_one(np.asarray(direction, float), self._next(), "sequential")

    def evaluate_many(self, directions):
        D = np.atleast_2d(directions); nos = [self._next() for _ in D]
        with ThreadPoolExecutor(min(self.max_parallel, len(D))) as ex:
            return list(ex.map(lambda a: self._run_one(a[0], a[1], "initial_batch"), zip(D, nos)))


class ResumeOracle(PsbenchMooseOracle):
    """resume a campaign: a requested direction already recorded (same specimen) is served from the record (logged with reused_from);
    anything else runs MOOSE.  Used only to continue an earlier run whose recorded prefix the controller reproduces exactly."""

    def __init__(self, specimen_id, run_name, log_path, reuse_paths, max_parallel=8):
        super().__init__(specimen_id, run_name, log_path, max_parallel); self.table = {}
        for p in reuse_paths:
            for l in open(p):
                r = json.loads(l)
                if r["specimen_id"] == specimen_id:
                    self.table.setdefault(_key(r["direction"]), (str(p), r))

    def _run_one(self, d, call_no, batch):
        k = _key(d)
        if k not in self.table:
            return super()._run_one(d, call_no, batch)
        src, r = self.table[k]
        rec = dict(r, call_no=call_no, batch=batch, reused_from=src, original_tag=r["tag"])
        with self.lock:
            self.records.append(rec)
            with open(self.log_path, "a") as fh:
                fh.write(json.dumps(rec, default=float) + "\n")
        print(f"  trial {call_no:2d} d={np.round(d, 4)} reused from {Path(src).parent.name}/{r['tag']}  r={r['observation'].get('yield_radius')}", flush=True)
        return to_yield_observation(d, r["observation"])


class ReplayOracle(YieldOracle):
    """re-serves recorded MOOSE observations (keyed by requested direction) so the BO side can be re-run without MOOSE."""

    def __init__(self, trials_path):
        self.table = {}
        for l in open(trials_path):
            r = json.loads(l); self.table[_key(r["direction"])] = r
        self.records = []

    def evaluate(self, direction):
        k = _key(direction)
        if k not in self.table:
            raise KeyError(f"replay: direction {direction} was never run (the BO sequence diverged from the recorded run)")
        r = self.table[k]; self.records.append(r); return to_yield_observation(np.asarray(direction, float), r["observation"])


def sha256(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def bo_source_hashes():
    man = pd.read_csv(BO_PKG.parents[1] / "MANIFEST.csv").set_index("relative_path")["sha256"].to_dict()
    out = {}
    for f in BO_CORE:
        rel = f"software/yield_surface_discovery/src/yield_surface_discovery/{f}"; h = sha256(BO_PKG / "src/yield_surface_discovery" / f)
        out[f] = dict(sha256=h, release_manifest_sha256=man.get(rel), matches_release=(h == man.get(rel)))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=str(HERE / "configs/S2_blind.yaml")); ap.add_argument("--out", required=True)
    ap.add_argument("--replay", help="trials.jsonl of a previous run: re-run the BO with recorded observations, no MOOSE")
    ap.add_argument("--max-parallel", type=int, default=8)
    ap.add_argument("--start", type=int, default=0, help="initialization start (0: candidate first_index; s>=1: seeded random first candidate)")
    ap.add_argument("--reuse", nargs="*", default=[], help="trials.jsonl files whose recorded observations may be reused (resume)")
    a = ap.parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    cfg = load_config(a.config); spec = cfg["oracle"]["specimen_id"]
    hashes = bo_source_hashes()
    if not all(v["matches_release"] for v in hashes.values()):
        raise SystemExit(f"BO source differs from the release manifest: {hashes}")

    if a.replay:
        base = ReplayOracle(a.replay)
    else:
        if (out / "trials.jsonl").exists():
            raise SystemExit(f"{out}/trials.jsonl exists: refusing to overwrite results (use a new --out, or --reuse it from a new --out)")
        base = (ResumeOracle(spec, out.name, out / "trials.jsonl", a.reuse, a.max_parallel) if a.reuse
                else PsbenchMooseOracle(spec, out.name, out / "trials.jsonl", a.max_parallel))
    oracle = CachedOracle(base)                                       # as the package's own build_oracle() does
    ctrl = build_controller(cfg, oracle=oracle)
    ctrl.settings.record_scores = True                                # output only: keep the full acquisition vector of every step
    L = ctrl.library

    print(f"BO package {ysd.__version__} ({BO_PKG}); library {L.name}: {len(L)} candidates; symmetry {ctrl.symmetry.name}; "
          f"budget {ctrl.settings.budget} ({ctrl.settings.n_initial} initial); specimen {spec}", flush=True)
    t0 = time.time()
    R = ctrl.run(start=a.start)                                             # <- the package's closed loop, unmodified
    wall = time.time() - t0
    assert "psbench" not in sys.modules and not any(m.startswith("psbench") for m in sys.modules), "BO process imported psbench"

    # ---------------- outputs ----------------
    tests = R.tests_frame(); tests["log_yield_radius"] = np.log(tests.yield_radius)
    meta = pd.DataFrame([o.metadata for o in R.observations])
    tests = pd.concat([tests, meta.add_prefix("moose_")], axis=1); tests.to_csv(out / "tests.csv", index=False)
    steps = pd.DataFrame(R.steps); steps.to_csv(out / "steps.csv", index=False)
    write_json([o.to_dict() for o in R.observations], out / "observations.json")
    sc = np.array(R.scores)                                           # (n_adaptive, n_candidates), -inf at tested/unavailable candidates
    pd.DataFrame(np.where(np.isfinite(sc), sc, np.nan).T, columns=[f"step{s['test_number']}" for s in R.steps if "score" in s][:len(sc)]
                 ).assign(candidate_index=np.arange(len(L)), label=L.labels, z1=L.directions[:, 0], z2=L.directions[:, 1], z3=L.directions[:, 2]
                          ).to_csv(out / "acquisition_scores.csv", index=False)

    # GP posterior after every update (n = 8 ... 13 real tests), refitted with the controller's own fit()
    V, F = icosphere(4); n_tot = len(R.observations); post = {}
    th = np.linspace(0, 2 * np.pi, 361)
    sec = {"sxy0": np.c_[np.cos(th), np.sin(th), 0 * th], "syy0": np.c_[np.cos(th), 0 * th, np.sin(th)], "sxx0": np.c_[0 * th, np.cos(th), np.sin(th)]}
    rows, resp = [], []
    for n in range(ctrl.settings.n_initial, n_tot + 1):
        gp = R.fit(ctrl, n); hp = gp.hyperparameters(); rec = Reconstruction(gp)
        mu_c, var_c = gp.predict_log(L.directions); mu_v, var_v = gp.predict_log(V)
        post[f"n{n}_cand_mean_log"], post[f"n{n}_cand_var_log"] = mu_c, var_c
        post[f"n{n}_sphere_mean_log"], post[f"n{n}_sphere_var_log"] = mu_v, var_v
        for name, Ds in sec.items():
            m, v = gp.predict_log(Ds); post[f"n{n}_sec_{name}_mean_log"], post[f"n{n}_sec_{name}_var_log"] = m, v
        Dn = L.directions[R.sequence[:n]]; fit_res = gp.predict_log(Dn, return_var=False) - np.log(R.radii()[:n])
        rows.append(dict(n_real_tests=n, **{k: v for k, v in hp.items() if k != "kernel"}, kernel=hp["kernel"], gp_mean_constant=gp.mean_,
                         max_abs_train_residual_log=float(np.abs(fit_res).max()), sphere_sd_log_max=float(np.sqrt(var_v.max())),
                         sphere_sd_log_mean=float(np.sqrt(var_v).mean())))
        if n < n_tot:                                                 # prediction at the NEXT test before it was measured
            j = R.sequence[n]; m_b, v_b = gp.predict_log(L.directions[[j]])
            gp2 = R.fit(ctrl, n + 1); m_a, v_a = gp2.predict_log(L.directions[[j]])
            resp.append(dict(test_number=n + 1, candidate_index=j, label=L.labels[j], measured_log_r=float(np.log(R.radii()[n])),
                             pred_log_r_before=float(m_b[0]), sd_log_before=float(np.sqrt(v_b[0])), pred_log_r_after=float(m_a[0]),
                             sd_log_after=float(np.sqrt(v_a[0])), z_score_before=float((np.log(R.radii()[n]) - m_b[0]) / np.sqrt(v_b[0])),
                             sphere_max_abs_change_log=float(np.abs(gp2.predict_log(V, return_var=False) - gp.predict_log(V, return_var=False)).max())))
        if n == n_tot:
            rec.export(out / f"surface_n{n}.obj", level=4)
    np.savez_compressed(out / "posterior.npz", sphere_vertices=V, sphere_faces=F, candidates=L.directions, section_theta=th,
                        **{f"section_dirs_{k}": v for k, v in sec.items()}, **post)
    pd.DataFrame(rows).to_csv(out / "gp_updates.csv", index=False); pd.DataFrame(resp).to_csv(out / "gp_response.csv", index=False)

    write_json(dict(config=cfg, specimen_id=spec, replay=a.replay, start=a.start, reuse=a.reuse, wall_s=wall, environment=environment_info(), bo_package=str(BO_PKG),
                    bo_version=ysd.__version__, bo_source_hashes=hashes, library=dict(name=L.name, size=len(L), removed_duplicates=list(L.removed_duplicates)),
                    settings=R.settings, symmetry=R.symmetry, sequence=R.sequence, labels=R.labels, oracle_calls=oracle.calls,
                    failed_events=[s for s in R.steps if "event" in s], psbench_imported_by_bo_process=False), out / "run_info.json")
    print(tests[["test_number", "phase", "label", "req_z1", "req_z2", "req_z3", "yield_radius", "status"]].to_string(index=False))
    print(f"wall {wall:.1f}s, oracle calls {oracle.calls}; outputs in {out}")


if __name__ == "__main__":
    main()
