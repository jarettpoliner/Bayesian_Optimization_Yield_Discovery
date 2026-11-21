"""Public building blocks: candidate libraries, method registry, oracles, compatibility validation and checkpointing.
All selection mathematics is delegated to the frozen controllers (bomes_v2.V2Controller, v3_lrv.V3Controller, v3s_lrv.V3SController,
V1 angular_maximin, bomes_v2.axial_maximin); nothing here changes a kernel, acquisition, schedule or tie rule."""
import json, os, shutil, subprocess, time
from dataclasses import dataclass, field
from pathlib import Path
import numpy as np
from . import _frozen as F

# ----------------------------------------------------------------------------------------------------------------- method registry
METHODS = {
    "v2_mes": dict(kind="bo", maturity="validated", label="V2-A20 BO+MES (maximum predictive entropy = maximum latent GP variance)",
                   note="Conservative default for live MOOSE experiments; validated on B1, B2 and TL1 (beats random; superiority over angular maximin not established)."),
    "v3_weighted": dict(kind="bo", maturity="research", label="V3 LOO-error-weighted predictive variance (frozen v3.0-LRV)",
                        note="Leading research acquisition: statistically supported gain over angular maximin and V2 on the independent analytic H1 Yld91 holdout. "
                             "Superiority over angular maximin has NOT been established on physical MOOSE specimens (TL1: no supported difference); "
                             "the V3_ISO ablation does not support attributing its gain to material-specific responses."),
    "v3_signed": dict(kind="bo", maturity="experimental", label="V3-Signed (v3.0-LRV-S1): signed LOO-risk weighting",
                      note="Experimental extension for signed loading directions; not a validated replacement for frozen V3 (TL1: worse than frozen V3 at n = 33)."),
    "angular_maximin": dict(kind="design", maturity="baseline", label="angular maximin (signed angle; axial angle only in validated antipodal mode)",
                            note="Geometric space-filling baseline; never beaten with statistical support on a physical MOOSE specimen."),
    "random": dict(kind="design", maturity="baseline", label="random unique candidates (benchmarking only)", note="Benchmark baseline."),
}
DEFAULT_LIVE_METHOD = "v2_mes"


def method_info(name):
    if name not in METHODS:
        raise ValueError(f"unknown method {name!r}; choose one of {sorted(METHODS)}")
    return METHODS[name]


# ----------------------------------------------------------------------------------------------------------------- libraries
BUILTIN_LIBRARIES = {
    "tl1_signed_8192": "8,192 signed directions (TLIB-4096 v1: 4,096 axes and their antipodes), TL1 campaigns",
    "analytic_axes_4096": "4,096 equal-area half-sphere axes (cs_common.library), analytic studies incl. H1; antipodal mode only",
    "b1b2_axes_207": "V1 expanded physical loading library (414 directions = 207 axes), B1/B2 campaigns",
}


def load_library(spec):
    """spec: {'builtin': name} or {'csv': path} (V1 CSV format: id,label,family,z1,z2,z3).  returns a frozen V1 CandidateLibrary."""
    M = F.load(); CL_ = M["candidates"].CandidateLibrary
    if "csv" in spec:
        return CL_.from_csv(Path(spec["csv"]).expanduser(), deduplicate=False)
    name = spec.get("builtin")
    if name == "tl1_signed_8192":
        return CL_.from_csv(F.TRUSS / "library" / "TLIB4096_v1.csv", deduplicate=False)
    if name == "analytic_axes_4096":
        import cs_common
        return cs_common.library()
    if name == "b1b2_axes_207":
        from yield_surface_discovery.io import load_config, build_library
        return build_library(load_config(F.CL / "configs" / "S2_campaign.yaml"))
    raise ValueError(f"unknown library {spec}; builtins: {sorted(BUILTIN_LIBRARIES)}")


def library_report(lib, tol=1e-9):
    """uniqueness, unit norm and antipodal closure of a candidate library (O(n^2) memory-safe in blocks)."""
    D = np.asarray(lib.directions, float); n = len(D); norm_dev = float(np.abs(np.linalg.norm(D, axis=1) - 1).max())
    dup, anti_present, min_sep = 0, 0, np.inf
    for i0 in range(0, n, 1024):
        blk = D[i0:i0 + 1024]; dd = np.linalg.norm(blk[:, None, :] - D[None, :, :], axis=2); idx = np.arange(i0, min(i0 + 1024, n))
        dd[np.arange(len(idx)), idx] = np.inf; dup += int((dd < tol).sum()); min_sep = min(min_sep, float(dd.min()))
        anti_present += int((np.linalg.norm(blk[:, None, :] + D[None, :, :], axis=2) < tol).any(1).sum())
    return dict(n=n, unit_norm_max_dev=norm_dev, duplicate_pairs=dup // 2, min_chord_separation=min_sep,
                min_angle_deg=float(np.degrees(2 * np.arcsin(min(min_sep / 2, 1.0)))), antipode_fraction=anti_present / n,
                antipodally_closed=anti_present == n)


# ----------------------------------------------------------------------------------------------------------------- compatibility
@dataclass
class Compat:
    errors: list = field(default_factory=list)
    warnings: list = field(default_factory=list)
    info: dict = field(default_factory=dict)

    @property
    def ok(self):
        return not self.errors


def check_compat(cfg, lib=None, lib_rep=None):
    """validate a campaign configuration BEFORE any oracle call.  never applies a symmetry the specimen has not declared validated."""
    c = Compat(); sp = cfg.get("specimen", {}); sym = sp.get("symmetry", "none"); meth = cfg.get("method", DEFAULT_LIVE_METHOD)
    live = sp.get("kind") in ("moose_truss", "external_command")
    if meth not in METHODS:
        c.errors.append(f"unknown method {meth!r} (choose from {sorted(METHODS)})"); return c
    c.info["method"] = dict(name=meth, **METHODS[meth])
    if sym not in ("none", "antipodal"):
        c.errors.append(f"symmetry must be 'none' or 'antipodal', got {sym!r}")
    if sym == "antipodal" and not sp.get("symmetry_validated", False):
        c.errors.append("antipodal symmetry requested but specimen.symmetry_validated is not true: r(d) = r(-d) must be verified for this specimen "
                        "(e.g. paired +/-d measurements within a stated tolerance) before antipodal pseudo-observations or axial distances may be used")
    if sym == "antipodal" and sp.get("symmetry_validated") and not sp.get("symmetry_evidence"):
        c.warnings.append("symmetry_validated is true but no specimen.symmetry_evidence was given; record how r(d) = r(-d) was verified")
    if meth == "v3_signed" and sym == "antipodal":
        c.errors.append("v3_signed is defined for signed campaigns (symmetry 'none'); use v3_weighted or v2_mes with validated antipodal symmetry")
    if meth == "v3_weighted" and sym == "none":
        c.warnings.append("frozen v3_weighted keeps its antipodally symmetric LOO-risk weighting (rho(d) = rho(-d)) even with a signed GP; this is the "
                          "documented TL1 transfer-test configuration, not a claim that the specimen is symmetric (consider v2_mes, or v3_signed as an experiment)")
    if live and meth in ("v3_weighted", "v3_signed"):
        c.warnings.append(f"{meth} for live MOOSE: superiority over angular maximin has not been established on physical specimens; v2_mes is the default")
    if live and meth == "random":
        c.warnings.append("random sampling is a benchmarking baseline")
    n0, nb = int(cfg.get("n_initial", 8)), int(cfg.get("budget", 33))
    if not 1 <= n0 < nb:
        c.errors.append(f"need 1 <= n_initial < budget (got {n0}, {nb})")
    if lib is not None:
        rep = lib_rep or library_report(lib); c.info["library"] = rep
        if rep["duplicate_pairs"]:
            c.errors.append(f"candidate library has {rep['duplicate_pairs']} duplicate direction pairs")
        if rep["unit_norm_max_dev"] > 1e-9:
            c.errors.append("candidate directions are not unit vectors")
        if sym == "none" and rep["antipode_fraction"] < 0.5:
            c.errors.append("signed campaign (symmetry 'none') with a half-sphere axis library: half of the loading space would be unreachable; "
                            "supply both signs of every direction (e.g. builtin tl1_signed_8192)")
        n_cand = rep["n"] if sym == "none" else int(round(rep["n"] - rep["antipode_fraction"] * rep["n"] / 2))
        c.info["n_candidates"] = n_cand
        if nb > n_cand:
            c.errors.append(f"budget {nb} exceeds the number of candidates {n_cand}")
    return c


# ----------------------------------------------------------------------------------------------------------------- observations
def _obs(d, radius, success, stress_mandel=None, status=None, meta=None):
    from yield_surface_discovery.oracle import YieldObservation
    d = np.asarray(d, float); ok = bool(success) and radius is not None and np.isfinite(radius) and radius > 0
    sm = np.asarray(stress_mandel, float) if (ok and stress_mandel is not None) else (d * radius if ok else np.full(3, np.nan))
    return YieldObservation(requested_direction=d, achieved_direction=sm / np.linalg.norm(sm) if ok else d, stress_at_yield=sm * np.array([1, 1, 1 / np.sqrt(2)]),
                            yield_radius=float(radius) if ok else float("nan"), success=ok, status=status or ("ok" if ok else "rejected"), metadata=meta or {})


def dkey(d):
    d = np.asarray(d, float); d = d / np.linalg.norm(d); return tuple(np.round(d, 12).tolist())


class StopCampaign(RuntimeError):
    pass


class BaseOracle:
    truth = None                                                         # optional analytic radius function (evaluation only)

    def evaluate(self, d):
        raise NotImplementedError

    def evaluate_many(self, X):
        return [self.evaluate(x) for x in np.atleast_2d(X)]


class AnalyticOracle(BaseOracle):
    """analytic yield surfaces (no MOOSE): von_mises, cb2004 (V1 package), yld91_h1 (frozen H1 holdout realisation), or 'module:function'."""

    def __init__(self, surface="yld91_h1", realisation=0, scale=1.0, function=None):
        F.load(); self.name = surface
        if surface == "von_mises":
            from yield_surface_discovery.analytical import von_mises_radius as fn; self.fn = lambda D: scale * fn(D)
        elif surface == "cb2004":
            from yield_surface_discovery.analytical import cb2004_radius as fn; self.fn = fn
        elif surface == "yld91_h1":
            import run_h1 as RH                                           # frozen H1 population and Yld91 surfaces (imports only)
            self.rec = RH.POP[int(realisation)]; self.fn = RH.H.Yld91(self.rec).radius; self.name = f"yld91_h1_r{int(realisation):02d}"
        elif surface == "python":
            mod, attr = function.split(":"); import importlib; self.fn = getattr(importlib.import_module(mod), attr)
        else:
            raise ValueError(f"unknown analytic surface {surface!r}")
        self.truth = lambda D: np.atleast_1d(self.fn(np.atleast_2d(D)))

    def evaluate(self, d):
        r = float(self.truth(np.atleast_2d(d))[0]); return _obs(d, r, True, meta=dict(source=self.name))


class ReleaseTrialsReplay(BaseOracle):
    """recorded real MOOSE observations in the B1/B2 release format (trials.jsonl lines with 'direction' and 'observation').
    Uses the frozen run_bo_moose.ReplayOracle conversion (to_yield_observation), so replays are identical to the original."""

    def __init__(self, path):
        F.load(); import run_bo_moose as RB; self.inner = RB.ReplayOracle(Path(path)); self.calls = 0

    def evaluate(self, d):
        self.calls += 1
        try:
            return self.inner.evaluate(d)
        except KeyError as e:
            raise StopCampaign(f"replay diverged: {e}")


class StoreReplay(BaseOracle):
    """recorded real MOOSE observations in the TL1 trial-store format (trials/index.jsonl + trials/records/*.json), one or more stores."""

    def __init__(self, stores):
        self.table = {}; self.calls = 0
        for st in stores:
            st = Path(st)
            for l in (st / "trials" / "index.jsonl").read_text().splitlines():
                x = json.loads(l); r = json.loads((st / x["record"]).read_text()); self.table.setdefault(tuple(r["key"]), r)

    def evaluate(self, d):
        self.calls += 1; r = self.table.get(dkey(d))
        if r is None:
            raise StopCampaign(f"replay diverged: direction {np.round(d, 6).tolist()} was never measured")
        ok = bool(r.get("accepted")) and bool(r.get("valid"))
        return _obs(d, r.get("yield_radius"), ok, r.get("stress_mandel"), meta=dict(keyhash=r.get("keyhash"), cpu_s=r.get("cpu_s", 0.0)))


class CheckpointOracle(BaseOracle):
    """persists every observation (JSONL) and serves stored ones on resume; enforces a CPU-hour cap and a STOP file before any NEW run.
    Resume = re-run the deterministic controller from the start: stored observations are replayed exactly, new ones are measured."""

    def __init__(self, inner, path, cpu_cap_h=None, stop_file=None):
        self.inner, self.path, self.cap, self.stop = inner, Path(path), cpu_cap_h, (Path(stop_file) if stop_file else None)
        self.store, self.n_new, self.n_reused, self.cpu_s = {}, 0, 0, 0.0; self.truth = getattr(inner, "truth", None)
        if self.path.exists():
            for l in self.path.read_text().splitlines():
                x = json.loads(l); self.store[tuple(x["key"])] = x; self.cpu_s += float(x.get("cpu_s", 0.0))

    def evaluate(self, d):
        k = dkey(d)
        if k in self.store:
            x = self.store[k]; self.n_reused += 1
            return _obs(d, x["yield_radius"], x["success"], x.get("stress_mandel"), x.get("status"), dict(x.get("meta", {}), resumed=True))
        if self.stop is not None and self.stop.exists():
            raise StopCampaign(f"stop file {self.stop} present")
        if self.cap is not None and self.cpu_s / 3600 >= self.cap:
            raise StopCampaign(f"CPU cap {self.cap} h reached ({self.cpu_s / 3600:.3f} h)")
        o = self.inner.evaluate(d); cpu = float((o.metadata or {}).get("cpu_s", 0.0)); self.cpu_s += cpu; self.n_new += 1
        rec = dict(key=list(k), direction=[float(v) for v in np.asarray(d, float)], yield_radius=None if not o.success else float(o.yield_radius), success=bool(o.success),
                   stress_mandel=None if not o.success else [float(v) for v in np.asarray(o.stress_at_yield) * np.array([1, 1, np.sqrt(2)])], status=o.status,
                   cpu_s=cpu, meta={k_: v for k_, v in (o.metadata or {}).items() if isinstance(v, (str, int, float, bool, type(None)))}, time=time.time())
        self.path.parent.mkdir(parents=True, exist_ok=True)
        with open(self.path, "a") as f:
            f.write(json.dumps(rec) + "\n")
        self.store[k] = rec; return o


class ExternalCommandOracle(BaseOracle):
    """generic live adapter for ANY simulation (e.g. a user's MOOSE model): runs a command per direction in its own work directory.
    The command template may use {d0} {d1} {d2} (unit Mandel direction), {direction_json} and {workdir}; it must write <workdir>/result.json
    with {"yield_radius": float, "accepted": bool, optional "stress_mandel": [3], optional "cpu_s": float}."""

    def __init__(self, command, workroot, timeout_s=None, env=None):
        self.cmd, self.root, self.timeout, self.env, self.n = command, Path(workroot), timeout_s, env or {}, 0

    def evaluate(self, d):
        self.n += 1; d = np.asarray(d, float); wd = self.root / f"trial_{self.n:05d}"
        while wd.exists():
            self.n += 1; wd = self.root / f"trial_{self.n:05d}"
        wd.mkdir(parents=True)
        cmd = self.cmd.format(d0=repr(float(d[0])), d1=repr(float(d[1])), d2=repr(float(d[2])), direction_json=json.dumps([float(v) for v in d]), workdir=str(wd))
        t0 = time.time(); p = subprocess.run(cmd, shell=True, cwd=wd, capture_output=True, text=True, timeout=self.timeout, env=dict(os.environ, **self.env))
        (wd / "command.log").write_text(cmd + "\n" + p.stdout + "\n--- stderr ---\n" + p.stderr)
        res = json.loads((wd / "result.json").read_text()) if (wd / "result.json").exists() else dict(accepted=False, error=f"exit {p.returncode}, no result.json")
        return _obs(d, res.get("yield_radius"), bool(res.get("accepted")), res.get("stress_mandel"),
                    meta=dict(workdir=str(wd), cpu_s=float(res.get("cpu_s", time.time() - t0)), error=res.get("error")))
