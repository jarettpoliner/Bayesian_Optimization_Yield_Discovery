"""TV2C 1.0 -- versioned compatibility wrapper: the FROZEN V2-A20 controller (closed_loop/v2/bomes_v2.py, sha256 fbfcf412...) on the
TLIB-4096 v1 candidate library and the TL1 MOOSE truss oracle.

Nothing in the controller is changed or re-implemented.  The wrapper only
  * loads a candidate library CSV with the frozen V1 `CandidateLibrary.from_csv` (no dedup needed: TLIB has no duplicates, verified at build);
  * checks the frozen source hashes (V2 file + V1 package vs its release MANIFEST) before use;
  * provides a BLIND oracle adapter (`TrussYieldOracle`): every MOOSE observation runs in a separate worker process; the controller receives
    only (yield_radius, achieved stress direction, success); exact (specimen, direction) reuse within this benchmark only; every request logged;
  * provides the non-adaptive baselines on the same candidate set with the frozen selection rules (angular/axial maximin, random unique);
  * exposes V2Settings unchanged.  symmetry='none' is required for TL1 (antipodal symmetry failed its pre-stated admissibility test).
Kernel (Matern-3/2), acquisition (max latent predictive variance = max entropy), ML-II schedule A20, GP maths: all from the frozen module.
"""
import hashlib, json, sys, time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import numpy as np

VERSION = "TV2C-1.0"
ROOT = Path(__file__).resolve().parents[1]; CL = ROOT.parent
for p in (CL / "v2", CL, ROOT / "oracle"):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))
import bomes_v2 as B                                                                 # noqa: E402  (frozen)
from bomes_v2 import V2Controller, V2Settings, axial_maximin, AxisLibrary            # noqa: E402,F401
from yield_surface_discovery.candidates import CandidateLibrary                      # noqa: E402  (frozen V1)
from yield_surface_discovery.initialization import angular_maximin                   # noqa: E402  (frozen V1)
from yield_surface_discovery import YieldOracle                                      # noqa: E402
from run_bo_moose import bo_source_hashes, _key                                      # noqa: E402  (frozen hash check + reuse key)

FROZEN_V2_SHA = "fbfcf412a63a636b2b097eec2b57918fb95e837d0ec0c10ac1ee45a9480e453d"


def check_frozen():
    h = hashlib.sha256(Path(B.__file__).read_bytes()).hexdigest(); v1 = bo_source_hashes()
    ok = h == FROZEN_V2_SHA and all(v["matches_release"] for v in v1.values())
    if not ok:
        raise RuntimeError(f"frozen sources changed: bomes_v2 {h}, V1 {[k for k, v in v1.items() if not v['matches_release']]}")
    return dict(bomes_v2=h, v1={k: v["sha256"] for k, v in v1.items()})


def load_library(csv_path):
    lib = CandidateLibrary.from_csv(csv_path, deduplicate=False); return lib


def first_index(n_lib, start):
    """frozen V1/V2 rule for the first initial test (V2Controller.first)."""
    return 0 if start == 0 else int(np.random.default_rng(7000 + start).choice(n_lib))


def maximin_design(lib, start, n_total, symmetry="none"):
    """response-independent angular maximin (frozen rules; same first test as the BO start)."""
    s = V2Settings(symmetry=symmetry)
    if symmetry == "antipodal":
        ax = AxisLibrary.from_library(lib); return axial_maximin(ax.directions, int(ax.class_of[first_index(len(lib.directions), start)]), n_total, s.maximin_tol_rad)
    return angular_maximin(lib.directions, first_index(len(lib.directions), start), n_total, "deterministic", s.maximin_tol_rad)


def random_design(n_cand, seed, n_total):
    """random unique candidates (frozen B1/B2 rule: default_rng(9000 + k).choice(n, n_total, replace=False))."""
    return [int(c) for c in np.random.default_rng(9000 + seed).choice(n_cand, n_total, replace=False)]


class TrussYieldOracle(YieldOracle):
    """blind adapter: MOOSE runs in worker processes (oracle/runner._one); controller sees only the public observation."""

    def __init__(self, cal, run_dir, log_path, reuse_logs=(), workers=4, specimen="TL1/oracle-levin-v2"):
        self.cal, self.run_dir, self.log, self.workers, self.specimen = cal, Path(run_dir), Path(log_path), workers, specimen
        self.table, self.n = {}, 0
        for p in reuse_logs:
            for l in open(p):
                r = json.loads(l)
                if r.get("specimen") == specimen and "reused_from" not in r:
                    self.table.setdefault(_key(r["direction"]), (str(p), r))

    def _obs(self, d, rec):
        from yield_surface_discovery.oracle import YieldObservation
        ok = bool(rec.get("accepted")) and np.isfinite(rec.get("yield_radius", np.nan)) and rec["yield_radius"] > 0
        sm = np.asarray(rec["stress_mandel"], float) if ok else np.full(3, np.nan)
        return YieldObservation(requested_direction=np.asarray(d, float), achieved_direction=sm / np.linalg.norm(sm) if ok else np.asarray(d, float),
                                yield_radius=float(rec["yield_radius"]) if ok else float("nan"), stress_at_yield=sm * np.array([1, 1, 1 / np.sqrt(2)]) if ok else sm,
                                success=ok, status="ok" if ok else "rejected", metadata=dict(trial_id=rec.get("trial_id")))

    def _log(self, rec):
        with open(self.log, "a") as fh:
            fh.write(json.dumps(rec, default=float) + "\n")

    def evaluate_many(self, directions):
        import runner as R
        D = np.atleast_2d(directions); out = [None] * len(D); todo = []
        for i, d in enumerate(D):
            k = _key(d)
            if k in self.table:
                src, r = self.table[k]; self._log(dict(r, reused_from=src)); out[i] = self._obs(d, r)
            else:
                self.n += 1; todo.append((i, d, self.run_dir / f"t{self.n:05d}"))
        if todo:
            recs = R.run_batch([(d, wd) for _, d, wd in todo], self.cal, workers=self.workers, label=f"campaign:{self.run_dir.name}")
            for (i, d, wd), r in zip(todo, recs):
                r = dict(r, specimen=self.specimen, trial_id=f"{self.run_dir.name}/{Path(wd).name}", direction=[float(x) for x in d])
                self._log(r); self.table[_key(d)] = (str(self.log), r); out[i] = self._obs(d, r)
        return out

    def evaluate(self, direction):
        return self.evaluate_many(np.atleast_2d(direction))[0]
