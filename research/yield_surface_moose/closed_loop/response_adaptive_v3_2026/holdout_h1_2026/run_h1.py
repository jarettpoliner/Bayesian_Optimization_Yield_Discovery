#!/usr/bin/env python3
"""Holdout H1 campaigns (analytic Yld91 oracle; no MOOSE).  All arms: 4,096-axis library, start = population start, budget 64,
identical initial 8 (= frozen V2Controller.initial_design(start) = first 8 of frozen axial maximin; asserted).
  V3      frozen V3Controller (response_adaptive_v3_2026/v3_lrv.py, hash-checked), actual radii
  V2      frozen V2Controller (bomes_v2, A20, antipodal), actual radii
  AXMM    frozen bomes_v2.axial_maximin from first(start)
  RANDOM  initial 8, then numpy default_rng(9000 + start) permutation of the remaining axes
  V3_ISO  DIAGNOSTIC response ablation: frozen V3Controller observing ISOTROPIC von Mises radii (no information about the material);
          its design is re-observed with the actual radii for evaluation (evaluate_h1.py)
Checkpointed: a campaign with run_info.json is skipped; an incomplete directory is moved to runs_incomplete/ (never deleted).
    nice -n 19 python run_h1.py all [workers]     |     python run_h1.py <ARM> <i>"""
import os, sys
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "1")
import hashlib, json, shutil, time
from pathlib import Path
import numpy as np, pandas as pd
HERE = Path(__file__).resolve().parent; V3D = HERE.parent; CL = V3D.parent; CS = CL / "complexity_scaling_2026"; sys.dont_write_bytecode = True
for p in (str(HERE), str(V3D), str(CS)):
    sys.path.insert(0, p)
import cs_common as C                                                                # noqa: E402
import bomes_v2 as B                                                                 # noqa: E402  frozen
from yield_surface_discovery.oracle import AnalyticalOracle                          # noqa: E402  frozen V1
import v3_lrv as V3                                                                  # noqa: E402  frozen V3
import surfaces_h1 as H                                                              # noqa: E402
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
V3_FROZEN = dict(l.split()[::-1] for l in open(V3D / "PROTOCOL.sha256"))["v3_lrv.py"]; assert sha(V3D / "v3_lrv.py") == V3_FROZEN, "V3 must be the frozen v3.0"
POP = json.load(open(HERE / "population.json"))["realisations"]; BUDGET = 64; ARMS = ("V3", "V2", "AXMM", "RANDOM", "V3_ISO")
J2 = lambda D: 1.0 / np.sqrt(np.einsum("ij,jk,ik->i", np.atleast_2d(D), C.M_J2, np.atleast_2d(D)))


class LoggingOracle(AnalyticalOracle):
    def __init__(self, fn, name, log):
        super().__init__(fn, name); self.src = name; self.log = open(log, "w"); self.n = 0

    def evaluate_many(self, directions):
        out = super().evaluate_many(directions)
        for o in out:
            self.n += 1; self.log.write(json.dumps(dict(call=self.n, source=self.src, direction=[float(x) for x in o.requested_direction], yield_radius=o.yield_radius)) + "\n")
        self.log.flush(); return out


def campaign(arm, i):
    out = HERE / "runs" / arm / f"r{i:02d}"
    if (out / "run_info.json").exists():
        return
    if out.exists():
        dst = HERE / "runs_incomplete" / arm / f"r{i:02d}_{int(time.time())}"; dst.parent.mkdir(parents=True, exist_ok=True); shutil.move(str(out), str(dst))
    out.mkdir(parents=True); t0 = time.time(); rec = POP[i]; F = H.Yld91(rec); s = rec["start"]; lib = C.library()
    oracle = LoggingOracle(J2 if arm == "V3_ISO" else F.radius, f"H1_r{i:02d}_{'isotropic_J2' if arm == 'V3_ISO' else 'actual'}", out / "evaluations.jsonl")
    st = B.V2Settings(schedule="A20", symmetry="antipodal", budget=BUDGET)
    ctrl = (V3.V3Controller if arm in ("V3", "V3_ISO") else B.V2Controller)(lib, oracle, st); init = ctrl.initial_design(s)
    steps = None
    if arm in ("V3", "V2", "V3_ISO"):
        res = ctrl.run(s); seq = res["sequence"]; radii = [o.yield_radius for o in res["observations"]]; steps = res["steps"]
        np.save(out / "scores_f32.npy", np.array([np.where(np.isfinite(x), x, np.nan) for x in res["scores"]], np.float32))
    else:
        if arm == "AXMM":
            seq = [int(c) for c in B.axial_maximin(ctrl.C, ctrl.first(s), BUDGET)]
        else:
            rest = [c for c in np.random.default_rng(9000 + s).permutation(len(ctrl.C)) if c not in set(init)]; seq = list(init) + [int(c) for c in rest[: BUDGET - 8]]
        radii = [o.yield_radius for o in oracle.evaluate_many(ctrl.C[seq])]
    assert len(seq) == BUDGET and len(set(seq)) == BUDGET and list(seq[:8]) == list(init)
    if steps is not None:
        pd.DataFrame([{k: v for k, v in x.items() if k != "tie_set"} for x in steps]).to_csv(out / "steps.csv", index=False)
    json.dump(dict(axes=[int(c) for c in seq], radii=[float(r) for r in radii]), open(out / "design.json", "w"))
    json.dump(dict(arm=arm, realisation=i, start=s, budget=BUDGET, wall_s=time.time() - t0, n_oracle_calls=oracle.n, v3_sha256=sha(V3D / "v3_lrv.py"),
                   library_hash=C.library_hash(lib), settings=st.__dict__,
                   frozen_sources=C.file_hashes([C.CL / "v2/bomes_v2.py", *sorted((C.REL / "yield_surface_discovery").glob("*.py")), CS / "cs_common.py", V3D / "v3_lrv.py"])),
              open(out / "run_info.json", "w"), indent=1, default=str)
    print(f"done {arm} r{i:02d}: {time.time() - t0:.1f} s", flush=True)


def _job(a):
    campaign(*a)


if __name__ == "__main__":
    if sys.argv[1] == "all":
        import multiprocessing as mp
        with mp.Pool(int(sys.argv[2]) if len(sys.argv) > 2 else 4) as P:
            P.map(_job, [(a, i) for i in range(H.N_REAL) for a in ARMS], chunksize=1)
    else:
        campaign(sys.argv[1], int(sys.argv[2]))
