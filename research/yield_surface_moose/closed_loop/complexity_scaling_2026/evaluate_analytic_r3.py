#!/usr/bin/env python3
"""Evaluator r3 (final closeout): extends evaluation_r2 (n <= 512) to the full protocol grid up to n = 1000.
EXACT same code path: evaluate_analytic_r2.task / design are re-used unchanged with only NMAX = 1000, the grid restricted to the
points > 512 (640, 768, 896, 1000 for A20; 1000 for FIXED) and a separate output directory.  BO designs are evaluated only once the
campaign has >= 1000 recorded observations (design() returns None before that; re-run to pick them up).
Rows for n <= 512 are taken from evaluation_r2/designs (bitwise identical to the original evaluator on the overlap; see
evaluation_r3/consistency.json).  Merged per-design CSVs -> evaluation_r3/designs/, all rows -> evaluation_r3/metrics_all.csv.
    python evaluate_analytic_r3.py [PAR] [ARMS...]"""
import os, sys
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "1")
import json, multiprocessing as mp
from pathlib import Path
import pandas as pd
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import cs_common as C                                                              # noqa: E402
import evaluate_analytic_r2 as R2                                                  # noqa: E402
R2.NMAX = 1000; R2.GRID = (640, 768, 896, 1000); R2.BUD = (1000,)
EXT = HERE / "evaluation_r3" / "ext"; EXT.mkdir(parents=True, exist_ok=True); R2.OUT = EXT
MERGED = HERE / "evaluation_r3" / "designs"; MERGED.mkdir(parents=True, exist_ok=True); ORIG = HERE / "evaluation" / "designs"


def seed_from_original():
    """the 6 designs the ORIGINAL evaluator (evaluate_analytic.py) already took to 1000: reuse its n > 512 rows after checking that its
    n <= 512 rows equal evaluation_r2 bitwise (same code path, deterministic)."""
    rep = {}
    for f in sorted(ORIG.glob("*.csv")):
        o, r2 = pd.read_csv(f), HERE / "evaluation_r2" / "designs" / f.name
        if not r2.exists():
            continue
        b = pd.read_csv(r2); m = o.merge(b, on=["surrogate", "n"], suffixes=("_o", "_r2")); num = [c for c in b.columns if c not in ("level", "arm", "start", "surrogate", "n", "at_bounds", "gp_mode", "fit_s")]
        same = bool(len(m) == len(b) and all((m[c + "_o"] == m[c + "_r2"]).all() for c in num))
        rep[f.name] = dict(overlap_rows=len(m), bitwise_equal=same)
        if same and not (EXT / f.name).exists():
            o[o.n > 512].to_csv(EXT / f.name, index=False)
    json.dump(rep, open(HERE / "evaluation_r3" / "consistency.json", "w"), indent=1)


def merge():
    n = 0
    for f in sorted(EXT.glob("*.csv")):
        r2 = HERE / "evaluation_r2" / "designs" / f.name
        if r2.exists():
            pd.concat([pd.read_csv(r2), pd.read_csv(f)]).sort_values(["surrogate", "n"]).to_csv(MERGED / f.name, index=False); n += 1
    fs = sorted(MERGED.glob("*.csv"))
    if fs:
        M = pd.concat([pd.read_csv(f) for f in fs]); M.to_csv(HERE / "evaluation_r3" / "metrics_all.csv", index=False); print(len(M), "rows from", len(fs), "designs", flush=True)


if __name__ == "__main__":
    par = int(sys.argv[1]) if len(sys.argv) > 1 else 3; arms = sys.argv[2:] or ["AXMM", "RANDOM", "BOFIXED", "BO"]
    seed_from_original()
    jobs = [(L, a, k) for L in C.LEVELS for k in range(4) for a in arms]
    with mp.Pool(par) as p:
        list(p.imap_unordered(R2.task, jobs))
    merge()
