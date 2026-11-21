#!/usr/bin/env python3
"""Evaluator (pre-specified in PROTOCOL.md) for the analytic complexity-scaling campaigns.  Resumable: one CSV per design in evaluation/designs/.
For every (level, arm, start) and every n on the grid: frozen GPv2(ANTIPODAL, A20) refit on the first n observations -> 20,000 held-out
Fibonacci references with EXACT r.  FIXED-schedule surrogate as sensitivity at the 8 budgets.  BOFIXED designs are evaluated on every level.
    python evaluate_analytic.py [PAR]"""
import os, sys
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "1")
import json, time, warnings, multiprocessing as mp
from pathlib import Path
import numpy as np, pandas as pd
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import cs_common as C                                                              # noqa: E402
from bomes_v2 import GPv2, SCHEDULES, ANTIPODAL                                     # noqa: E402
GRID = (8, 12, 16, 20, 24, 28, 33, 40, 48, 64, 80, 96, 128, 160, 192, 256, 320, 384, 512, 640, 768, 896, 1000)
OUT = HERE / "evaluation" / "designs"; OUT.mkdir(parents=True, exist_ok=True)
_REF = {}


def refset(level):
    if level not in _REF:
        D = C.references(); f = C.radius_fn(level); r = f(D); k1, _ = C.curvatures(f, D); kn = k1 * r
        hi = kn >= np.percentile(kn, 90); _REF[level] = (D, r, hi)
    return _REF[level]


def task(a):
    level, arm, start = a; src = "ANY" if arm == "BOFIXED" else level; out = OUT / f"{level}__{arm}_s{start}.csv"
    if out.exists():
        return str(out)
    d = json.load(open(HERE / "campaigns" / src / f"{arm}_s{start}" / "design.json")); lib = C.library(); X = lib.directions[d["axes"]]
    r = np.asarray(d["radii"]) if arm != "BOFIXED" else C.radius_fn(level)(X)          # BOFIXED design is response-independent: re-observe on this level
    D, rr, hi = refset(level); rows = []
    for sur in ("A20", "FIXED"):
        for n in (GRID if sur == "A20" else C.BUDGETS):
            t = time.perf_counter()
            with warnings.catch_warnings():
                warnings.simplefilter("ignore"); g = GPv2(ANTIPODAL, SCHEDULES[sur]).fit(X[:n], r[:n])
            amp, ls, wn = g.theta(); mu, vl = g.predict_log(D); s = np.sqrt(vl + wn); e = np.exp(mu) / rr - 1; ae = np.abs(e); z = (np.log(rr) - mu) / s
            hp = g.hyperparameters()
            rows.append(dict(level=level, arm=arm, start=start, surrogate=sur, n=n, mean=ae.mean(), rms=np.sqrt(np.mean(e ** 2)), p95=np.percentile(ae, 95), max=ae.max(),
                             mean_highcurv=ae[hi].mean(), max_highcurv=ae[hi].max(), mean_lowcurv=ae[~hi].mean(), cov90=np.mean(np.abs(z) <= 1.645), rms_z=np.sqrt(np.mean(z ** 2)),
                             mean_sd=s.mean(), mean_sd_highcurv=s[hi].mean(), corr_sd_abs_err=float(np.corrcoef(s, ae)[0, 1]), amplitude=amp, length_scale=ls, white_noise=wn,
                             at_bounds=";".join(hp["at_bounds"]), gp_mode=g.mode_, fit_s=time.perf_counter() - t))
    pd.DataFrame(rows).to_csv(out, index=False); print("evaluated", out.name, flush=True); return str(out)


if __name__ == "__main__":
    par = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    jobs = []
    for level in C.LEVELS:
        for k in range(4):
            for arm in ("AXMM", "RANDOM", "BOFIXED", "BO"):
                src = "ANY" if arm == "BOFIXED" else level
                if (HERE / "campaigns" / src / f"{arm}_s{k}" / "run_info.json").exists():
                    jobs.append((level, arm, k))
    with mp.Pool(par) as p:
        list(p.imap_unordered(task, jobs))
    M = pd.concat([pd.read_csv(f) for f in sorted(OUT.glob("*.csv"))]); M.to_csv(HERE / "evaluation" / "metrics_all.csv", index=False)
    print(len(M), "rows from", len(list(OUT.glob("*.csv"))), "designs")
