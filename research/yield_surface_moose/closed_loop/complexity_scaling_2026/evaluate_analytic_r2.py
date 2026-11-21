#!/usr/bin/env python3
"""Versioned evaluator for the post-reboot closeout (r2).  Same metrics, references, surrogates and code path as evaluate_analytic.py
(its refset() is imported unchanged), restricted to the protocol grid n <= 512: the BO campaigns were stopped at 512 recorded observations
(stop_at_512.py; resumable to 1000).  BO designs = original evaluations.jsonl + campaigns_r2 continuation (prefix verified in
logs/resume_selftest.json).  Resumable: one CSV per design in evaluation_r2/designs/.     python evaluate_analytic_r2.py [PAR] [ARMS...]"""
import os, sys
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "1")
import json, time, warnings, multiprocessing as mp
from pathlib import Path
import numpy as np, pandas as pd
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import cs_common as C                                                              # noqa: E402
import evaluate_analytic as E                                                      # noqa: E402  (refset, GRID)
from bomes_v2 import GPv2, SCHEDULES, ANTIPODAL                                     # noqa: E402
NMAX = 512; GRID = tuple(n for n in E.GRID if n <= NMAX); BUD = tuple(n for n in C.BUDGETS if n <= NMAX)
OUT = HERE / "evaluation_r2" / "designs"; OUT.mkdir(parents=True, exist_ok=True)


def design(level, arm, start):
    lib = C.library()
    if arm == "BO":
        recs = [json.loads(l) for p in (HERE / "campaigns" / level / f"BO_s{start}" / "evaluations.jsonl", HERE / "campaigns_r2" / level / f"BO_s{start}" / "evaluations_cont.jsonl")
                if p.exists() for l in open(p)]
        if len(recs) < NMAX:
            return None
        D = np.array([r["direction"] for r in recs[:NMAX]]); idx = [int(np.argmin(np.abs(lib.directions - d).sum(1))) for d in D]
        assert np.array_equal(lib.directions[idx], D); return lib.directions[idx], np.array([r["yield_radius"] for r in recs[:NMAX]])
    src = "ANY" if arm == "BOFIXED" else level; d = json.load(open(HERE / "campaigns" / src / f"{arm}_s{start}" / "design.json")); X = lib.directions[d["axes"][:NMAX]]
    return X, (np.asarray(d["radii"][:NMAX]) if arm != "BOFIXED" else C.radius_fn(level)(X))


def task(a):
    level, arm, start = a; out = OUT / f"{level}__{arm}_s{start}.csv"
    if out.exists():
        return
    dz = design(level, arm, start)
    if dz is None:
        return
    X, r = dz; D, rr, hi = E.refset(level); rows = []
    for sur in ("A20", "FIXED"):
        for n in (GRID if sur == "A20" else BUD):
            t = time.perf_counter()
            with warnings.catch_warnings():
                warnings.simplefilter("ignore"); g = GPv2(ANTIPODAL, SCHEDULES[sur]).fit(X[:n], r[:n])
            amp, ls, wn = g.theta(); mu, vl = g.predict_log(D); s = np.sqrt(vl + wn); e = np.exp(mu) / rr - 1; ae = np.abs(e); z = (np.log(rr) - mu) / s; hp = g.hyperparameters()
            rows.append(dict(level=level, arm=arm, start=start, surrogate=sur, n=n, mean=ae.mean(), rms=np.sqrt(np.mean(e ** 2)), p95=np.percentile(ae, 95), max=ae.max(),
                             mean_highcurv=ae[hi].mean(), max_highcurv=ae[hi].max(), mean_lowcurv=ae[~hi].mean(), cov90=np.mean(np.abs(z) <= 1.645), rms_z=np.sqrt(np.mean(z ** 2)),
                             mean_sd=s.mean(), mean_sd_highcurv=s[hi].mean(), corr_sd_abs_err=float(np.corrcoef(s, ae)[0, 1]), amplitude=amp, length_scale=ls, white_noise=wn,
                             at_bounds=";".join(hp["at_bounds"]), gp_mode=g.mode_, fit_s=time.perf_counter() - t))
    pd.DataFrame(rows).to_csv(out, index=False); print("evaluated", out.name, time.strftime("%H:%M:%S"), flush=True)


if __name__ == "__main__":
    par = int(sys.argv[1]) if len(sys.argv) > 1 else 4; arms = sys.argv[2:] or ["AXMM", "RANDOM", "BOFIXED", "BO"]
    jobs = [(L, a, k) for L in C.LEVELS for k in range(4) for a in arms]
    with mp.Pool(par) as p:
        list(p.imap_unordered(task, jobs))
    fs = sorted(OUT.glob("*.csv")); M = pd.concat([pd.read_csv(f) for f in fs]); M.to_csv(HERE / "evaluation_r2" / "metrics_all.csv", index=False)
    print(len(M), "rows from", len(fs), "designs", flush=True)
