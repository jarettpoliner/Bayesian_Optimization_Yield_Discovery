#!/usr/bin/env python3
"""H1 evaluator.  Same surrogate protocol as every earlier study: the frozen GPv2(ANTIPODAL, A20) is refit from scratch on each design's first n
observations (ACTUAL Yld91 radii at the design's axes -- for V3_ISO this re-observes the ablation design) and predicts the 20,000 rotated
equal-area H1 references (references_h1.npy).  Equal-area => the plain mean is the area-weighted mean.
High-curvature region HC: references with kappa_1 r in the top 10 % for that surface (cs_common.curvatures; the complexity-study definition).
`selfcheck` first proves that this generic metric code reproduces the frozen LSF-v2 metrics (glob/feat mean & max, cov90) of BO/AXMM r00.
    nice -n 19 python evaluate_h1.py selfcheck | all [workers]   ->  evaluation/{selfcheck.json, metrics_all.csv}"""
import os, sys
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "1")
import json, time, warnings
from pathlib import Path
import numpy as np, pandas as pd
HERE = Path(__file__).resolve().parent; V3D = HERE.parent; CL = V3D.parent; CS = CL / "complexity_scaling_2026"; LSF = CL / "adaptivity_stress_test_2026"; sys.dont_write_bytecode = True
for p in (str(HERE), str(CS), str(LSF)):
    sys.path.insert(0, p)
import cs_common as C                                                                # noqa: E402
from bomes_v2 import GPv2, SCHEDULES, ANTIPODAL                                      # noqa: E402  frozen
import surfaces_h1 as H                                                              # noqa: E402
GRID = (8, 12, 16, 20, 24, 28, 33, 40, 48, 56, 64); EV = HERE / "evaluation"; EV.mkdir(exist_ok=True)


def metrics(X, r, R, rr, hc, grid, sch="A20"):
    rows = []
    for n in grid:
        t = time.perf_counter()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore"); g = GPv2(ANTIPODAL, SCHEDULES[sch]).fit(X[:n], r[:n])
        amp, ls, wn = g.theta(); mu, vl = g.predict_log(R); sd = np.sqrt(vl + wn); ae = np.abs(np.exp(mu) / rr - 1); z = (np.log(rr) - mu) / sd
        rows.append(dict(n=n, glob_mean=ae.mean(), glob_max=ae.max(), hc_mean=ae[hc].mean(), hc_max=ae[hc].max(), out_mean=ae[~hc].mean(),
                         glob_cov90=np.mean(np.abs(z) <= 1.645), glob_cov95=np.mean(np.abs(z) <= 1.96), hc_cov90=np.mean(np.abs(z[hc]) <= 1.645), hc_cov95=np.mean(np.abs(z[hc]) <= 1.96),
                         glob_rms_z=np.sqrt(np.mean(z ** 2)), length_scale=ls, amplitude=amp, white_noise=wn, gp_mode=g.mode_, fit_s=time.perf_counter() - t))
    return rows


def selfcheck():
    import surface_lsf as S
    M = pd.read_csv(LSF / "evaluation/metrics_all.csv"); out = {}; LIB = C.library().directions; R = C.references()
    for arm in ("BO", "AXMM"):
        k = json.load(open(LSF / f"runs/{arm}/r00/run_info.json"))["orientation_k"]; F = S.LSF(k); d = json.load(open(LSF / f"runs/{arm}/r00/design.json"))
        a = pd.DataFrame(metrics(LIB[d["axes"]], np.asarray(d["radii"]), R, F.radius(R), F.feature_region(R), (8, 20, 33, 64)))
        b = M[(M.arm == arm) & (M.realisation == 0) & M.n.isin([8, 20, 33, 64])].sort_values("n").reset_index(drop=True)
        out[arm] = float(max(np.abs(a.glob_mean - b.glob_mean).max(), np.abs(a.hc_mean - b.feat_mean).max(), np.abs(a.glob_max - b.glob_max).max(),
                             np.abs(a.hc_max - b.feat_max).max(), np.abs(a.glob_cov90 - b.glob_cov90).max(), np.abs(a.hc_cov90 - b.feat_cov90).max()))
    out["PASS"] = bool(max(out["BO"], out["AXMM"]) <= 1e-14); json.dump(out, open(EV / "selfcheck.json", "w"), indent=1); print(out)
    if not out["PASS"]:
        sys.exit("selfcheck FAILED")


def task(a):
    arm, i = a; rec = H.POP[i] if hasattr(H, "POP") else json.load(open(HERE / "population.json"))["realisations"][i]; F = H.Yld91(rec)
    R = np.load(HERE / "references_h1.npy"); rr = F.radius(R); k1, _ = C.curvatures(F.radius, R); kn = k1 * rr; hc = kn >= np.percentile(kn, 90)
    d = json.load(open(HERE / "runs" / arm / f"r{i:02d}" / "design.json")); X = C.library().directions[d["axes"]]; r = F.radius(X)
    if arm != "V3_ISO":
        assert np.array_equal(r, np.asarray(d["radii"])) or np.allclose(r, d["radii"], rtol=1e-15, atol=0), "recorded radii must equal the exact oracle"
    return [dict(arm=arm, realisation=i, level=rec["level"], m=rec["m"], s=rec["s"], **x) for x in metrics(X, r, R, rr, hc, GRID)]


if __name__ == "__main__":
    if sys.argv[1] == "selfcheck":
        selfcheck()
    else:
        import multiprocessing as mp
        jobs = [(a, i) for a in ("V3", "V2", "AXMM", "RANDOM", "V3_ISO") for i in range(24)]
        assert all((HERE / "runs" / a / f"r{i:02d}" / "run_info.json").exists() for a, i in jobs), "all 120 campaigns must be complete"
        with mp.Pool(int(sys.argv[2]) if len(sys.argv) > 2 else 4) as P:
            rows = [r for rr in P.map(task, jobs, chunksize=1) for r in rr]
        pd.DataFrame(rows).to_csv(EV / "metrics_all.csv", index=False); print(len(rows), "rows")
