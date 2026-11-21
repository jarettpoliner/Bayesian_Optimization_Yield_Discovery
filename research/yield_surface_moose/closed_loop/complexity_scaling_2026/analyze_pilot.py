#!/usr/bin/env python3
"""Physical pilot analysis: measured complexity of P2-P6 (+ P1 analytic J2 / HOMOG2), mesh check, cost and full-campaign resource estimate.
Complexity class from the frozen GPv2 fitted to the 24 pilot observations (coarse: 24 axes resolve features of ~20 deg only; curvature is a
LOWER bound where the surface has features narrower than the axis spacing).  -> pilot/pilot_summary.json, pilot/pilot_*.png"""
import json, sys, warnings
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import cs_common as C                                                              # noqa: E402
from bomes_v2 import GPv2, SCHEDULES, ANTIPODAL                                     # noqa: E402
CONFIGS = ["P2_single", "P3_regular", "P4_clustered_s4401", "P5_multiscale_s5501", "P6_aligned_hiAR_s6601"]


def load(cfg, nel=40):
    f = HERE / "pilot" / cfg / "trials.jsonl"
    R = [json.loads(l) for l in open(f)] if f.exists() else []
    return [r for r in R if r["nel"] == nel]


def main():
    out, rows = {}, []
    J2 = C.descriptors(C.radius_fn("A1_J2")); J2.pop("corr_curve")
    out["P1_homog"] = dict(source="analytic J2 (HOMOG2 MOOSE oracle validated against it in the homog benchmark)", **{k: J2[k] for k in ("anisotropy_max_over_min", "kappa_max_dimless_p99", "curvature_top10pct_share")},
                           cls=C.complexity_class(J2["kappa_max_dimless_p99"]))
    for cfg in CONFIGS:
        R = load(cfg); ok = [r for r in R if r["observation"].get("accepted")]
        if len(ok) < 8:
            out[cfg] = dict(n_runs=len(R), n_accepted=len(ok), note="insufficient pilot data"); continue
        D = np.array([r["direction"] for r in ok]); r = np.array([r["observation"]["yield_radius"] for r in ok])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore"); g = GPv2(ANTIPODAL, SCHEDULES["A20"]).fit(D, r)
        fn = lambda Z: np.exp(g.predict_log(np.atleast_2d(Z))[0]); d = C.descriptors(fn, n=6000); d.pop("corr_curve")
        wall = np.array([x["wall_s"] for x in R]); rt = np.array([x["observation"].get("runtime_seconds") or np.nan for x in R], float)
        out[cfg] = dict(n_runs=len(R), n_accepted=len(ok), r_min_measured=float(r.min()), r_max_measured=float(r.max()), anisotropy_measured=float(r.max() / r.min()),
                        std_log_r_measured=float(np.log(r).std()), gp_hyper=[float(x) for x in g.theta()], gp_mode=g.mode_,
                        descriptors_from_gp_fit=d, cls=C.complexity_class(d["kappa_max_dimless_p99"]), wall_s_median=float(np.median(wall)), runtime_s_median=float(np.nanmedian(rt)),
                        cpu_h=float(np.nansum(rt) / 3600))
        rows.append(dict(config=cfg, **{k: out[cfg][k] for k in ("n_accepted", "anisotropy_measured", "std_log_r_measured", "cls", "runtime_s_median")},
                         kappa_p99=d["kappa_max_dimless_p99"], top10=d["curvature_top10pct_share"], convex=d["convex"]))
    # mesh check P6 h80 vs h40 (same axes)
    h80 = {x["axis_index"]: x for x in load("P6_aligned_hiAR_s6601", 80)}; h40 = {x["axis_index"]: x for x in load("P6_aligned_hiAR_s6601", 40)}
    out["mesh_check_P6"] = [dict(axis=i, r_h40=h40[i]["observation"]["yield_radius"], r_h80=x["observation"]["yield_radius"],
                                 rel_diff=(x["observation"]["yield_radius"] / h40[i]["observation"]["yield_radius"] - 1) if (i in h40 and x["observation"].get("yield_radius") and h40[i]["observation"].get("yield_radius")) else None,
                                 runtime_h80_s=x["observation"].get("runtime_seconds")) for i, x in h80.items()]
    # full-campaign resource estimate (PROTOCOL arms, 1,000 observations, 4 starts, 5 configs x 2 geometry seeds where predefined)
    rt = np.nanmedian([out[c]["runtime_s_median"] for c in CONFIGS if "runtime_s_median" in out[c]])
    n_cfg = 5 + 3                                                           # P2, P3, P4a, P5a, P6a + second seeds of P4, P5, P6
    per_cfg_runs = 4 * 1000 + 4 * 992 + 4 * 992 + 200                       # BO 4 starts; AXMM / RANDOM (shared first 8 axes); 200 references
    reuse_est = 0.0                                                         # exact (config, direction) reuse between arms is rare on a 4,096-axis library
    cpu_h = n_cfg * per_cfg_runs * rt / 3600
    bo_wall_h = 1000 * rt / 3600 + 2.5                                       # BO is sequential: 1,000 MOOSE solves one after another + ~2.5 h of GP refits
    out["full_campaign_estimate"] = dict(median_trial_runtime_s=float(rt), configs=n_cfg, runs_per_config=per_cfg_runs, total_runs=n_cfg * per_cfg_runs,
                                         total_cpu_h=float(cpu_h), wall_h_at_12way=float(cpu_h / 12), bo_campaign_wall_h_each_sequential=float(bo_wall_h),
                                         storage_gb_est=float(n_cfg * per_cfg_runs * 0.9e-3), reference_cpu_h_only=float(n_cfg * 200 * rt / 3600),
                                         approved_reference_budget_cpu_h=23.2, fits_in_approved_budget=bool(cpu_h <= 23.2))
    json.dump(out, open(HERE / "pilot/pilot_summary.json", "w"), indent=1, default=float)
    T = pd.DataFrame(rows); T.to_csv(HERE / "pilot/pilot_complexity.csv", index=False); print(T.round(3).to_string()); print(json.dumps(out["full_campaign_estimate"], indent=1)); print(out["mesh_check_P6"])


if __name__ == "__main__":
    main()
