#!/usr/bin/env python3
"""Stage 5 evaluation (EVALUATOR side) of the final B1/B2 benchmark, exactly as pre-specified in PROTOCOL_FINAL.md.
Outputs: evaluation/*.csv|json"""
import os, sys
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import json, multiprocessing as mp, warnings
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
FB = Path(__file__).resolve().parent; BB = FB.parent; CL = BB.parent
for p in (FB, CL / "v2", CL, CL / "homog_benchmark"):
    sys.path.insert(0, str(p))
import designs_common as D                                                          # noqa: E402
from bomes_v2 import GPv2, SCHEDULES, ANTIPODAL, axial_angle                         # noqa: E402
import analytic_reference as AR                                                     # noqa: E402  (region definition, Fibonacci sphere)
OUT = FB / "evaluation"; OUT.mkdir(exist_ok=True)
BUDGETS = (8, 12, 16, 20, 24, 28, 33); REG = ("uniaxial", "pure_shear", "tension_compression", "biaxial_same_sign", "equal_biaxial")
SHARE = pd.Series(AR.principal_region(AR.fibonacci_sphere(20000))[0]).value_counts(normalize=True).to_dict()
DELTA_NUM = {"B1": 0.010e-2, "B2": 0.034e-2}           # S1b-min load-step |dr/r| (fraction), from s1bmin/results/checks.json (rounded up to 2 sig. digits)
SUR = {"A20": SCHEDULES["A20"], "FIXED": SCHEDULES["FIXED"]}


def references(v):
    des = pd.read_csv(FB / "refs/reference_design.csv"); recs = {tuple(np.round(json.loads(l)["direction"], 12)): json.loads(l) for l in open(FB / "refs" / v / "trials.jsonl")}
    rows = []
    for r in des.itertuples():
        x = recs[tuple(np.round([r.z1, r.z2, r.z3], 12))]; R = json.load(open(BB / "runs" / v / "h40" / x["tag"] / "record.json"))
        rows.append(dict(ref_id=r.ref_id, set=r.set, region=r.region, z1=r.z1, z2=r.z2, z3=r.z3, r=x["observation"]["yield_radius"], accepted=x["observation"]["accepted"],
                         sweep_t_deg=getattr(r, "sweep_t_deg", np.nan), inclusion_share=R["volumes"]["f_i"] * R["at_yield"]["eqps_i_avg"] / R["at_yield"]["eqps_rve_avg"],
                         r1_matrix=R["first_local_yield_radius_probe"]["matrix"], r1_inclusion=R["first_local_yield_radius_probe"]["inclusion"],
                         first_plastic_matrix=R["first_plastic_ramp_step_radius"]["matrix"], first_plastic_inclusion=R["first_plastic_ramp_step_radius"]["inclusion"],
                         eqps_m=R["at_yield"]["eqps_m_avg"], eqps_i=R["at_yield"]["eqps_i_avg"], nl_its=R["nonlinear_iterations"], runtime_s=R["runtime_seconds"], flags=";".join(R["flags"])))
    return pd.DataFrame(rows)


def designs(v):
    out = {}
    for k in D.STARTS:
        t = pd.read_csv(FB / "campaigns" / v / f"bo_s{k}/tests.csv"); out[("BO", k)] = (t[["z1", "z2", "z3"]].to_numpy(), t.yield_radius.to_numpy(), t.axis.to_numpy())
    d = json.load(open(FB / "campaigns" / v / "nonadaptive/designs.json"))["designs"]
    for k in D.STARTS:
        x = d[f"axmm_s{k}"]; out[("AXMM", k)] = (np.array(x["directions"]), np.array(x["radii"]), np.array(x["axes"]))
    for k in D.RANDOM_SEEDS:
        x = d[f"random_{k}"]; out[("RANDOM", k)] = (np.array(x["directions"]), np.array(x["radii"]), np.array(x["axes"]))
    return out


def task(a):
    v, (meth, k), (X, r, ax), ref = a; P = ref[["z1", "z2", "z3"]].to_numpy(); rr = ref.r.to_numpy(); rows, preds = [], []
    for sur, sch in SUR.items():
        for n in BUDGETS:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore"); g = GPv2(ANTIPODAL, sch).fit(X[:n], r[:n])
            amp, ls, wn = g.theta(); mu, vl = g.predict_log(P); s = np.sqrt(vl + wn); e = np.exp(mu) / rr - 1; z = (np.log(rr) - mu) / s
            for st in ("area_stratified", "continuity_sweep"):
                m = (ref.set == st).to_numpy(); ae = np.abs(e[m]); reg = ref.region.to_numpy()[m]
                row = dict(variant=v, method=meth, design=k, surrogate=sur, n=n, ref_set=st, mean=ae.mean(), rms=np.sqrt(np.mean(e[m] ** 2)), p95=np.percentile(ae, 95), max=ae.max(),
                           cov90=np.mean(np.abs(z[m]) <= 1.645), cov95=np.mean(np.abs(z[m]) <= 1.96), rms_z=np.sqrt(np.mean(z[m] ** 2)), mean_signed=e[m].mean(),
                           length_scale=ls, amplitude=amp, white_noise=wn, gp_mode=g.mode_)
                if st == "area_stratified":
                    rm = {rg: ae[reg == rg].mean() for rg in REG}; row.update({f"mean_{rg}": rm[rg] for rg in REG})
                    row["area_weighted"] = sum(rm[rg] * SHARE[rg] for rg in REG) / sum(SHARE[rg] for rg in REG)
                rows.append(row)
            if n == 33 or n == 20 or n == 12:
                preds += [dict(variant=v, method=meth, design=k, surrogate=sur, n=n, ref_id=i, pred_r=float(np.exp(m_)), sd_log=float(s_), meas_r=float(r_), z=float(z_))
                          for i, m_, s_, r_, z_ in zip(ref.ref_id, mu, s, rr, z)]
    return rows, preds


def first_div(a, b):
    return next((i + 1 for i, (x, y) in enumerate(zip(a, b)) if x != y), None)


def main():
    refs = {v: references(v) for v in D.VARIANTS}; des = {v: designs(v) for v in D.VARIANTS}
    for v in D.VARIANTS:
        refs[v].to_csv(OUT / f"references_{v}.csv", index=False)
    jobs = [(v, key, val, refs[v]) for v in D.VARIANTS for key, val in des[v].items()]
    with mp.Pool(22) as p:
        res = p.map(task, jobs, chunksize=1)
    M = pd.DataFrame([x for r_, _ in res for x in r_]); M.to_csv(OUT / "metrics_all.csv", index=False)
    pd.DataFrame([x for _, p_ in res for x in p_]).to_csv(OUT / "predictions_at_references.csv", index=False)
    A = M[M.ref_set == "area_stratified"]
    summ = A.groupby(["variant", "surrogate", "method", "n"])[["area_weighted", "mean", "rms", "p95", "max", "cov90", "cov95", "rms_z"] + [f"mean_{r}" for r in REG]].agg(["mean", "median", "min", "max"])
    summ.columns = ["_".join(c) for c in summ.columns]; summ.reset_index().to_csv(OUT / "summary.csv", index=False)
    # paired BO vs AXMM; BO vs RANDOM
    pr, claims = [], {}
    for v in D.VARIANTS:
        for sur in SUR:
            for n in BUDGETS:
                for met in ("area_weighted", "mean", "rms", "p95", "max"):
                    g = A[(A.variant == v) & (A.surrogate == sur) & (A.n == n)]
                    bo = g[g.method == "BO"].set_index("design")[met].reindex(D.STARTS).to_numpy(); am = g[g.method == "AXMM"].set_index("design")[met].reindex(D.STARTS).to_numpy()
                    ra = g[g.method == "RANDOM"][met].to_numpy(); d = am - bo
                    w = stats.wilcoxon(d).pvalue if np.any(d != 0) else np.nan
                    pr.append(dict(variant=v, surrogate=sur, n=n, metric=met, bo_mean=bo.mean(), axmm_mean=am.mean(), random_mean=ra.mean(), bo_better_in=int((d > 0).sum()),
                                   axmm_better_in=int((d < 0).sum()), ties=int((d == 0).sum()), paired_median_improvement=float(np.median(d)), paired_mean_improvement=float(d.mean()),
                                   paired_iqr=float(np.subtract(*np.percentile(d, [75, 25]))), wilcoxon_p=w, mw_p_bo_vs_random=float(stats.mannwhitneyu(bo, ra).pvalue),
                                   bo_better_than_random_median=bool(np.median(bo) < np.median(ra))))
        P_ = pd.DataFrame(pr); q = lambda sur, met: P_[(P_.variant == v) & (P_.surrogate == sur) & (P_.n == 33) & (P_.metric == met)].iloc[0]
        c1 = q("A20", "area_weighted").bo_better_in >= 8; c2 = q("A20", "area_weighted").paired_median_improvement > DELTA_NUM[v]
        c3 = q("FIXED", "area_weighted").bo_better_in >= 8 and q("FIXED", "area_weighted").paired_median_improvement > DELTA_NUM[v]; c4 = q("A20", "mean").paired_median_improvement > 0
        claims[v] = dict(c1_8_of_9=bool(c1), c2_exceeds_delta_num=bool(c2), c3_robust_FIXED=bool(c3), c4_mean_supports=bool(c4), BO_ADVANTAGE_CLAIMED=bool(c1 and c2 and c3 and c4),
                         delta_num=DELTA_NUM[v], primary=q("A20", "area_weighted").to_dict(), fixed=q("FIXED", "area_weighted").to_dict())
    pd.DataFrame(pr).to_csv(OUT / "paired_comparison.csv", index=False); json.dump(claims, open(OUT / "claims.json", "w"), indent=1, default=float)
    # selections, response dependence, hyper-parameters, acquisition
    sel = {}
    for k in D.STARTS:
        s = {v: list(des[v][("BO", k)][2]) for v in D.VARIANTS}; am = list(des["B1"][("AXMM", k)][2])
        hb = pd.read_csv(CL / f"homog_benchmark/results/campaigns/HB_V2_A20_start{k}/tests.csv").axis.tolist()
        sel[k] = dict(B1_vs_B2_first_div=first_div(s["B1"], s["B2"]), B1_vs_B2_shared=len(set(s["B1"]) & set(s["B2"])), B1_vs_homog_first_div=first_div(s["B1"], hb),
                      B2_vs_homog_first_div=first_div(s["B2"], hb), B1_vs_axmm_first_div=first_div(s["B1"], am), B2_vs_axmm_first_div=first_div(s["B2"], am),
                      B1_vs_axmm_shared=len(set(s["B1"]) & set(am)), B2_vs_axmm_shared=len(set(s["B2"]) & set(am)))
    sel["region_of_adaptive_picks"] = {v: pd.Series(np.concatenate([AR.principal_region(des[v][("BO", k)][0][8:])[0] for k in D.STARTS])).value_counts().to_dict() for v in D.VARIANTS}
    sel["region_of_axmm_picks"] = pd.Series(np.concatenate([AR.principal_region(des["B1"][("AXMM", k)][0][8:])[0] for k in D.STARTS])).value_counts().to_dict()
    sel["region_share"] = SHARE; json.dump(sel, open(OUT / "selection_analysis.json", "w"), indent=1, default=int)
    st = []
    for v in D.VARIANTS:
        for k in D.STARTS:
            s = pd.read_csv(FB / "campaigns" / v / f"bo_s{k}/steps.csv"); sc = np.load(FB / "campaigns" / v / f"bo_s{k}/acquisition_scores.npy")
            fin = [np.sort(x[np.isfinite(x)])[::-1] for x in sc]; s["variant"], s["start"] = v, k
            s["runner_up_rel_gap"] = [(f[0] - f[1]) / f[0] for f in fin]; s["score_median"] = [float(np.median(f)) for f in fin]; st.append(s)
    st = pd.concat(st); st.to_csv(OUT / "bo_steps_all.csv", index=False)
    # cost, failures, repeatability of concurrent duplicates
    cost = {}
    for v in D.VARIANTS:
        logs = {"references": FB / "refs" / v / "trials.jsonl", "nonadaptive": FB / "campaigns" / v / "nonadaptive/trials.jsonl"} | {f"bo_s{k}": FB / "campaigns" / v / f"bo_s{k}/trials.jsonl" for k in D.STARTS}
        recs = {n: [json.loads(l) for l in open(p)] for n, p in logs.items()}; new = [r for x in recs.values() for r in x if "reused_from" not in r]
        dup = {}
        for r in new:
            dup.setdefault(tuple(np.round(r["direction"], 12)), set()).add(r["observation"].get("yield_radius"))
        cost[v] = dict(new_trials={n: sum("reused_from" not in r for r in x) for n, x in recs.items()}, reused={n: sum("reused_from" in r for r in x) for n, x in recs.items()},
                       failed={n: sum(not r["observation"].get("accepted") for r in x) for n, x in recs.items()}, total_new=len(new),
                       distinct_new_directions=len(dup), duplicate_runs_bitwise_identical=all(len(s_) == 1 for s_ in dup.values()), n_duplicate_directions=sum(1 for r in new) - len(dup),
                       cpu_h=sum(r["worker_wall_s"] for r in new) / 3600, mean_trial_s=float(np.mean([r["worker_wall_s"] for r in new])))
    json.dump(cost, open(OUT / "cost.json", "w"), indent=1, default=float)
    pd.set_option("display.width", 250)
    print((summ[["area_weighted_mean", "mean_mean", "p95_mean", "max_mean"]].unstack("method") * 100).round(3).to_string())
    print(json.dumps(claims, indent=1, default=float)[:3000]); print(json.dumps(cost, indent=1, default=float))


if __name__ == "__main__":
    main()
