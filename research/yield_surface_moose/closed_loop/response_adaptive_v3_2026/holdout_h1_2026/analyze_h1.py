#!/usr/bin/env python3
"""Pre-registered analysis of holdout H1 (PROTOCOL.md section 5; frozen before any controller ran).
Paired statistic per realisation: L = log(err_other / err_V3)  (> 0: V3 better).  Statistics (identical code to the frozen LSF-v2 analysis):
exact two-sided Wilcoxon signed-rank, Hodges-Lehmann ratio exp(HL) with its exact 95 % CI, wins.
CONFIRMATORY FAMILY (Holm, alpha = 0.05):
  H1  V3 vs AXMM,  glob_mean, n = 33   (PRIMARY)
  H2  V3 vs V2,    glob_mean, n = 33
  H3  V3 vs AXMM,  glob_mean, n = 64
  Hk CONFIRMED (practical superiority) iff Holm-adjusted p < 0.05 AND exp(HL) >= 1.10 AND lower 95 % CI > 1 AND wins >= 16/24.
  Hk INFERIOR iff upper 95 % CI < 1.  Otherwise NOT ESTABLISHED.
TRADE-OFF GUARDS (V3 vs AXMM, n = 33): high-curvature-region mean error and global max error, lower 95 % CI of exp(HL) >= 0.90 each;
  failing either = 'unacceptable trade-off'.
DESCRIPTIVE (nominal p, no claim): all metrics x {33, 64} x {AXMM, V2, RANDOM, V3_ISO}; response attribution V3 vs V3_ISO design;
  moderators (Spearman of L_H1 with anisotropy max/min and kappa_1 r p99; per-level medians); error vs n; coverage; selection
  differences; cost.
    python analyze_h1.py  ->  analysis/"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
HERE = Path(__file__).resolve().parent; V3D = HERE.parent; sys.dont_write_bytecode = True; sys.path.insert(0, str(V3D)); sys.path.insert(0, str(HERE))
from evaluate_v3 import paired                                                       # noqa: E402  (== frozen analyze_stress.paired; self-checked in the V3 pilot)
AN = HERE / "analysis"; AN.mkdir(exist_ok=True); M = pd.read_csv(HERE / "evaluation/metrics_all.csv"); NJ = 24
POP = json.load(open(HERE / "population.json"))["realisations"]; VER = {v["realisation"]: v for v in json.load(open(HERE / "surface_verification.json"))["checks"] if v["ok"]}


def L(n, metric, other):
    P = M[M.n == n].pivot_table(index="realisation", columns="arm", values=metric); return np.log(P[other] / P["V3"]), P


def holm(ps):
    o = np.argsort(ps); adj = np.empty(len(ps)); run = 0.0
    for r, i in enumerate(o):
        run = max(run, min(1.0, (len(ps) - r) * ps[i])); adj[i] = run
    return adj


fam = [("H1", "AXMM", "glob_mean", 33), ("H2", "V2", "glob_mean", 33), ("H3", "AXMM", "glob_mean", 64)]; res = {}
for h, o, met, n in fam:
    res[h] = dict(comparison=f"V3 vs {o}", metric=met, n=n, **paired(L(n, met, o)[0]))
adj = holm(np.array([res[h]["wilcoxon_p"] for h, *_ in fam]))
for (h, *_), a in zip(fam, adj):
    r = res[h]; r["holm_p"] = float(a)
    r["verdict"] = ("CONFIRMED" if (a < 0.05 and r["hl_ratio"] >= 1.10 and r["hl_ci95"][0] > 1 and r["wins"] >= 16) else "INFERIOR" if r["hl_ci95"][1] < 1 else "NOT ESTABLISHED")
guards = {}
for met in ("hc_mean", "glob_max"):
    g = paired(L(33, met, "AXMM")[0]); guards[met] = dict(**g, met=bool(g["hl_ci95"][0] >= 0.90))
rows = []
for n in (33, 64):
    for met in ("glob_mean", "glob_max", "hc_mean", "hc_max", "out_mean"):
        for o in ("AXMM", "V2", "RANDOM", "V3_ISO"):
            x, Q = L(n, met, o); s = paired(x)
            rows.append(dict(n=n, metric=met, vs=o, **{k: v for k, v in s.items() if k != "hl_ci95"}, hl_ci95_lo=s["hl_ci95"][0], hl_ci95_hi=s["hl_ci95"][1],
                             median_V3_pct=100 * Q.V3.median(), median_other_pct=100 * Q[o].median()))
CMP = pd.DataFrame(rows); CMP.to_csv(AN / "comparisons.csv", index=False)
L1 = L(33, "glob_mean", "AXMM")[0]; desc = pd.DataFrame([dict(realisation=i, level=POP[i]["level"], m=POP[i]["m"], s=POP[i]["s"], aniso=VER[i]["anisotropy_max_over_min"],
                                                              kappa_p99=VER[i]["kappa1_dimless_p99"], L_H1=float(L1.loc[i]), L_vs_V2=float(L(33, "glob_mean", "V2")[0].loc[i]),
                                                              L_vs_ISO=float(L(33, "glob_mean", "V3_ISO")[0].loc[i])) for i in range(NJ)])
desc.to_csv(AN / "per_realisation.csv", index=False)
mod = {k: dict(rho=float(stats.spearmanr(desc[k], desc.L_H1).statistic), p=float(stats.spearmanr(desc[k], desc.L_H1).pvalue)) for k in ("aniso", "kappa_p99", "m")}
lev = desc.groupby("level").agg(median_ratio_AXMM_over_V3=("L_H1", lambda x: float(np.exp(np.median(x)))), wins=("L_H1", lambda x: int((x > 0).sum())),
                                median_ratio_V2_over_V3=("L_vs_V2", lambda x: float(np.exp(np.median(x))))).reset_index()
lev.to_csv(AN / "by_level.csv", index=False)
LIB_ = None; sel = []
for i in range(NJ):
    a = json.load(open(HERE / f"runs/V3/r{i:02d}/design.json"))["axes"]; b = json.load(open(HERE / f"runs/V3_ISO/r{i:02d}/design.json"))["axes"]; c = json.load(open(HERE / f"runs/V2/r{i:02d}/design.json"))["axes"]
    st = pd.read_csv(HERE / f"runs/V3/r{i:02d}/steps.csv")
    sel.append(dict(realisation=i, D33_V3_vs_ISO=len(set(a[8:33]) - set(b[8:33])), D33_V3_vs_V2=len(set(a[8:33]) - set(c[8:33])), first_div_V3_vs_ISO=next((k + 1 for k, (x, y) in enumerate(zip(a, b)) if x != y), None),
                    share_differs_from_MES_9_33=float(st[st.test_number <= 33].differs_from_mes.mean()), fallback=int(st.fallback.sum())))
SEL = pd.DataFrame(sel); SEL.to_csv(AN / "selection_differences.csv", index=False)
cov = M[M.n.isin([33, 64])].groupby(["n", "arm"])[["glob_cov90", "glob_cov95", "hc_cov90", "hc_cov95", "glob_rms_z"]].median().reset_index(); cov.to_csv(AN / "coverage.csv", index=False)
curve = M.groupby(["arm", "n"])[["glob_mean", "hc_mean", "glob_max"]].median().reset_index(); curve.to_csv(AN / "error_vs_n_median.csv", index=False)
wall = {a: [json.load(open(HERE / f"runs/{a}/r{i:02d}/run_info.json"))["wall_s"] for i in range(NJ)] for a in ("V3", "V2", "AXMM", "RANDOM", "V3_ISO")}
cost = dict(wall_s_median={a: float(np.median(v)) for a, v in wall.items()}, wall_s_total={a: float(np.sum(v)) for a, v in wall.items()},
            V3_over_V2_median=float(np.median(np.array(wall["V3"]) / np.array(wall["V2"]))))
attr = paired(L(33, "glob_mean", "V3_ISO")[0])
out = dict(status="independent holdout (Yld91 family, fresh seeds/starts/references; frozen before any controller run)", confirmatory=res,
           primary_verdict=res["H1"]["verdict"], tradeoff_guards=guards, unacceptable_tradeoff=not all(g["met"] for g in guards.values()),
           response_attribution_V3_vs_V3ISO_glob33=attr, moderators_spearman_L_H1=mod, by_level=lev.to_dict("records"),
           selection=dict(D33_V3_vs_ISO_median=float(SEL.D33_V3_vs_ISO.median()), D33_V3_vs_ISO_range=[int(SEL.D33_V3_vs_ISO.min()), int(SEL.D33_V3_vs_ISO.max())],
                          D33_V3_vs_V2_median=float(SEL.D33_V3_vs_V2.median()), share_differs_from_MES_median=float(SEL.share_differs_from_MES_9_33.median()), fallback_steps=int(SEL.fallback.sum())),
           cost=cost)
json.dump(out, open(AN / "results.json", "w"), indent=1, default=float); print(json.dumps(out, indent=1, default=float))
