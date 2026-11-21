#!/usr/bin/env python3
"""Pre-specified analysis (PROTOCOL.md) of the analytic complexity-scaling study, r2 closeout: evaluation_r2/metrics_all.csv (grid n <= 512;
BO stopped at 512 after the reboot, resumable to 1000).  Statistics exactly as PROTOCOL.md:
  per level and budget: BO vs AXMM / RANDOM over 4 matched starts (wins, median improvement; no per-level significance claimed)
  pooled over 6 levels x 4 starts = 24 pairs: Wilcoxon signed-rank; Spearman(kappa1 r p99 of the level, log(err_other / err_BO))
  'BO's relative advantage increases with complexity' only if, at EVERY evaluated budget n >= 256 (here 256 and 512):
      Spearman rho > 0 with p < 0.05  AND  C3-level BO wins >= 7 of 8 pairs
  n to reach 1 %, 0.5 %, 0.2 %, 0.1 % (median over starts, first grid n)
Version: python analyze_analytic_r2.py r2   (evaluation_r2, n <= 512 -> analysis_r2/)
         python analyze_analytic_r2.py r3   (evaluation_r3, n <= 1000 -> analysis_r3/)
Secondary (descriptive, no new claims): the same pooled tests on mean_highcurv (top-10 % curvature references) and max; calibration (cov90, rms_z);
controller cost from the recorded step wall times (campaigns_r2 steps_cont.jsonl) and evaluation fit times.
-> analysis_<ver>/*.csv|json|png|pdf|svg"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
VER = sys.argv[1] if len(sys.argv) > 1 else "r2"; NMAX = {"r2": 512, "r3": 1000}[VER]
OUT = HERE / f"analysis_{VER}"; OUT.mkdir(exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11, "axes.spines.top": False, "axes.spines.right": False, "savefig.dpi": 300})
COL = {"BO": "#2a78d6", "AXMM": "#eb6834", "RANDOM": "#1baf7a", "BOFIXED": "#6b6b6b"}
LAB = {"BO": "BO+MES (A20)", "AXMM": "axial maximin", "RANDOM": "random", "BOFIXED": "BO-FIXED (diagnostic)"}
DESC = json.load(open(HERE / "analytic_descriptors.json")); KAP = {k: v["kappa_max_dimless_p99"] for k, v in DESC.items()}; ORDER = sorted(KAP, key=KAP.get)
BUD = tuple(b for b in (8, 16, 33, 64, 128, 256, 512, 1000) if b <= NMAX)


def save(fig, name):
    for e in ("png", "pdf", "svg"):
        fig.savefig(OUT / f"{name}.{e}", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main():
    M = pd.read_csv(HERE / f"evaluation_{VER}/metrics_all.csv"); M = M[M.n <= NMAX]; A = M[M.surrogate == "A20"].copy()
    cov = A.groupby(["level", "arm"]).start.nunique().unstack(); cov.to_csv(OUT / "design_coverage.csv")
    rows, pooled = [], []
    for n in sorted(A.n.unique()):
        for L in ORDER:
            g = A[(A.level == L) & (A.n == n)].pivot_table(index="start", columns="arm", values="mean")
            for o in ("AXMM", "RANDOM", "BOFIXED"):
                if "BO" not in g or o not in g:
                    continue
                p = g[["BO", o]].dropna(); d = p[o] - p.BO
                rows.append(dict(n=n, level=L, cls=DESC[L]["class"], kappa_p99=KAP[L], vs=o, pairs=len(p), bo_wins=int((d > 0).sum()),
                                 median_improvement_pp=100 * float(np.median(d)), median_log_ratio=float(np.median(np.log(p[o] / p.BO))),
                                 bo_median_pct=100 * float(p.BO.median()), other_median_pct=100 * float(p[o].median())))
        for o, met in (("AXMM", "mean"), ("RANDOM", "mean"), ("AXMM", "mean_highcurv"), ("RANDOM", "mean_highcurv"), ("AXMM", "max"), ("RANDOM", "max")):
            g = A[A.n == n].pivot_table(index=["level", "start"], columns="arm", values=met)
            if "BO" not in g or o not in g:
                continue
            g = g[["BO", o]].dropna(); d = g[o] - g.BO; lr = np.log(g[o] / g.BO); kp = np.array([KAP[l] for l, _ in g.index]); rho, pr = stats.spearmanr(kp, lr)
            c3 = g[[DESC[l]["class"] == "C3_localized" for l, _ in g.index]]
            pooled.append(dict(n=n, vs=o, metric=met, pairs=len(g), bo_wins=int((d > 0).sum()), wilcoxon_p=float(stats.wilcoxon(d).pvalue) if np.any(d != 0) else np.nan,
                               median_log_ratio=float(np.median(lr)), median_ratio_other_over_bo=float(np.exp(np.median(lr))), spearman_rho=float(rho), spearman_p=float(pr),
                               c3_pairs=len(c3), c3_bo_wins=int((c3[o] > c3.BO).sum())))
    P = pd.DataFrame(rows); P.to_csv(OUT / "paired_by_level.csv", index=False); PO = pd.DataFrame(pooled); PO.to_csv(OUT / "paired_pooled_and_trend.csv", index=False)
    crit = {}
    for o in ("AXMM", "RANDOM"):
        big = PO[(PO.vs == o) & (PO.metric == "mean") & (PO.n >= 256)]
        crit[o] = dict(budgets_evaluated=[int(x) for x in big.n], rows=big.to_dict("records"),
                       SUPPORTED=bool(len(big) > 0 and all(r.spearman_rho > 0 and r.spearman_p < 0.05 and r.c3_bo_wins >= 7 and r.c3_pairs == 8 for r in big.itertuples())),
                       note=("PROTOCOL.md specifies n >= 256 including 1000; this is the n <= 512 interim checkpoint (not the final test)." if NMAX < 1000 else "all protocol budgets n >= 256 (256, 512, 1000)."))
    th = []
    for L in ORDER:
        for arm in ("BO", "AXMM", "RANDOM", "BOFIXED"):
            g = A[(A.level == L) & (A.arm == arm)].groupby("n")["mean"].median()
            for t in (0.01, 0.005, 0.002, 0.001):
                hit = g[g < t]; th.append(dict(level=L, arm=arm, threshold_pct=100 * t, n_first=int(hit.index.min()) if len(hit) else None))
    T = pd.DataFrame(th); T.to_csv(OUT / "n_to_threshold.csv", index=False)
    summ = A[A.n.isin(BUD)].groupby(["level", "arm", "n"])[["mean", "max", "mean_highcurv", "max_highcurv", "cov90", "rms_z", "mean_sd"]].median().reset_index()
    summ.to_csv(OUT / "summary_median.csv", index=False)
    cal = A[A.n.isin(BUD)].groupby(["arm", "n"])[["cov90", "rms_z", "corr_sd_abs_err"]].agg(["median", "min", "max"]); cal.to_csv(OUT / "calibration.csv")
    ct = []
    for f in sorted(HERE.glob("campaigns_r2/*/BO_s*/steps_cont.jsonl")):
        S = [json.loads(l) for l in open(f)]
        for a, b in zip(S, S[1:]):
            ct.append(dict(level=f.parent.parent.name, start=int(f.parent.name[4:]), n=b["test_number"], step_wall_s=b["wall_unix"] - a["wall_unix"]))
    CT = pd.DataFrame(ct); CT = CT[(CT.step_wall_s > 0) & (CT.step_wall_s < 3600)]
    CT["bin"] = pd.cut(CT.n, [0, 256, 384, 512, 640, 768, 896, 1000]); cb = CT.groupby("bin", observed=True).step_wall_s.agg(["count", "median", "max"]); cb.to_csv(OUT / "bo_step_cost.csv")
    m_ = CT.n > 150; slope, icpt = np.polyfit(np.log(CT.n[m_]), np.log(CT.step_wall_s[m_]), 1)
    ft = M.groupby(["surrogate", "n"]).fit_s.median(); ft.to_csv(OUT / "evaluation_fit_cost.csv")
    cost = dict(bo_step_wall_s_fit=dict(exponent=float(slope), prefactor=float(np.exp(icpt)), note="wall s per BO step vs n, 24 concurrent campaigns on 12 cores / 24 threads"),
                bo_step_wall_s_by_n=cb.reset_index().astype({"bin": str}).to_dict("records"), bo_total_recorded_step_wall_h=float(CT.step_wall_s.sum() / 3600),
                nonadaptive_design_cost="negligible (geometric designs; analytic oracle)")
    json.dump(dict(version=VER, nmax=NMAX, cost=cost, criterion=crit, descriptors={L: dict(kappa_p99=KAP[L], cls=DESC[L]["class"]) for L in ORDER}), open(OUT / "conclusion.json", "w"), indent=1, default=float)
    # figures
    fig, ax = plt.subplots(2, 3, figsize=(16, 9), sharex=True)
    for a, L in zip(ax.ravel(), ORDER):
        for arm in ("BO", "AXMM", "RANDOM", "BOFIXED"):
            g = A[(A.level == L) & (A.arm == arm)].groupby("n")["mean"]
            if g.ngroups == 0:
                continue
            m, lo, hi = g.median() * 100, g.min() * 100, g.max() * 100
            a.fill_between(m.index, lo, hi, color=COL[arm], alpha=0.12, lw=0); a.plot(m.index, m, "-o", ms=3, color=COL[arm], lw=1.8 if arm != "BOFIXED" else 1.0,
                                                                                     ls="-" if arm != "BOFIXED" else "--", label=LAB[arm])
        a.set(xscale="log", yscale="log", title=f"{L}  ({DESC[L]['class']}, kappa1*r p99 = {KAP[L]:.1f})", xticks=list(BUD)); a.set_xticklabels([str(b) for b in BUD])
        a.grid(alpha=0.3)
    for a in ax[1]:
        a.set_xlabel("observations n")
    for a in ax[:, 0]:
        a.set_ylabel("area-weighted mean |rel. error| (%)")
    ax[0, 0].legend(frameon=False, fontsize=9)
    fig.suptitle(f"EXPLORATORY analytic complexity-scaling study: median over 4 matched starts (band: min-max), 20,000 exact references, n <= {NMAX}", fontsize=12)
    fig.tight_layout(); save(fig, "error_vs_n_by_level")
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.6))
    for a, o in zip(ax, ("AXMM", "RANDOM")):
        for n, c in [x for x in ((33, "#cccccc"), (128, "#999999"), (256, "#666666"), (512, "#333333"), (1000, "#000000")) if x[0] <= NMAX]:
            s = P[(P.vs == o) & (P.n == n)].set_index("level").reindex(ORDER)
            a.plot([KAP[l] for l in ORDER], np.exp(s.median_log_ratio), "-o", color=c, label=f"n = {n}")
        a.axhline(1, color="#c00", lw=1); a.set(xlabel="complexity: kappa1*r p99 of the level", ylabel=f"median error ratio {o} / BO", title=f"relative advantage of BO over {LAB[o]} (>1 = BO better)")
        a.legend(frameon=False, fontsize=9); a.grid(alpha=0.3)
    fig.tight_layout(); save(fig, "advantage_vs_complexity")
    print(cov.to_string()); print(PO.round(4).to_string()); print(json.dumps({o: c["SUPPORTED"] for o, c in crit.items()})); print(json.dumps(cost, indent=1, default=str)[:1500])


if __name__ == "__main__":
    main()
