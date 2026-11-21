#!/usr/bin/env python3
"""Pre-specified analysis (PROTOCOL.md) of the analytic complexity-scaling study -> analysis/*.csv|json|png|pdf|svg."""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import cs_common as C                                                              # noqa: E402
OUT = HERE / "analysis"; OUT.mkdir(exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11, "axes.spines.top": False, "axes.spines.right": False, "savefig.dpi": 300})
COL = {"BO": "#2a78d6", "AXMM": "#eb6834", "RANDOM": "#1baf7a", "BOFIXED": "#6b6b6b"}; LAB = {"BO": "BO+MES (A20)", "AXMM": "axial maximin", "RANDOM": "random", "BOFIXED": "BO-FIXED (diagnostic)"}
DESC = json.load(open(HERE / "analytic_descriptors.json")); KAP = {k: v["kappa_max_dimless_p99"] for k, v in DESC.items()}
ORDER = sorted(KAP, key=KAP.get); BUD = C.BUDGETS


def save(fig, name):
    for e in ("png", "pdf", "svg"):
        fig.savefig(OUT / f"{name}.{e}", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main():
    M = pd.read_csv(HERE / "evaluation/metrics_all.csv"); A = M[M.surrogate == "A20"].copy(); A["kappa_p99"] = A.level.map(KAP); A["cls"] = A.level.map(lambda l: DESC[l]["class"])
    # ---------------- paired statistics
    rows, pooled, trend = [], [], []
    for n in sorted(A.n.unique()):
        for L in ORDER:
            g = A[(A.level == L) & (A.n == n)].pivot_table(index="start", columns="arm", values="mean")
            for other in ("AXMM", "RANDOM", "BOFIXED"):
                if other not in g or "BO" not in g:
                    continue
                p = g[["BO", other]].dropna(); d = p[other] - p.BO; lr = np.log(p[other] / p.BO)
                rows.append(dict(n=n, level=L, cls=DESC[L]["class"], kappa_p99=KAP[L], vs=other, pairs=len(p), bo_wins=int((d > 0).sum()), median_improvement_pp=100 * float(np.median(d)),
                                 median_log_ratio=float(np.median(lr)), bo_median=100 * float(p.BO.median()), other_median=100 * float(p[other].median()),
                                 wilcoxon_p=float(stats.wilcoxon(d).pvalue) if len(p) >= 2 and np.any(d != 0) else np.nan))
    P = pd.DataFrame(rows); P.to_csv(OUT / "paired_by_level.csv", index=False)
    for n in sorted(A.n.unique()):
        for other in ("AXMM", "RANDOM"):
            s = P[(P.n == n) & (P.vs == other)]
            if s.empty:
                continue
            g = A[A.n == n].pivot_table(index=["level", "start"], columns="arm", values="mean")[["BO", other]].dropna(); d = g[other] - g.BO; lr = np.log(g[other] / g.BO)
            kp = np.array([KAP[l] for l, _ in g.index]); rho, pr = stats.spearmanr(kp, lr)
            c3 = g[[DESC[l]["class"] == "C3_localized" for l, _ in g.index]]
            pooled.append(dict(n=n, vs=other, pairs=len(g), bo_wins=int((d > 0).sum()), wilcoxon_p=float(stats.wilcoxon(d).pvalue) if np.any(d != 0) else np.nan,
                               median_log_ratio=float(np.median(lr)), spearman_rho_kappa_vs_logratio=float(rho), spearman_p=float(pr),
                               c3_pairs=len(c3), c3_bo_wins=int((c3[other] > c3.BO).sum())))
    PO = pd.DataFrame(pooled); PO.to_csv(OUT / "paired_pooled_and_trend.csv", index=False)
    crit = {}
    for other in ("AXMM", "RANDOM"):
        big = PO[(PO.vs == other) & (PO.n >= 256)]
        crit[other] = dict(per_budget=big[["n", "spearman_rho_kappa_vs_logratio", "spearman_p", "c3_bo_wins", "c3_pairs"]].to_dict("records"),
                           SUPPORTED=bool(len(big) > 0 and all((r.spearman_rho_kappa_vs_logratio > 0) and (r.spearman_p < 0.05) and (r.c3_bo_wins >= 7) for r in big.itertuples())))
    # ---------------- n to reach thresholds (median over starts)
    th = []
    for L in ORDER:
        for arm in ("BO", "AXMM", "RANDOM", "BOFIXED"):
            med = A[(A.level == L) & (A.arm == arm)].groupby("n")["mean"].median()
            if med.empty:
                continue
            for t in (0.01, 0.005, 0.002, 0.001):
                ok = med[med < t]; th.append(dict(level=L, arm=arm, threshold_pct=100 * t, n_required=int(ok.index.min()) if len(ok) else None))
    TH = pd.DataFrame(th); TH.to_csv(OUT / "n_to_threshold.csv", index=False)
    # ---------------- response dependence
    seqs = {}
    for L in ORDER:
        for k in range(4):
            for arm in ("BO", "AXMM", "BOFIXED"):
                f = HERE / "campaigns" / ("ANY" if arm == "BOFIXED" else L) / f"{arm}_s{k}" / "design.json"
                if f.exists():
                    seqs[(L if arm != "BOFIXED" else "ANY", arm, k)] = json.load(open(f))["axes"]
    fd = lambda a, b: next((i + 1 for i, (x, y) in enumerate(zip(a, b)) if x != y), None)
    resp = []
    for k in range(4):
        for i, L1 in enumerate(ORDER):
            if (L1, "BO", k) not in seqs:
                continue
            r = dict(start=k, level=L1, first_div_vs_BOFIXED=fd(seqs[(L1, "BO", k)], seqs.get(("ANY", "BOFIXED", k), [])),
                     first_div_vs_AXMM=fd(seqs[(L1, "BO", k)], seqs.get((L1, "AXMM", k), [])), shared_axes_with_AXMM_at_1000=len(set(seqs[(L1, "BO", k)]) & set(seqs.get((L1, "AXMM", k), []))))
            for L2 in ORDER:
                if L2 != L1 and (L2, "BO", k) in seqs:
                    r[f"first_div_vs_{L2}"] = fd(seqs[(L1, "BO", k)], seqs[(L2, "BO", k)])
            resp.append(r)
    R = pd.DataFrame(resp); R.to_csv(OUT / "response_dependence.csv", index=False)
    # ---------------- hyper-parameters, acquisition, cost (from BO fits)
    fits = []
    for L in ORDER:
        for k in range(4):
            f = HERE / "campaigns" / L / f"BO_s{k}" / "fits.csv"
            if f.exists():
                x = pd.read_csv(f); x["level"], x["start"] = L, k
                s = pd.read_csv(HERE / "campaigns" / L / f"BO_s{k}" / "steps.csv"); x["selected_latent_var"] = s.latent_var.to_numpy()[: len(x)]
                fits.append(x)
    FI = pd.concat(fits) if fits else pd.DataFrame(); FI.to_csv(OUT / "bo_fits_all.csv", index=False)
    cost = {L: dict(bo_wall_s=[json.load(open(HERE / "campaigns" / L / f"BO_s{k}/run_info.json"))["wall_s"] for k in range(4) if (HERE / "campaigns" / L / f"BO_s{k}/run_info.json").exists()]) for L in ORDER}
    # ---------------- figures
    fig, axs = plt.subplots(2, 3, figsize=(16, 9), sharex=True, sharey=True)
    for ax, L in zip(axs.ravel(), ORDER):
        for arm in ("RANDOM", "AXMM", "BOFIXED", "BO"):
            g = A[(A.level == L) & (A.arm == arm)].groupby("n")["mean"]
            if g.ngroups == 0:
                continue
            med, lo, hi = g.median() * 100, g.min() * 100, g.max() * 100
            ax.fill_between(med.index, lo, hi, color=COL[arm], alpha=0.12, lw=0); ax.plot(med.index, med, "-" if arm != "BOFIXED" else "--", color=COL[arm], lw=2 if arm == "BO" else 1.6, label=LAB[arm])
        ax.set_xscale("log"); ax.set_yscale("log"); ax.set_title(f"{L}\n{DESC[L]['class']}, κ₁r p99 = {KAP[L]:.2f}", fontsize=11.5)
        for b in BUD:
            ax.axvline(b, color="#eeeeee", lw=0.8, zorder=0)
    for ax in axs[1]:
        ax.set_xlabel("number of observations n")
    for ax in axs[:, 0]:
        ax.set_ylabel("mean |relative error| at 20,000 refs [%]")
    axs[0, 0].legend(frameon=False, fontsize=9.5); fig.suptitle("Analytic stress tests (exact ground truth, synthetic oracle): error vs budget, median over 4 starts (band = min–max)", fontsize=13.5)
    save(fig, "error_vs_budget_by_level")
    fig, axs = plt.subplots(1, 3, figsize=(16, 4.8))
    for ax, n in zip(axs, (33, 256, 1000)):
        for arm in ("RANDOM", "AXMM", "BO"):
            g = A[(A.n == n) & (A.arm == arm)].groupby("level")["mean"].median() * 100
            ax.plot([KAP[l] for l in ORDER if l in g], [g[l] for l in ORDER if l in g], "-o", color=COL[arm], label=LAB[arm])
        ax.set_yscale("log"); ax.set_xlabel("measured complexity: κ₁·r p99"); ax.set_title(f"n = {n}"); ax.set_ylabel("median mean |error| [%]")
    axs[0].legend(frameon=False); fig.suptitle("Error vs measured yield-surface complexity", fontsize=13.5); save(fig, "error_vs_complexity")
    fig, axs = plt.subplots(1, 2, figsize=(14, 5))
    for ax, other in zip(axs, ("AXMM", "RANDOM")):
        for n, mk in zip((33, 128, 256, 512, 1000), ("o", "s", "^", "D", "*")):
            g = A[A.n == n].pivot_table(index=["level", "start"], columns="arm", values="mean")
            if other not in g or "BO" not in g:
                continue
            g = g[["BO", other]].dropna(); x = np.array([KAP[l] for l, _ in g.index]); y = np.log(g[other] / g.BO).to_numpy()
            ax.scatter(x + (np.log(n) - 5) * 0.04, y, marker=mk, s=40, label=f"n = {n}", alpha=0.8)
        ax.axhline(0, color="#52514e", lw=1); ax.set_xlabel("measured complexity κ₁·r p99"); ax.set_ylabel(f"log(err_{other} / err_BO)   (> 0: BO better)")
        ax.set_title(f"BO+MES vs {LAB[other]}: paired log error ratio per start"); ax.legend(frameon=False, fontsize=9)
    save(fig, "paired_ratio_vs_complexity")
    fig, axs = plt.subplots(1, 3, figsize=(16, 4.8))
    for arm in ("RANDOM", "AXMM", "BO"):
        g = A[A.arm == arm].groupby(["level", "n"])[["mean_highcurv", "mean_sd", "cov90"]].median().reset_index()
        for L in ORDER[-2:]:
            h = g[g.level == L]; axs[0].plot(h.n, 100 * h.mean_highcurv, ("-" if L == ORDER[-1] else ":"), color=COL[arm], label=f"{LAB[arm]} {L}")
            axs[1].plot(h.n, 100 * h.mean_sd, ("-" if L == ORDER[-1] else ":"), color=COL[arm]); axs[2].plot(h.n, h.cov90, ("-" if L == ORDER[-1] else ":"), color=COL[arm])
    axs[0].set_title("ACTUAL error in the top-10% curvature region [%]"); axs[1].set_title("PREDICTED mean GP std at the references [%]"); axs[2].set_title("90% interval coverage (nominal 0.90)")
    axs[2].axhline(0.9, color="#52514e", lw=1)
    for a in axs:
        a.set_xscale("log"); a.set_xlabel("n")
    axs[0].set_yscale("log"); axs[1].set_yscale("log"); axs[0].legend(frameon=False, fontsize=8); fig.suptitle(f"Most complex levels ({ORDER[-2]} dotted, {ORDER[-1]} solid): error vs uncertainty vs calibration", fontsize=13)
    save(fig, "highcurv_uncertainty_calibration")
    if len(FI):
        fig, axs = plt.subplots(1, 4, figsize=(19, 4.6))
        for L in ORDER:
            h = FI[FI.level == L].groupby("n_real")[["length_scale", "amplitude", "white_noise", "fit_s", "selected_latent_var"]].median()
            axs[0].plot(h.index, h.length_scale, label=L); axs[1].plot(h.index, h.white_noise); axs[2].plot(h.index, np.sqrt(h.selected_latent_var)); axs[3].plot(h.index, h.fit_s)
        axs[0].axhline(5, color="#52514e", ls="--", lw=1); axs[0].axhline(0.05, color="#52514e", ls="--", lw=1); axs[1].axhline(1e-6, color="#52514e", ls="--", lw=1)
        axs[0].set_title("ML-II length scale ℓ (bounds dashed)"); axs[1].set_title("white noise (lower bound dashed)"); axs[2].set_title("max σ_f at decision (acquisition)"); axs[3].set_title("GP fit time per step [s]")
        for a in axs:
            a.set_xscale("log"); a.set_xlabel("n")
        axs[1].set_yscale("log"); axs[2].set_yscale("log"); axs[3].set_yscale("log"); axs[0].legend(frameon=False, fontsize=8); save(fig, "hyperparameters_acquisition_cost")
    json.dump(dict(criterion=crit, descriptors={k: {a: b for a, b in v.items() if a != "corr_curve"} for k, v in DESC.items()}, cost=cost), open(OUT / "summary.json", "w"), indent=1, default=float)
    pd.set_option("display.width", 220)
    print(PO.round(4).to_string()); print(json.dumps(crit, indent=1, default=float))
    print(P[(P.vs == "AXMM") & P.n.isin(BUD)].pivot_table(index="level", columns="n", values="median_log_ratio").round(3).loc[ORDER].to_string())
    print(TH.pivot_table(index=["level", "arm"], columns="threshold_pct", values="n_required").to_string())
    print(R.to_string())


if __name__ == "__main__":
    main()
