#!/usr/bin/env python3
"""H1 figures from saved outputs only.  Palette (validated): V2 #2a78d6, maximin #eb6834, V3 #1f9e74, random neutral grey (own marker);
V3_ISO (diagnostic) neutral dark dashed / rings.   ->  figures/fig_h1_*.png|pdf"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = Path(__file__).resolve().parent; sys.dont_write_bytecode = True; sys.path.insert(0, str(HERE.parents[1] / "complexity_scaling_2026"))
import cs_common as C                                                                # noqa: E402
FIG = HERE / "figures"; FIG.mkdir(exist_ok=True); INK, INK2, GRID, BG = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"; BO, AX, V3c, RA = "#2a78d6", "#eb6834", "#1f9e74", "#8a8984"
plt.rcParams.update({"font.size": 8.5, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "figure.facecolor": BG, "axes.facecolor": BG, "savefig.facecolor": BG, "legend.frameon": False})
M = pd.read_csv(HERE / "evaluation/metrics_all.csv"); RES = json.load(open(HERE / "analysis/results.json")); PR = pd.read_csv(HERE / "analysis/per_realisation.csv")
STY = [("V3", V3c, "-", "o", "V3 (frozen v3.0-LRV)"), ("V3_ISO", INK, "--", "o", "V3 driven by isotropic J2 responses (diagnostic; re-observed)"), ("V2", BO, "-", "x", "V2-A20 BO+MES"),
       ("AXMM", AX, "-", "s", "angular maximin"), ("RANDOM", RA, "-", "D", "random")]


def save(fig, nm):
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"{nm}.{ext}", dpi=170, bbox_inches="tight")
    plt.close(fig); print(nm)


fig, axs = plt.subplots(1, 3, figsize=(14.5, 4.3), gridspec_kw=dict(wspace=0.28))
for a, (met, lab) in zip(axs, (("glob_mean", "global mean |error| (%)"), ("hc_mean", "high-curvature-region mean |error| (%)"), ("glob_max", "global max |error| (%)"))):
    for arm, col, ls, mk, l in STY:
        g = M[M.arm == arm].groupby("n")[met]; md, q1, q3 = 100 * g.median(), 100 * g.quantile(0.25), 100 * g.quantile(0.75)
        a.plot(md.index, md, color=col, ls=ls, marker=mk, ms=3.5, lw=1.6, label=l); a.fill_between(md.index, q1, q3, color=col, alpha=0.10, lw=0)
    a.set_yscale("log"); a.set_xlabel("observations n"); a.set_ylabel(lab); a.axvline(33, color=GRID, lw=1.2); a.axvline(64, color=GRID, lw=1.2)
h, l = axs[0].get_legend_handles_labels(); fig.legend(h, l, loc="lower center", ncol=5, bbox_to_anchor=(0.5, -0.08), fontsize=7.8)
fig.suptitle("Holdout H1 (24 independent Yld91 surfaces): median and interquartile band of the reconstruction error vs number of experiments", fontsize=9.5)
save(fig, "fig_h1_error_vs_n")
fig, axs = plt.subplots(1, 3, figsize=(14.5, 4.2), gridspec_kw=dict(wspace=0.3))
P = M[M.n.isin([33, 64])].pivot_table(index=["n", "realisation"], columns="arm", values="glob_mean")
comps = [("AXMM", 33, AX, "H1  vs maximin, n = 33"), ("V2", 33, BO, "H2  vs V2, n = 33"), ("AXMM", 64, AX, "H3  vs maximin, n = 64"), ("V3_ISO", 33, INK, "attribution: vs V3_ISO, n = 33")]
for k, (o, n, col, lab) in enumerate(comps):
    r = (P.loc[n][o] / P.loc[n]["V3"]).to_numpy(); xx = k + np.random.default_rng(k).uniform(-0.12, 0.12, len(r))
    axs[0].scatter(xx, r, s=16, color=col, facecolor=col if o != "V3_ISO" else "none", zorder=3)
    key = {"H1": "H1", "H2": "H2", "H3": "H3"}; h = ["H1", "H2", "H3", None][k]
    s = RES["confirmatory"][h] if h else RES["response_attribution_V3_vs_V3ISO_glob33"]
    axs[0].errorbar(k + 0.3, s["hl_ratio"], yerr=[[s["hl_ratio"] - s["hl_ci95"][0]], [s["hl_ci95"][1] - s["hl_ratio"]]], fmt="D", color=INK, ms=4, capsize=3)
axs[0].axhline(1, color=INK2, lw=0.8); axs[0].axhline(1.10, color=INK2, lw=0.6, ls=":"); axs[0].set_yscale("log"); axs[0].set_xticks(range(4))
axs[0].set_xticklabels([c[3] for c in comps], fontsize=7, rotation=12); axs[0].set_ylabel("error ratio other / V3 (global mean; > 1: V3 better)")
axs[0].set_title("Paired ratios (dots); HL ratio + exact 95 % CI\n(dotted: practical threshold 1.10)", fontsize=9)
lv = PR.copy(); lv["mlab"] = lv.m.map({2: "m = 2", 6: "m = 6", 8: "m = 8"})
for k, m in enumerate((2, 6, 8)):
    for s_, mk in ((0.15, "o"), (0.30, "^")):
        d = lv[(lv.m == m) & (lv.s == s_)]; axs[1].scatter(k + (0.15 if s_ == 0.30 else -0.15) + np.zeros(len(d)), np.exp(d.L_H1), marker=mk, color=V3c, s=26,
                                                            facecolor=V3c if s_ == 0.15 else "none", label=(f"anisotropy spread s = {s_}" if k == 0 else None))
axs[1].axhline(1, color=INK2, lw=0.8); axs[1].set_yscale("log"); axs[1].set_xticks(range(3)); axs[1].set_xticklabels(["m = 2 (quadratic)", "m = 6", "m = 8"])
axs[1].set_ylabel("maximin / V3, global mean, n = 33"); axs[1].legend(fontsize=7.5); axs[1].set_title(f"Advantage depends on Hosford m\n(Spearman ρ = {RES['moderators_spearman_L_H1']['m']['rho']:.2f}); anisotropy: ρ = {RES['moderators_spearman_L_H1']['aniso']['rho']:.2f}", fontsize=9)
axs[2].scatter(PR.aniso, np.exp(PR.L_H1), color=V3c, s=22, label="vs maximin"); axs[2].scatter(PR.aniso, np.exp(PR.L_vs_ISO), facecolor="none", edgecolor=INK, s=26, label="vs V3_ISO (attribution)")
axs[2].axhline(1, color=INK2, lw=0.8); axs[2].set_yscale("log"); axs[2].set_xlabel("anisotropy r max / r min"); axs[2].set_ylabel("other / V3, global mean, n = 33"); axs[2].legend(fontsize=7.5)
axs[2].set_title("V3 beats maximin, but not the same rule driven by\nisotropic responses: the gain is not material-specific", fontsize=9)
save(fig, "fig_h1_paired_and_moderators")
# selection map, realisation 0 (Lambert, axis hemisphere)
def proj(D):
    D = np.atleast_2d(D) * np.where(np.atleast_2d(D)[:, [2]] < 0, -1.0, 1.0); f = np.sqrt(2.0 / (1.0 + D[:, 2])); return D[:, 0] * f, D[:, 1] * f
LIB = C.library().directions; fig, axs = plt.subplots(1, 3, figsize=(13, 4.6))
for a, i in zip(axs, (0, 8, 16)):
    for arm, col, mk, fc, s, lab in (("AXMM", AX, "s", "none", 22, "angular maximin"), ("V2", BO, "x", BO, 18, "V2 BO+MES"), ("V3_ISO", INK, "o", "none", 60, "V3 on isotropic responses"), ("V3", V3c, "o", V3c, 18, "V3")):
        X = LIB[json.load(open(HERE / f"runs/{arm}/r{i:02d}/design.json"))["axes"]]
        a.scatter(*proj(X[8:33]), marker=mk, s=s, facecolor=fc, edgecolor=col, color=col if mk == "x" else None, lw=1.0, label=lab + " 9–33", zorder=4 if arm == "V3" else 3)
    a.scatter(*proj(X[:8]), marker="s", s=22, color=INK, label="shared initial 8", zorder=5); eb = proj(np.array([[1, 1, 0]]) / np.sqrt(2)); a.scatter(*eb, marker="*", s=120, color=INK2, zorder=6, label="equal biaxial")
    t = np.linspace(0, 2 * np.pi, 300); a.plot(np.sqrt(2) * np.cos(t), np.sqrt(2) * np.sin(t), color=INK2, lw=0.8); a.set_aspect("equal"); a.axis("off")
    pop = json.load(open(HERE / "population.json"))["realisations"][i]; a.set_title(f"realisation {i} ({pop['level']})", fontsize=9)
h, l = axs[0].get_legend_handles_labels(); fig.legend(h, l, loc="lower center", ncol=6, bbox_to_anchor=(0.5, -0.02), fontsize=7.8)
fig.suptitle("Selections 9–33 (Lambert map of loading axes; centre: pure σxy shear; rim: σxy = 0, incl. equal biaxial ★).  V3 and V3_ISO pick different axes but share a bias toward the\n"
             "biaxial cap (24 % of picks within 30° of equal biaxial vs 12 % for maximin; area share 13 %)", fontsize=9, y=1.02)
save(fig, "fig_h1_selection")
