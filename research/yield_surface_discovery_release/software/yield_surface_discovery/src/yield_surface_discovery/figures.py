"""Release figures, rebuilt from the bundled frozen results (data/frozen_results).  Nothing is re-run.

Terminology: 'test' = one real physical experiment / oracle call (symmetry pseudo-observations never count); global p95 relative error;
conventional loading library (uniaxial + aligned biaxial, Pool 1); expanded physical loading library (+ pure shear + combined normal-shear, Pool 2);
arbitrary-direction reference (Pool 3; a relaxed reference, not an upper bound); generic blind (no analytical template, no symmetry);
physics-informed (same method + known shear-mirror symmetry); predictive uncertainty = maximum predictive-variance (entropy) sampling."""
import json
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from .geometry import angle_rad, lonlat_deg
from .symmetry import MIRROR_MATRIX

INK, INK2, MUTED, SURF = "#0b0b0b", "#52514e", "#c3c2b7", "#fcfcfb"
LIBC = {"conventional": "#2a78d6", "expanded physical": "#eb6834", "arbitrary-direction reference": "#8a8984"}
LIBN = {"conventional": "conventional loading library", "expanded physical": "expanded physical loading library", "arbitrary-direction reference": "arbitrary-direction reference"}
MODE = {"physics-informed": "#6250d6", "generic blind": "#008300"}
POS, NEG = "#e34948", "#2a78d6"                                                   # diverging pair (worse / better)
FAMC = {"uniaxial": "#2a78d6", "biaxial (aligned cruciform)": "#eb6834", "in-plane shear": "#1baf7a", "combined normal+shear (axis-aligned)": "#eda100",
        "combined normal+shear (rotated)": "#e87ba4"}
BUD = [8, 12, 16, 20, 24, 28]
plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": "#e6e5e1", "grid.linewidth": 0.6,
                     "axes.edgecolor": INK2, "xtick.color": INK2, "ytick.color": INK2, "legend.frameon": False, "figure.facecolor": SURF,
                     "axes.facecolor": SURF, "savefig.facecolor": SURF, "pdf.fonttype": 42})


class Bundle:
    def __init__(self, d):
        self.d = Path(d); r = pd.read_csv(self.d / "runs.csv"); self.runs = r[r.method == "mesh"].copy()

    def csv(self, name):
        return pd.read_csv(self.d / name)

    def sel(self, arm, lib=None):
        x = self.runs[self.runs.arm_id == arm]; return x if lib is None else x[x.library == lib]


def _save(fig, out, name, vals):
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out / f"{name}.pdf", bbox_inches="tight"); fig.savefig(out / f"{name}.png", dpi=200, bbox_inches="tight"); plt.close(fig)
    pd.DataFrame(vals).to_csv(out / f"{name}_plotted_values.csv", index=False); return out / f"{name}.png"


def _series(ax, g, col, color, label, vals, fid, ls="-", marker="o", off=0.0, faint=True, scale=100, lw=2.0):
    if faint:
        for s, gg in g.groupby("start"):
            gg = gg.sort_values("n_tests"); ax.plot(gg.n_tests + off, scale * gg[col], color=color, lw=0.6, alpha=0.3, ls=ls)
    m = g.groupby("n_tests")[col].median()
    ax.plot(m.index + off, scale * m.values, color=color, lw=lw, ls=ls, marker=marker, ms=4, mfc=color if marker != "s" else SURF, mec=color, label=label)
    for _, r in g.iterrows():
        vals.append(dict(figure=fid, series=label, start=r.start, n_real_tests=r.n_tests, metric=col, value=r[col]))
    for n, v in m.items():
        vals.append(dict(figure=fid, series=label, start="median", n_real_tests=n, metric=col, value=v))
    return m


def _tests_axis(ax, log=True, ylabel="global p95 relative error (%)"):
    ax.set_xticks(BUD); ax.set_xlabel("number of real tests"); ax.set_ylabel(ylabel)
    if log:
        ax.set_yscale("log")


# ---------------------------------------------------------------------------------------------------------------- Figure 1
def fig1(B, out):
    vals = []; F = B.csv("fuhg_same_points.csv"); F = F[F.information_matched_no_symmetry & ~F.model.str.contains("continuous")]
    S = [("Fuhg et al. GP (released checkpoint)", "Fuhg GP", "#52514e", "o"), ("Fuhg et al. SVR (released checkpoint)", "Fuhg SVR", "#8a8984", "s"),
         ("Fuhg et al. ICNN (released checkpoint)", "Fuhg ICNN", "#b5b4ad", "^"), ("direct Matern-3/2 GP (generic blind)", "direct Matérn-3/2 GP (ours, generic blind)", MODE["generic blind"], "D"),
         ("template-informed GP (Hill + Matern residual), no symmetry", "template-informed GP (ours, earlier; secondary)", "#9085e9", "v")]
    fig, axs = plt.subplots(1, 3, figsize=(10.4, 3.4), constrained_layout=True)
    for ax, (c, t) in zip(axs, (("abs_p95", "(a) global p95 relative error"), ("weak10_max_over", "(b) weak-10% region: max overprediction"), ("sharp10_max_over", "(c) sharp-10% region: max overprediction"))):
        for k, (m, lab, col, mk) in enumerate(S):
            g = F[F.model == m].sort_values("N"); x = np.arange(3) + (k - 2) * 0.09
            ax.plot(x, 100 * g[c], color=col, marker=mk, ms=5.5 if "ours" in lab else 4.5, lw=2.0 if "generic blind" in lab else 1.0, ls="-" if "secondary" not in lab else ":", label=lab,
                    mfc=col if mk != "v" else SURF)
            for _, r in g.iterrows():
                vals.append(dict(figure="fig1", model=lab, N=int(r.N), metric=c, value=r[c]))
        ax.set_xticks(range(3)); ax.set_xticklabels(["8", "16", "28"]); ax.set_xlabel("real tests (published fixed design)")
        ax.set_ylabel("%"); ax.set_yscale("log"); ax.set_title(t, loc="left", fontsize=8)
    h, l = axs[0].get_legend_handles_labels(); fig.legend(h, l, loc="upper center", bbox_to_anchor=(0.5, 1.12), ncol=3, fontsize=7)
    fig.suptitle("Figure 1  Reconstruction on identical published test points, information-matched (no symmetry for any model; single run per design)",
                 x=0.01, ha="left", y=1.19, fontsize=9)
    return _save(fig, out, "fig1_fuhg_same_points", vals)


# ---------------------------------------------------------------------------------------------------------------- Figure 2
def _gapmap(ax, C, title, cmap="Blues"):
    lo, la = np.meshgrid(np.linspace(-180, 180, 241), np.linspace(-90, 90, 121))
    G = np.stack([np.cos(np.radians(la)) * np.cos(np.radians(lo)), np.cos(np.radians(la)) * np.sin(np.radians(lo)), np.sin(np.radians(la))], -1).reshape(-1, 3)
    CC = np.vstack([C, C @ MIRROR_MATRIX.T]); gap = np.degrees(angle_rad(G, CC).min(1)).reshape(lo.shape)
    im = ax.pcolormesh(lo, la, gap, cmap=cmap, vmin=0, vmax=40, shading="auto", rasterized=True); L, Bt = lonlat_deg(C)
    ax.scatter(L, Bt, s=4, color=INK, lw=0); ax.axhline(45, color=INK2, lw=0.7, ls="--"); ax.axhline(-45, color=INK2, lw=0.7, ls="--")
    ax.set_xlim(-180, 180); ax.set_ylim(-90, 90); ax.set_xticks([-180, -90, 0, 90, 180]); ax.set_yticks([-90, -45, 0, 45, 90]); ax.grid(False)
    ax.set_xlabel("in-plane angle atan2(z2, z1) (deg)"); ax.set_ylabel("shear latitude asin(z3) (deg)"); ax.set_title(title, loc="left", fontsize=8); return im


def fig2(B, out, cand_dir):
    vals = []; fig = plt.figure(figsize=(11.0, 6.9), constrained_layout=True); gs = fig.add_gridspec(2, 3, width_ratios=[1, 1, 0.95]); libs = list(LIBC)
    lib1 = pd.read_csv(Path(cand_dir) / "conventional_loading_library.csv", float_precision="round_trip"); lib2 = pd.read_csv(Path(cand_dir) / "expanded_physical_loading_library.csv", float_precision="round_trip")
    im = _gapmap(fig.add_subplot(gs[0, 0]), lib1[["z1", "z2", "z3"]].values, "(a) conventional library, 140 tests:\nangular gap to the nearest test or its mirror")
    _gapmap(fig.add_subplot(gs[0, 1]), lib2[["z1", "z2", "z3"]].values, "(a′) expanded physical library, 416 tests: same")
    cb = fig.colorbar(im, ax=fig.axes[:2], location="bottom", shrink=0.45, pad=0.02); cb.set_label("angular gap (deg); dashed: |z3|/|z| = 1/√2, the shear limit of uniaxial/biaxial tests")
    ax = axb = fig.add_subplot(gs[0, 2]); bands = ["band_0-0.25_abs_p95", "band_0.25-0.5_abs_p95", "band_0.5-0.707_abs_p95", "band_0.707-0.85_abs_p95", "band_0.85-1_abs_p95"]
    for i, l in enumerate(libs):
        g = B.sel("blind_study:directGP_mirror_on_frozen_sequences", l); g = g[g.n_tests == 28]
        med = [100 * g[b].median() for b in bands]; ax.bar(np.arange(5) + (i - 1) * 0.27, med, width=0.25, color=LIBC[l], edgecolor=SURF, label=LIBN[l])
        for b, v in zip(bands, med):
            vals.append(dict(figure="fig2", panel="B", library=l, band=b, median_p95_pct=v))
    ax.axvline(2.5, color=INK2, lw=0.8, ls="--"); ax.text(2.42, 0.97, "beyond the reach of\nuniaxial/biaxial tests →", transform=ax.get_xaxis_transform(), fontsize=6.3, color=INK2, ha="right", va="top")
    ax.set_xticks(range(5)); ax.set_xticklabels(["0–0.25", "0.25–0.5", "0.5–0.71", "0.71–0.85", "0.85–1"], fontsize=7); ax.set_xlabel("shear content |z3|/|z| of held-out direction")
    ax.set_ylabel("p95 relative error (%)"); ax.set_title("(b) error by shear content, 28 tests\n(median of 5 starts)", loc="left", fontsize=8)
    z = np.load(Path(B.d) / "per_ray_example_directGP_mirror_s0_n28.npz"); L, Bt = lonlat_deg(z["eval"])
    for j, (p, l) in enumerate((("pool1", "conventional"), ("pool2", "expanded physical"))):
        ax = fig.add_subplot(gs[1, j]); e = 100 * np.abs(z[f"rel_{p}"])
        sc = ax.scatter(L, Bt, c=e, s=5, cmap="Reds", vmin=0, vmax=12, lw=0, rasterized=True); ax.axhline(45, color=INK2, lw=0.7, ls="--"); ax.axhline(-45, color=INK2, lw=0.7, ls="--")
        ax.set_xlim(-180, 180); ax.set_ylim(-90, 90); ax.set_xticks([-180, -90, 0, 90, 180]); ax.set_yticks([-90, -45, 0, 45, 90]); ax.grid(False)
        ax.set_xlabel("in-plane angle atan2(z2, z1) (deg)"); ax.set_ylabel("shear latitude asin(z3) (deg)")
        ax.set_title(f"(c{'′' if j else ''}) |relative error| at 28 tests, {LIBN[l]}\n(one start; held-out directions; global p95 {np.quantile(e, .95):.1f}%)", loc="left", fontsize=7.5)
    fig.colorbar(sc, ax=fig.axes[-2:], location="bottom", shrink=0.45, pad=0.02).set_label("|relative radial error| (%)")
    ax = fig.add_subplot(gs[1, 2]); cov = B.csv("library_coverage.csv")
    W = [("eval_reachable_within_10deg", "≤ 10°"), ("eval_reachable_within_5deg", "≤ 5°"), ("eval_reachable_envelope", "within shear\nenvelope")]
    for k, (c, lab) in enumerate(W):
        for i, l in enumerate(libs):
            v = float(cov[cov.library == l][c].iloc[0]); ax.bar(k + (i - 1) * 0.27, 100 * v, width=0.25, color=LIBC[l], edgecolor=SURF)
            ax.annotate(f"{100 * v:.0f}", (k + (i - 1) * 0.27, 100 * v), xytext=(0, 2), textcoords="offset points", ha="center", fontsize=6)
            vals.append(dict(figure="fig2", panel="D", library=l, statistic=c, value=v))
    ax.set_xticks(range(3)); ax.set_xticklabels([w[1] for w in W], fontsize=7); ax.set_ylabel("% of held-out directions"); ax.set_ylim(0, 112)
    ax.set_xlabel("held-out directions near a test (or its mirror)"); ax.set_title("(d) coverage statistics", loc="left", fontsize=8)
    h, l = axb.get_legend_handles_labels(); fig.legend(h, l, loc="upper center", bbox_to_anchor=(0.5, 1.035), ncol=3, fontsize=7.2)
    fig.suptitle("Figure 2  Why conventional tests leave an error floor: they cannot reach high-shear directions (physics-informed direct GP on frozen sequences; mesh-scored)",
                 x=0.01, ha="left", y=1.075, fontsize=9)
    return _save(fig, out, "fig2_action_space_mechanism", vals)


# ---------------------------------------------------------------------------------------------------------------- Figure 3
def fig3(B, out):
    vals = []; fig, axs = plt.subplots(1, 2, figsize=(10.0, 3.7), constrained_layout=True)
    ax = axs[0]; g = B.sel("same_surrogate:C", "expanded physical")
    _series(ax, g, "abs_p95", MODE["physics-informed"], "frozen practical method, expanded physical library (5 starts; 4 distinct designs)", vals, "fig3")
    c28 = B.sel("blind_study:directGP_mirror_on_frozen_sequences", "conventional"); f = 100 * c28[c28.n_tests == 28].abs_p95.median()
    ax.axhline(f, color=LIBC["conventional"], lw=1.0, ls=":"); ax.annotate(f"conventional-library floor ({f:.1f}%, panel b)", (8, f), xytext=(0, 3), textcoords="offset points", fontsize=6.5, color=INK2)
    _tests_axis(ax); ax.set_title("(a) frozen practical (physics-informed) method — registered", loc="left", fontsize=8); ax.legend(fontsize=6.6, loc="upper right")
    ax = axs[1]
    for l, off, ls in (("conventional", -0.35, "-"), ("expanded physical", 0.0, "-"), ("arbitrary-direction reference", 0.35, "--")):
        m = _series(ax, B.sel("blind_study:directGP_mirror_on_frozen_sequences", l), "abs_p95", LIBC[l], LIBN[l], vals, "fig3", ls=ls, off=off)
        ax.annotate(f"{100 * m.iloc[-1]:.1f}%", (28 + off, 100 * m.iloc[-1]), xytext=(5, 0), textcoords="offset points", fontsize=6.5, va="center", color=INK2)
    _tests_axis(ax); ax.set_title("(b) effect of the loading library (same direct GP + mirror; frozen symmetric sequences)", loc="left", fontsize=8); ax.legend(fontsize=6.6, loc="upper right")
    fig.suptitle("Figure 3  Adaptive discovery of the Mg-4%Li yield surface: direct Matérn-3/2 GP on log radius + known shear-mirror symmetry (thin lines = starts, thick = median)",
                 x=0.01, ha="left", y=1.05, fontsize=9)
    return _save(fig, out, "fig3_main_adaptive_result", vals)


# ---------------------------------------------------------------------------------------------------------------- Figure 4
def effects(B):
    """four ingredient changes, each from ONE registered experiment (difference of medians, global p95, percentage points)."""
    def med(arm, lib):
        return B.sel(arm, lib).groupby("n_tests").abs_p95.median()
    E = [("conventional instead of expanded loading library", med("blind_study:directGP_mirror_on_frozen_sequences", "conventional") - med("blind_study:directGP_mirror_on_frozen_sequences", "expanded physical"),
          "blind study, registered decomposition: direct GP + mirror on frozen symmetric sequences"),
         ("no shear-mirror symmetry", med("same_surrogate:A", "expanded physical") - med("same_surrogate:C", "expanded physical"), "same-surrogate ablation (registered): Arm A − Arm C"),
         ("IVR instead of predictive uncertainty", med("same_surrogate:D", "expanded physical") - med("same_surrogate:C", "expanded physical"), "same-surrogate ablation (registered): Arm D − Arm C"),
         ("Hill-1948 template added (residual GP)", med("pool_experiment:template_informed", "expanded physical") - med("blind_study:directGP_mirror_on_frozen_sequences", "expanded physical"),
          "blind study, registered decomposition: identical tests and symmetry")]
    return E


def fig4(B, out):
    vals = []; E = effects(B); fig = plt.figure(figsize=(10.4, 3.9), constrained_layout=True); gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.15])
    for k, n in enumerate((16, 28)):
        ax = fig.add_subplot(gs[0, k]); y = np.arange(len(E))[::-1]
        for yi, (lab, s, src) in zip(y, E):
            v = 100 * float(s[n]); ax.barh(yi, v, color=POS if v > 0 else NEG, height=0.6, edgecolor=SURF)
            ax.annotate(f"{v:+.1f}", (v, yi), xytext=(4 if v >= 0 else -4, 0), textcoords="offset points", ha="left" if v >= 0 else "right", va="center", fontsize=7)
            vals.append(dict(figure="fig4", n_real_tests=n, change=lab, delta_global_p95_pp=v, source=src))
        ax.axvline(0, color=INK2, lw=0.8); ax.set_yticks(y); ax.set_yticklabels([e[0] for e in E] if k == 0 else [], fontsize=7.5)
        lim = max(abs(100 * float(s[n])) for _, s, _ in E) * 1.35; ax.set_xlim(-lim, lim); ax.grid(axis="y", visible=False)
        ax.set_xlabel("Δ median global p95 (pp)\n> 0: worse than the frozen method"); ax.set_title(f"({'ab'[k]}) at {n} real tests", loc="left", fontsize=8)
    ax = fig.add_subplot(gs[0, 2]); sty = ["-", "-", "--", ":"]; mk = ["s", "o", "^", "v"]; cols = [LIBC["conventional"], MODE["physics-informed"], INK2, "#9085e9"]
    for (lab, s, src), ls, m, c in zip(E, sty, mk, cols):
        s = s.reindex([12, 16, 20, 24, 28]); ax.plot(s.index, 100 * s.values, ls=ls, marker=m, color=c, lw=1.5, ms=4, label=lab)
    ax.axhline(0, color=INK2, lw=0.8); ax.set_xticks([12, 16, 20, 24, 28]); ax.set_xlabel("number of real tests"); ax.set_ylabel("Δ median global p95 (pp)")
    ax.set_title("(c) the same changes across budgets", loc="left", fontsize=8); ax.legend(fontsize=6.3, loc="upper right")
    fig.suptitle("Figure 4  What matters: loading space and known symmetry dominate; the analytical template and IVR-vs-predictive-uncertainty do not (each bar from one registered comparison)",
                 x=0.01, ha="left", y=1.06, fontsize=9)
    return _save(fig, out, "fig4_what_matters", vals)


# ---------------------------------------------------------------------------------------------------------------- Figure 5
def fig5(B, out):
    vals = []; fig, axs = plt.subplots(1, 4, figsize=(11.6, 3.3), constrained_layout=True, gridspec_kw=dict(width_ratios=[1.2, 1, 1, 0.75]))
    G = {"generic blind": B.sel("same_surrogate:A", "expanded physical"), "physics-informed": B.sel("same_surrogate:C", "expanded physical")}
    lab = {"generic blind": "generic blind (no symmetry)", "physics-informed": "physics-informed (+ shear mirror)"}
    for ax, (c, t, log) in zip(axs[:3], (("abs_p95", "(a) global p95 relative error", True), ("weak10_max_over", "(b) weak-10% region: max overprediction", True),
                                         ("sharp10_max_over", "(c) sharp-10% region: max overprediction", True))):
        for (k, g), off in zip(G.items(), (-0.25, 0.25)):
            _series(ax, g, c, MODE[k], lab[k], vals, "fig5", off=off)
        _tests_axis(ax, log, "%"); ax.set_title(t, loc="left", fontsize=8)
    axs[0].legend(fontsize=6.8, loc="upper right"); ax = axs[3]
    for i, (k, g) in enumerate(G.items()):
        v = 100 * g[g.n_tests == 28]["high_shear_0.85_1_abs_p95"]; ax.scatter(np.full(len(v), i) + np.linspace(-0.08, 0.08, len(v)), v, color=MODE[k], s=16, alpha=0.8, lw=0)
        ax.plot([i - 0.22, i + 0.22], [v.median()] * 2, color=MODE[k], lw=2.4); ax.annotate(f"median {v.median():.2f}%", (i + 0.26, v.median()), xytext=(4, 0), textcoords="offset points", fontsize=6.5, va="center")
        for s_, x in zip(g[g.n_tests == 28].start, v):
            vals.append(dict(figure="fig5", series=lab[k], start=s_, n_real_tests=28, metric="high_shear_p95", value=x / 100))
    ax.set_xticks([0, 1]); ax.set_xticklabels(["generic\nblind", "physics-\ninformed"]); ax.set_xlim(-0.5, 1.7); ax.set_ylabel("%"); ax.set_title("(d) high-shear p95, 28 tests", loc="left", fontsize=8)
    fig.suptitle("Figure 5  Generic blind vs physics-informed: identical algorithm and first 8 tests; only the shear-mirror symmetry differs (registered same-surrogate Arms A and C)",
                 x=0.01, ha="left", y=1.06, fontsize=9)
    return _save(fig, out, "fig5_blind_vs_physics_informed", vals)


# ---------------------------------------------------------------------------------------------------------------- Figure 6
def fig6(B, out):
    fig, ax = plt.subplots(figsize=(11.2, 4.9)); ax.set_xlim(0, 121); ax.set_ylim(0, 49); ax.axis("off")

    def box(x, y, w, h, text, fc="#f0efec", ec=INK2, ls="-", fs=7.4):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.35,rounding_size=1.2", fc=fc, ec=ec, lw=1.1, ls=ls))
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs, color=INK, linespacing=1.25)

    def arrow(x0, y0, x1, y1, txt=None, ls="-", col=INK2, rad=0.0, dy=1.3):
        ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=11, lw=1.2, color=col, ls=ls, connectionstyle=f"arc3,rad={rad}"))
        if txt:
            ax.text((x0 + x1) / 2, (y0 + y1) / 2 + dy, txt, ha="center", fontsize=6.6, color=INK2)
    V, VE = "#eceaf9", "#6250d6"
    box(1, 30, 20, 10, "candidate physical\nloading directions\n(conventional or\nexpanded library)")
    box(26, 30, 19, 10, "8 angular-maximin\ninitial tests\n(geometry only)")
    box(50, 30, 20, 10, "next loading path\n(requested\ndirection d)", fc=V, ec=VE)
    box(77, 29, 43, 12, "YieldOracle.evaluate(d)\nlaboratory test  |  MOOSE  |  other solver  |  analytical model\n"
        "returns: stress at first yield, yield radius,\nachieved direction, status, metadata", fc=SURF, ls="--", fs=7.0)
    box(1, 11, 20, 10, "select the maximum\n(deterministic ties);\nrepeat until the budget\nof REAL tests")
    box(26, 11, 19, 10, "predictive\nuncertainty over\nuntested candidates", fc=V, ec=VE)
    box(50, 11, 20, 10, "direct Matérn-3/2 GP\non log yield radius\n(constant mean,\nno analytical template)", fc=V, ec=VE)
    box(77, 11, 43, 10, "yield observation\nlog r(d) = log |z at first yield|")
    box(50, 1, 20, 6.5, "optional, physics-informed:\nknown shear-mirror symmetry\n(free pseudo-observations)", fc=SURF, ec=VE, ls=":", fs=6.6)
    arrow(21, 35, 26, 35); arrow(45, 35, 50, 35, "tests 1–8"); arrow(70, 35, 77, 35, "execute")
    arrow(98.5, 29, 98.5, 21); arrow(77, 16, 70, 16, "update"); arrow(50, 16, 45, 16); arrow(26, 16, 21, 16)
    arrow(11, 21, 55, 30, "tests 9, 10, …", rad=-0.15, dy=2.4); arrow(60, 7.5, 60, 11, ls=":", col=VE)
    ax.text(1, 47.3, "Figure 6  Frozen practical workflow: the controller sees only candidate directions and oracle observations; the solver interface is the only coupling.",
            fontsize=8.6, color=INK)
    ax.text(1, 44.6, "Generic blind mode = the same loop without the dotted symmetry box.  Test count always = real experiments (oracle calls); pseudo-observations are free.",
            fontsize=7.0, color=INK2)
    return _save(fig, out, "fig6_workflow", [dict(figure="fig6", note="schematic; no data")])


# ---------------------------------------------------------------------------------------------------------------- supplementary
def figS1(B, out):
    vals = []; fig, ax = plt.subplots(figsize=(6.2, 3.4), constrained_layout=True); cols = {"A": "#2a78d6", "B": "#eb6834", "C": "#1baf7a"}
    for k, ls in (("A", "-"), ("B", "--"), ("C", "-")):
        g = B.sel(f"init_experiment:{k}", "conventional"); name = g.initialization.iloc[0] if len(g) else k
        _series(ax, g, "abs_p95", cols[k], f"{k}: {name}".replace("-> MES", "-> predictive uncertainty"), vals, "figS1", ls=ls, faint=(k == "C"))
    _tests_axis(ax); ax.legend(fontsize=6.6); ax.set_title("S1  Initialisation (registered, template-informed era; conventional library;\n4 Hill-calibration tests + 4 tests by rule A/B/C, then predictive uncertainty)", loc="left", fontsize=8)
    return _save(fig, out, "figS1_initialization", vals)


def figS2(B, out):
    vals = []; fig, axs = plt.subplots(1, 2, figsize=(8.4, 3.3), constrained_layout=True, sharey=True)
    for ax, arms, t, col in ((axs[0], (("A", "predictive uncertainty"), ("B", "IVR")), "(a) generic blind (no symmetry)", MODE["generic blind"]),
                             (axs[1], (("C", "predictive uncertainty"), ("D", "IVR")), "(b) physics-informed (shear mirror)", MODE["physics-informed"])):
        for (a, lab), ls, mk, off in zip(arms, ("-", "--"), ("o", "s"), (-0.2, 0.2)):
            _series(ax, B.sel(f"same_surrogate:{a}", "expanded physical"), "abs_p95", col, f"{lab} (Arm {a})", vals, "figS2", ls=ls, marker=mk, off=off)
        _tests_axis(ax); ax.set_title(t, loc="left", fontsize=8); ax.legend(fontsize=7)
    fig.suptitle("S2  Predictive uncertainty vs IVR from the same surface GP (registered; tie by the pre-specified rule -> predictive uncertainty kept)", x=0.01, ha="left", y=1.06, fontsize=8.5)
    return _save(fig, out, "figS2_entropy_vs_ivr", vals)


def figS3(B, out):
    vals = []; S = B.csv("response_dependence_summary.csv"); S = S[S.variant == "shuffle"]; fig, axs = plt.subplots(1, 2, figsize=(8.6, 3.2), constrained_layout=True)
    sty = {"A_entropy_nosym": (MODE["generic blind"], "-", "o", "generic blind, predictive uncertainty"), "B_ivr_nosym": (MODE["generic blind"], "--", "s", "generic blind, IVR"),
           "C_entropy_sym": (MODE["physics-informed"], "-", "o", "physics-informed, predictive uncertainty"), "D_ivr_sym": (MODE["physics-informed"], "--", "s", "physics-informed, IVR")}
    for ax, (c, t) in zip(axs, (("share_genuine_changes", "(a) next test changes when measured values are shuffled"), ("spearman_median", "(b) rank correlation of candidate scores"))):
        for arm, (col, ls, mk, lab) in sty.items():
            g = S[S.arm == arm].sort_values("next_test"); ax.plot(g.next_test, g[c], color=col, ls=ls, marker=mk, ms=4, mfc=col if mk == "o" else SURF, lw=1.4, label=lab)
            for _, r in g.iterrows():
                vals.append(dict(figure="figS3", arm=arm, next_test=r.next_test, metric=c, value=r[c]))
        ax.axvline(21, color=INK2, lw=0.8, ls=":"); ax.set_xticks([9, 13, 17, 21, 25]); ax.set_xlabel("next test chosen"); ax.set_title(t, loc="left", fontsize=8)
    axs[0].legend(fontsize=6.4); axs[0].annotate("hyper-parameter learning\nstarts (no symmetry)", (21.3, 0.05), fontsize=6.3, color=INK2)
    fig.suptitle("S3  Selections respond to measured yield values only through GP hyper-parameter learning (diagnostic copies; registered prediction confirmed)", x=0.01, ha="left", y=1.06, fontsize=8.5)
    return _save(fig, out, "figS3_response_dependence", vals)


def figS4(B, out):
    vals = []; fig, axs = plt.subplots(1, 2, figsize=(9.0, 3.3), constrained_layout=True, gridspec_kw=dict(width_ratios=[1.2, 1]))
    _series(axs[0], B.sel("same_surrogate:C", "expanded physical"), "abs_p95", MODE["physics-informed"], "frozen practical method (Arm C; ordinary maximin; 5 starts)", vals, "figS4", off=-0.2)
    _series(axs[0], B.runs[B.runs.arm_id == "cleanup:cleaned"], "abs_p95", "#8a8984", "symmetry-aware maximin + deduplication + deterministic ties (10 starts; not adopted)", vals, "figS4", ls="--", off=0.2)
    _tests_axis(axs[0]); axs[0].legend(fontsize=6.4); axs[0].set_title("(a) global p95", loc="left", fontsize=8)
    V = B.csv("cleanup_vs_armC.csv"); ax = axs[1]; mets = [("weak10_max_over", "weak-10% max over"), ("sharp10_max_over", "sharp-10% max over"), ("high_shear_0.85_1_abs_p95", "high-shear p95")]
    for i, (c, lab) in enumerate(mets):
        for j, n in enumerate((16, 28)):
            r = V[(V.metric == c) & (V.n_tests == n)].iloc[0]; x = i * 3 + j
            ax.bar(x - 0.18, 100 * r.armC_median, width=0.34, color=MODE["physics-informed"], edgecolor=SURF, label="Arm C (frozen)" if x == 0 else None)
            ax.bar(x + 0.18, 100 * r.cleaned_median, width=0.34, color="#8a8984", edgecolor=SURF, label="cleaned (not adopted)" if x == 0 else None)
            vals.append(dict(figure="figS4", metric=c, n_real_tests=n, armC=r.armC_median, cleaned=r.cleaned_median))
    ax.set_xticks([0, 1, 3, 4, 6, 7]); ax.set_xticklabels(["16", "28"] * 3); ax.set_ylabel("% (median)")
    for i, (c, lab) in enumerate(mets):
        ax.text(i * 3 + 0.5, -0.17, lab, transform=ax.get_xaxis_transform(), ha="center", fontsize=6.6, color=INK2)
    ax.legend(fontsize=6.6); ax.set_title("(b) regional errors at 16 and 28 tests", loc="left", fontsize=8)
    fig.suptitle("S4  Symmetry-aware initialisation cleanup (registered): better at 8 tests, non-inferior globally, worse weak/high-shear at 16 -> Arm C kept", x=0.01, ha="left", y=1.06, fontsize=8.5)
    return _save(fig, out, "figS4_symmetry_aware_cleanup", vals)


def figS5(B, out):
    vals = []; Q = B.csv("same_surrogate_sequences.csv"); fams = list(FAMC); fig, axs = plt.subplots(1, 2, figsize=(9.6, 2.9), constrained_layout=True, sharey=True)
    for ax, (arm, t) in zip(axs, (("C_entropy_sym", "(a) physics-informed (Arm C)"), ("A_entropy_nosym", "(b) generic blind (Arm A)"))):
        g = Q[Q.arm == arm]
        for _, r in g.iterrows():
            f = r.family.replace("F1 ", "").replace("F2 ", "").replace("F3 ", ""); ax.add_patch(plt.Rectangle((r.test_number - 0.45, r.start - 0.4), 0.9, 0.8, color=FAMC[f], lw=0))
            vals.append(dict(figure="figS5", arm=arm, start=r.start, test_number=r.test_number, label=r.label, family=f))
        ax.axvline(8.5, color=INK, lw=1.0); ax.set_xlim(0.4, 28.6); ax.set_ylim(-0.6, 4.6); ax.set_yticks(range(5)); ax.set_yticklabels([f"start {s}" for s in range(5)])
        ax.set_xticks([1, 4, 8, 12, 16, 20, 24, 28]); ax.set_xlabel("test number (1-8 angular maximin | 9-28 predictive uncertainty)"); ax.grid(False); ax.set_title(t, loc="left", fontsize=8)
    for f in fams:
        axs[1].scatter([], [], marker="s", color=FAMC[f], label=f)
    fig.legend(*axs[1].get_legend_handles_labels(), loc="upper center", bbox_to_anchor=(0.5, 1.14), ncol=5, fontsize=6.6)
    fig.suptitle("S5  Selected test families, expanded physical library (registered same-surrogate runs)", x=0.01, ha="left", y=1.22, fontsize=8.5)
    return _save(fig, out, "figS5_test_family_sequences", vals)


def figS6(B, out):
    vals = []; Z = np.load(Path(B.d) / "tsai_wu_sections.npz"); D = B.csv("tsai_wu_diagnostic.csv"); fig, axs = plt.subplots(2, 2, figsize=(9.4, 6.0), constrained_layout=True)
    S = Z["section_dirs"]; th = np.degrees(np.arctan2(S[:, 1], S[:, 0])) % 360; o = np.argsort(th); th = th[o]
    for row, (t, nm) in enumerate((("tsai_wu_t0", "Tsai-Wu T300"), ("tsai_wu_t1", "Tsai-Wu E-glass"))):
        ax = axs[row, 0]; ref = Z[f"{t}__sec_ref"][o]; bl = Z[f"{t}__blind_sec_t"][o]; sd = Z[f"{t}__blind_sec_sd"][o]
        ax.fill_between(th, bl * np.exp(-1.96 * sd), bl * np.exp(1.96 * sd), color=MODE["generic blind"], alpha=0.13, lw=0)
        ax.plot(th, ref, color=INK, lw=1.7, label="reference"); ax.plot(th, Z[f"{t}__template_bank_unseen_sec_t"][o], color="#9085e9", lw=1.2, label="template-informed (held-out-family bank)")
        ax.plot(th, bl, color=MODE["generic blind"], lw=1.2, ls="--", label="direct GP, generic blind (±1.96 sd)")
        ax.plot(th, Z[f"{t}__posthoc_directGP_mirror_sec_t"][o], color=MODE["physics-informed"], lw=1.0, ls=":", label="direct GP + mirror (post hoc)")
        A = Z[f"{t}__acq_dirs"]; near = np.abs(A[:, 2]) <= np.sin(np.radians(10)); ang = np.degrees(np.arctan2(A[:, 1], A[:, 0])) % 360
        ax.scatter(ang[near], Z[f"{t}__acq_r"][near], marker="x", color=INK, s=20, zorder=5, label="tests ≤ 10° from the section")
        ax.set_yscale("log"); ax.set_xlim(0, 360); ax.set_xticks(range(0, 361, 45)); ax.set_xlabel("angle in the (s_xx, s_yy) plane (deg)"); ax.set_ylabel("normalised yield radius")
        ax.set_title(f"({'ac'[row]}) {nm}: s_xy = 0 section, 50 tests (one saved run)", loc="left", fontsize=8)
        if row == 0:
            ax.legend(fontsize=6.0, loc="upper right")
        ax = axs[row, 1]; d = D[(D.target == t) & D.setup.str.startswith("analytical")]
        for model, col, mk in (("blind", MODE["generic blind"], "o"), ("template_bank_unseen", "#9085e9", "s"), ("blindGP_sym (post-hoc)", MODE["physics-informed"], "^")):
            g = d[d.model == model].sort_values("N"); ax.plot(g.N, 100 * g.p95_median, color=col, marker=mk, lw=1.4, label={"blind": "direct GP, generic blind", "template_bank_unseen": "template-informed bank",
                                                                                                                       "blindGP_sym (post-hoc)": "direct GP + mirror (post hoc)"}[model])
            for _, r in g.iterrows():
                vals.append(dict(figure="figS6", target=t, model=model, N=r.N, p95_median=r.p95_median))
        ax.set_xscale("log"); ax.set_xticks([10, 50, 100]); ax.set_xticklabels(["10", "50", "100"]); ax.set_xlabel("number of tests (same saved points for all models)")
        ax.set_ylabel("global p95 relative error (%)"); ax.set_title(f"({'bd'[row]}) {nm}: identical points, median of 6 runs", loc="left", fontsize=8); ax.legend(fontsize=6.3)
    fig.suptitle("S6  Tsai-Wu limitation: a smooth stationary GP cannot follow needle-like surfaces (representation-limited; the template is not the cause)", x=0.01, ha="left", y=1.03, fontsize=8.5)
    return _save(fig, out, "figS6_tsai_wu_limitation", vals)


FIGS = dict(fig1=fig1, fig2=fig2, fig3=fig3, fig4=fig4, fig5=fig5, fig6=fig6, figS1=figS1, figS2=figS2, figS3=figS3, figS4=figS4, figS5=figS5, figS6=figS6)


def make_all(results_dir, out_dir, only=None):
    B = Bundle(results_dir); cand = Path(results_dir).parent / "candidates"; paths = []
    for k, f in FIGS.items():
        if only and k not in only:
            continue
        paths.append(f(B, out_dir, cand) if k == "fig2" else f(B, out_dir))
    return paths
