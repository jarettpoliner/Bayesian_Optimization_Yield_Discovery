#!/usr/bin/env python3
"""Stage 5 figures (descriptive; the pre-specified numbers are in evaluation/*, produced by the frozen evaluate_final.py).
Reads evaluation/*.csv, campaigns/*, refs/*; writes figures/*.png.  Surrogates are the frozen GPv2 (antipodal) exactly as in evaluate_final."""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import json, sys, warnings
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
FB = Path(__file__).resolve().parent; BB = FB.parent; CL = BB.parent
for p in (FB, CL / "v2", CL, CL / "homog_benchmark"):
    sys.path.insert(0, str(p))
import designs_common as D                                                          # noqa: E402
from bomes_v2 import GPv2, SCHEDULES, ANTIPODAL                                     # noqa: E402
EV = FB / "evaluation"; FIG = FB / "figures"; FIG.mkdir(exist_ok=True)
COL = {"BO": "#2a78d6", "AXMM": "#eb6834", "RANDOM": "#1baf7a"}; NAME = {"BO": "BO (V2-A20)", "AXMM": "Axial maximin", "RANDOM": "Random axes"}
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
REG = ("uniaxial", "pure_shear", "tension_compression", "biaxial_same_sign", "equal_biaxial")
plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2, "axes.grid": True,
                     "grid.color": GRID, "grid.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 130, "savefig.bbox": "tight",
                     "lines.linewidth": 2, "legend.frameon": False})
rt = lambda p: pd.read_csv(p, float_precision="round_trip")
M = rt(EV / "metrics_all.csv"); A = M[M.ref_set == "area_stratified"]; PR = rt(EV / "predictions_at_references.csv"); REF = {v: rt(EV / f"references_{v}.csv") for v in D.VARIANTS}


def save(fig, name):
    for ext, kw in ((".png", dict(dpi=300)), (".pdf", {}), (".svg", {})):   # presentation PNG + publication vectors
        fig.savefig((FIG / name).with_suffix(ext), **kw)
    plt.close(fig)


def campaign(v, meth, k):
    if meth == "BO":
        t = rt(FB / "campaigns" / v / f"bo_s{k}/tests.csv"); return t[["z1", "z2", "z3"]].to_numpy(), t.yield_radius.to_numpy()
    x = json.load(open(FB / "campaigns" / v / "nonadaptive/designs.json"))["designs"][f"axmm_s{k}" if meth == "AXMM" else f"random_{k}"]
    return np.array(x["directions"]), np.array(x["radii"])


def fitgp(X, r, sur="A20"):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore"); return GPv2(ANTIPODAL, SCHEDULES[sur]).fit(X, r)


# 1 error vs budget
for v in D.VARIANTS:
    fig, ax = plt.subplots(2, 2, figsize=(9, 6.2), sharex=True)
    for j, sur in enumerate(("A20", "FIXED")):
        for i, met in enumerate(("area_weighted", "mean")):
            a = ax[i, j]
            for m in COL:
                g = A[(A.variant == v) & (A.surrogate == sur) & (A.method == m)].groupby("n")[met]
                med, lo, hi = g.median() * 100, g.min() * 100, g.max() * 100
                a.fill_between(med.index, lo, hi, color=COL[m], alpha=0.13, lw=0); a.plot(med.index, med, color=COL[m], marker="o", ms=4, label=NAME[m])
            a.set_title(f"{v} - {sur} surrogate - {'area-weighted' if i == 0 else 'unweighted'} mean |rel. error|", fontsize=9, color=INK)
            a.set_ylabel("error (%)"); a.set_xticks((8, 12, 16, 20, 24, 28, 33))
    ax[1, 0].set_xlabel("physical observations n"); ax[1, 1].set_xlabel("physical observations n"); ax[0, 0].legend(loc="upper right")
    fig.suptitle(f"{v}: reconstruction error vs experiment budget (60 independent references; line = median over designs, band = min-max)", fontsize=10, color=INK)
    save(fig, f"error_vs_budget_{v}.png")
# 2 paired starts
fig, ax = plt.subplots(1, 4, figsize=(12, 3.6), sharey=False)
for i, (v, sur) in enumerate([(v, s) for v in D.VARIANTS for s in ("A20", "FIXED")]):
    g = A[(A.variant == v) & (A.surrogate == sur) & (A.n == 33)]; bo = g[g.method == "BO"].set_index("design").area_weighted.reindex(D.STARTS) * 100
    am = g[g.method == "AXMM"].set_index("design").area_weighted.reindex(D.STARTS) * 100; a = ax[i]
    for k in D.STARTS:
        a.plot([0, 1], [am[k], bo[k]], color=COL["BO"] if bo[k] < am[k] else COL["AXMM"], lw=1.2, alpha=0.8)
    a.scatter(np.zeros(9), am, color=COL["AXMM"], s=30, zorder=3, edgecolor="white", lw=1); a.scatter(np.ones(9), bo, color=COL["BO"], s=30, zorder=3, edgecolor="white", lw=1)
    a.set_xticks([0, 1], ["Axial maximin", "BO"]); a.set_xlim(-0.3, 1.3); a.set_title(f"{v}, {sur}: BO better in {(bo < am).sum()}/9", fontsize=9, color=INK)
    a.set_ylabel("area-weighted error at n=33 (%)" if i == 0 else "")
fig.suptitle("Paired initial designs (same first axis): primary endpoint", fontsize=10, color=INK); save(fig, "paired_starts.png")
# 3 regional
fig, ax = plt.subplots(1, 2, figsize=(11, 3.8))
for a, v in zip(ax, D.VARIANTS):
    g = A[(A.variant == v) & (A.surrogate == "A20") & (A.n == 33)]; w = 0.26
    for j, m in enumerate(COL):
        x = g[g.method == m][[f"mean_{r}" for r in REG]] * 100; med = x.median(); lo, hi = x.min(), x.max()
        a.bar(np.arange(5) + (j - 1) * w, med, w * 0.92, color=COL[m], label=NAME[m]); a.errorbar(np.arange(5) + (j - 1) * w, med, [med - lo, hi - med], fmt="none", ecolor=INK2, lw=0.8, capsize=2)
    a.set_xticks(range(5), [r.replace("_", "\n") for r in REG]); a.set_title(f"{v}: mean |rel. error| by loading region, n=33 (median, min-max)", fontsize=9, color=INK); a.set_ylabel("error (%)")
ax[0].legend(); save(fig, "regional_errors.png")
# 4 predicted vs measured
fig, ax = plt.subplots(1, 4, figsize=(13, 3.4))
for i, (v, m) in enumerate([(v, m) for v in D.VARIANTS for m in ("BO", "AXMM")]):
    p = PR[(PR.variant == v) & (PR.method == m) & (PR.surrogate == "A20") & (PR.n == 33)]; a = ax[i]; lo, hi = p.meas_r.min() * 0.95, p.meas_r.max() * 1.05
    a.plot([lo, hi], [lo, hi], color=INK2, lw=1); a.errorbar(p.meas_r, p.pred_r, yerr=1.96 * p.sd_log * p.pred_r, fmt="o", ms=2.5, color=COL[m], ecolor=COL[m], alpha=0.35, elinewidth=0.6)
    a.set_title(f"{v}, {NAME[m]} (9 designs, n=33)", fontsize=9, color=INK); a.set_xlabel("measured r (MPa)"); a.set_ylabel("GP predicted r (MPa)" if i == 0 else "")
fig.suptitle("Predicted vs measured yield radius at the 72 independent references (bars: 95% predictive interval)", fontsize=10, color=INK); save(fig, "predicted_vs_measured.png")
# 5 hyper-parameters / 6 acquisition
S = rt(EV / "bo_steps_all.csv")
fig, ax = plt.subplots(2, 3, figsize=(12, 6))
for i, v in enumerate(D.VARIANTS):
    for j, h in enumerate(("length_scale", "amplitude", "white_noise")):
        a = ax[i, j]
        for k in D.STARTS:
            s = S[(S.variant == v) & (S.start == k)]; a.plot(s.test_number, s[h], color=COL["BO"], alpha=0.45, lw=1.2)
        a.set_yscale("log"); a.axvline(20.5, color=INK2, lw=0.8, ls="--"); a.set_title(f"{v}: {h} (9 starts; ML-II from 20 real obs.)", fontsize=9, color=INK)
        a.set_xlabel("test number" if i == 1 else "")
save(fig, "hyperparameter_evolution.png")
fig, ax = plt.subplots(2, 2, figsize=(10, 6))
for i, v in enumerate(D.VARIANTS):
    for k in D.STARTS:
        s = S[(S.variant == v) & (S.start == k)]
        ax[i, 0].plot(s.test_number, s.latent_var, color=COL["BO"], alpha=0.45, lw=1.2); ax[i, 1].plot(s.test_number, s.runner_up_rel_gap, color=COL["BO"], alpha=0.45, lw=1.2)
    ax[i, 0].set_yscale("log"); ax[i, 1].set_yscale("log")
    ax[i, 0].set_title(f"{v}: selected max latent variance of log r", fontsize=9, color=INK); ax[i, 1].set_title(f"{v}: relative gap to runner-up candidate", fontsize=9, color=INK)
ax[1, 0].set_xlabel("test number"); ax[1, 1].set_xlabel("test number"); save(fig, "acquisition_scores.png")
# 7 GP maps (Lambert azimuthal, upper hemisphere z3 >= 0: antipodal symmetry makes it complete)
n = 161; u = np.linspace(-np.sqrt(2), np.sqrt(2), n); U, W = np.meshgrid(u, u); rho2 = U ** 2 + W ** 2; inside = rho2 <= 2
f = np.sqrt(np.clip(1 - rho2 / 4, 0, 1)); G = np.stack([U * f, W * f, 1 - rho2 / 2], -1)[inside]
lam = lambda P: (P[:, 0] / np.sqrt((1 + P[:, 2]) / 2), P[:, 1] / np.sqrt((1 + P[:, 2]) / 2))
for v in D.VARIANTS:
    fig, ax = plt.subplots(2, 3, figsize=(12, 7.4)); ref = REF[v]
    for i, m in enumerate(("BO", "AXMM")):
        X, r = campaign(v, m, 0); g = fitgp(X, r); mu, var = g.predict_log(G); Z = np.full(U.shape, np.nan); S_ = Z.copy(); Z[inside] = np.exp(mu); S_[inside] = np.sqrt(var)
        Xh = X * np.sign(X[:, 2:3] + 1e-15); xr, yr = lam(Xh); Pr = ref[["z1", "z2", "z3"]].to_numpy(); xq, yq = lam(Pr)
        mu_r, _ = g.predict_log(Pr); err = (np.exp(mu_r) / ref.r.to_numpy() - 1) * 100
        c0 = ax[i, 0].contourf(U, W, Z, 24, cmap="Blues"); fig.colorbar(c0, ax=ax[i, 0], label="predicted r (MPa)")
        c1 = ax[i, 1].contourf(U, W, S_, 24, cmap="Blues"); fig.colorbar(c1, ax=ax[i, 1], label="posterior sd of log r")
        for a in ax[i, :2]:
            a.scatter(xr, yr, s=14, c="white", edgecolor=INK, lw=0.8, label="observations (33)")
        lim = max(1, np.nanmax(np.abs(err))); c2 = ax[i, 2].scatter(xq, yq, c=err, cmap="RdBu_r", vmin=-lim, vmax=lim, s=26, edgecolor=INK2, lw=0.4)
        fig.colorbar(c2, ax=ax[i, 2], label="error at reference (%)")
        for a, t in zip(ax[i], ("mean", "uncertainty", "error at references")):
            a.set_aspect("equal"); a.set_xticks([]); a.set_yticks([]); a.grid(False); a.add_patch(plt.Circle((0, 0), np.sqrt(2), fill=False, color=INK2, lw=0.8))
            a.set_title(f"{NAME[m]} start 0: {t}", fontsize=9, color=INK)
    ax[0, 0].legend(loc="lower left", fontsize=7)
    fig.suptitle(f"{v}: GP posterior at n=33 on the stress-direction hemisphere (Lambert equal-area, centre = pure shear z3, rim = sigma12=0)", fontsize=10, color=INK)
    save(fig, f"gp_maps_{v}.png")
# 8 reconstruction along sigma12=0 (principal plane) + continuity sweep
fig, ax = plt.subplots(2, 2, figsize=(10.5, 8.6))
th = np.linspace(0, 2 * np.pi, 361); P0 = np.stack([np.cos(th), np.sin(th), 0 * th], 1)
for i, v in enumerate(D.VARIANTS):
    ref = REF[v]; sw = ref[ref.set == "continuity_sweep"].sort_values("sweep_t_deg"); a = ax[i, 0]
    t = np.radians(np.linspace(0, 90, 181)); e = np.array([1, 1, 0]) / np.sqrt(2); Q = np.cos(t)[:, None] * e + np.sin(t)[:, None] * np.array([0, 0, 1.0])
    for m in ("BO", "AXMM", "RANDOM"):
        ks = D.STARTS if m != "RANDOM" else D.RANDOM_SEEDS; curves = []; cq = []
        for k in ks:
            X, r = campaign(v, m, k); g = fitgp(X, r); curves.append(np.exp(g.predict_log(P0)[0])); cq.append(np.exp(g.predict_log(Q)[0]))
        c = np.median(curves, 0); a.plot(c * P0[:, 0] / np.sqrt(1), c * P0[:, 1], color=COL[m], lw=1.6, label=f"{NAME[m]} (median of designs)")
        ax[i, 1].plot(np.degrees(t), np.median(cq, 0), color=COL[m], lw=1.6, label=NAME[m]); ax[i, 1].fill_between(np.degrees(t), np.min(cq, 0), np.max(cq, 0), color=COL[m], alpha=0.12, lw=0)
    nr = ref[np.abs(ref.z3) < 0.2]; a.scatter(nr.r * nr.z1 / np.hypot(nr.z1, nr.z2), nr.r * nr.z2 / np.hypot(nr.z1, nr.z2), s=16, color=INK, zorder=4, label="references |z3|<0.2 (projected)")
    a.scatter(-nr.r * nr.z1 / np.hypot(nr.z1, nr.z2), -nr.r * nr.z2 / np.hypot(nr.z1, nr.z2), s=16, facecolor="none", edgecolor=INK, zorder=4, label="antipodal images")
    a.set_aspect("equal"); a.set_xlabel("z1 * r  (sigma11 direction, MPa)"); a.set_ylabel("z2 * r (MPa)"); a.set_title(f"{v}: n=33 reconstruction in the sigma12=0 plane", fontsize=9, color=INK)
    ax[i, 1].scatter(sw.sweep_t_deg, sw.r, color=INK, s=20, zorder=4, label="measured sweep references"); ax[i, 1].set_xlabel("t (deg): equal biaxial (0) -> pure shear (90)")
    ax[i, 1].set_ylabel("r (MPa)"); ax[i, 1].set_title(f"{v}: continuity sweep, n=33 (median, min-max over designs)", fontsize=9, color=INK)
ax[0, 0].legend(fontsize=7, loc="lower left"); ax[0, 1].legend(fontsize=7); save(fig, "reconstruction_principal_plane_and_sweep.png")
# 9 phase-specific plasticity diagnostics at references
fig, ax = plt.subplots(1, 3, figsize=(12.5, 3.6))
for j, v in enumerate(D.VARIANTS):
    ref = REF[v]; mk = "o" if v == "B1" else "s"; cc = COL["BO"] if v == "B1" else COL["AXMM"]
    xi = [REG.index(r) + (j - 0.5) * 0.25 + np.random.default_rng(1).uniform(-0.06, 0.06) for r in ref.region]
    ax[0].scatter(xi, ref.inclusion_share, s=14, marker=mk, color=cc, alpha=0.75, label=v)
    ax[1].scatter(ref.r1_matrix / ref.r, ref.r1_inclusion / ref.r, s=14, marker=mk, color=cc, alpha=0.75, label=v)
    ax[2].scatter(xi, ref.eqps_i / ref.eqps_m, s=14, marker=mk, color=cc, alpha=0.75, label=v)
ax[0].set_xticks(range(5), [r.replace("_", "\n") for r in REG], fontsize=7); ax[0].set_ylabel("inclusion share of RVE eqps at yield"); ax[0].set_title("phase share f_i eqps_i / eqps_RVE", fontsize=9)
ax[1].plot([0.5, 1.2], [0.5, 1.2], color=INK2, lw=0.8); ax[1].set_xlabel("first matrix yield r1_m / r"); ax[1].set_ylabel("first inclusion yield r1_i / r"); ax[1].set_title("first local yield per phase", fontsize=9)
ax[2].set_xticks(range(5), [r.replace("_", "\n") for r in REG], fontsize=7); ax[2].set_yscale("log"); ax[2].set_ylabel("eqps_i / eqps_m at yield"); ax[2].set_title("inclusion/matrix plastic strain ratio", fontsize=9)
nz = {v: int((REF[v].eqps_m == 0).sum()) for v in D.VARIANTS}
ax[2].text(0.02, 0.97, "\n".join(f"{v}: matrix eqps = 0 at yield in {nz[v]}/{len(REF[v])} refs (ratio undefined, not plotted)" for v in D.VARIANTS if nz[v]),
           transform=ax[2].transAxes, fontsize=7, va="top", color=INK2)
ax[0].legend(); fig.suptitle("Phase-specific plasticity at the macroscopic yield point (72 independent references per material)", fontsize=10); save(fig, "phase_diagnostics.png")
# 10 calibration
fig, ax = plt.subplots(1, 3, figsize=(12, 3.4))
for i, v in enumerate(D.VARIANTS):
    for m in COL:
        g = A[(A.variant == v) & (A.surrogate == "A20") & (A.method == m)].groupby("n")
        ax[0].plot(g.cov95.median().index, g.cov95.median(), color=COL[m], ls="-" if v == "B1" else "--", marker="o", ms=3, label=f"{v} {NAME[m]}")
        ax[1].plot(g.cov90.median().index, g.cov90.median(), color=COL[m], ls="-" if v == "B1" else "--", marker="o", ms=3)
        ax[2].plot(g.rms_z.median().index, g.rms_z.median(), color=COL[m], ls="-" if v == "B1" else "--", marker="o", ms=3)
ax[0].axhline(0.95, color=INK2, lw=0.8); ax[1].axhline(0.90, color=INK2, lw=0.8); ax[2].axhline(1, color=INK2, lw=0.8); ax[2].set_yscale("log")
for a, t in zip(ax, ("95% interval coverage", "90% interval coverage", "RMS standardised residual")):
    a.set_title(t + " (A20, median over designs)", fontsize=9); a.set_xlabel("n")
ax[0].legend(fontsize=6.5, ncol=2); save(fig, "calibration.png")
print("figures:", sorted(p.name for p in FIG.glob("*.png")))
