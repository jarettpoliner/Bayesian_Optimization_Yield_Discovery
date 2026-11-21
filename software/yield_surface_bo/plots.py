"""Campaign figures: GP reconstruction (labelled prediction, 95 % band) on two great circles with the measured observations near each circle,
the analytic truth or held-out references where available, and the error curve.  PNG (200 dpi) + SVG."""
from pathlib import Path
import numpy as np


def plot_campaign(out, D, r, gp, truth, err_rows, cfg, n_init):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    from .campaign import predict
    col = {"v2_mes": "#2a78d6", "v3_weighted": "#1f9e74", "v3_signed": "#e87ba4", "angular_maximin": "#eb6834", "random": "#8a8984"}.get(cfg["method"], "#2a78d6")
    circles = (("σxx–σyy plane (σxy = 0)", np.array([1.0, 0, 0]), np.array([0, 1.0, 0])),
               ("equal-biaxial → shear circle", np.array([1, 1, 0]) / np.sqrt(2), np.array([0, 0, 1.0])))
    ncol = 3 if err_rows else 2; fig, axs = plt.subplots(1, ncol, figsize=(5.2 * ncol, 4.2))
    t = np.radians(np.arange(0, 360.5, 1.0))
    ref = None
    if cfg.get("evaluation", {}).get("references_csv"):
        import pandas as pd
        ref = pd.read_csv(cfg["evaluation"]["references_csv"])
    for ax, (lab, e1, e2) in zip(axs, circles):
        G = np.cos(t)[:, None] * e1 + np.sin(t)[:, None] * e2; rp, sd = predict(gp, G); nrm = np.cross(e1, e2)
        ax.fill_between(np.degrees(t), rp * np.exp(-1.96 * sd), rp * np.exp(1.96 * sd), color=col, alpha=0.15, lw=0, label="GP 95 % band (prediction)")
        ax.plot(np.degrees(t), rp, color=col, lw=1.6, label="GP mean (prediction)")
        if truth is not None:
            ax.plot(np.degrees(t), truth(G), color="#0b0b0b", lw=1.0, ls="--", label="analytic truth")
        if ref is not None:
            P = ref[["z1", "z2", "z3"]].to_numpy(); near = np.abs(P @ nrm) < np.sin(np.radians(1.0))
            if near.any():
                q = P[near]; ax.scatter(np.degrees(np.arctan2(q @ e2, q @ e1)) % 360, ref.yield_radius.to_numpy()[near], s=6, color="#0b0b0b", zorder=3, label="held-out references")
        near = np.abs(D @ nrm) <= np.sin(np.radians(4.0))
        if near.any():
            q = D[near]; ax.scatter(np.degrees(np.arctan2(q @ e2, q @ e1)) % 360, r[near], s=28, marker="D", color=col, edgecolor="white", lw=0.6, zorder=4,
                                    label="observations within 4°")
        ax.set_title(lab, fontsize=9, loc="left"); ax.set_xlabel("angle along the circle [deg]"); ax.set_ylabel("yield radius r"); ax.set_xticks(range(0, 361, 90))
        ax.grid(alpha=0.3)
    axs[0].legend(fontsize=7, frameon=False)
    if err_rows:
        ax = axs[-1]; n = [x["n"] for x in err_rows]; e = [100 * x["mean_rel_error"] for x in err_rows]
        ax.plot(n, e, "-o", color=col); ax.set_xlabel("observations n"); ax.set_ylabel("mean |r̂/r − 1| [%]"); ax.set_yscale("log"); ax.grid(alpha=0.3)
        ax.set_title("reconstruction error vs n (evaluation only)", fontsize=9, loc="left")
    fig.suptitle(f"{cfg['method']}: {len(r)} observations ({n_init} initial)", x=0.01, ha="left", fontsize=10)
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(Path(out) / f"reconstruction.{ext}", dpi=200)
    plt.close(fig)
