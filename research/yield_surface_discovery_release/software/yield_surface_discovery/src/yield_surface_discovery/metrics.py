"""Evaluation metrics (frozen definitions of the benchmark), available whenever reference radii are known on evaluation directions.

global p95 relative error = 95th percentile of |r_hat / r_ref - 1| over the held-out directions (numpy linear quantile).
weak10 / sharp10 regions: weakest 10% of reference radii / 10% of directions with the largest discrete Laplacian of log r (9 nearest neighbours).
max overprediction = largest signed relative error in a region.  High shear: |z3|/|z| in [0.85, 1].  Pool-1 envelope: |z3|/|z| > 1/sqrt2."""
import numpy as np
from .geometry import angle_rad

SHEAR_BANDS = [(0, 0.25), (0.25, 0.5), (0.5, 1 / np.sqrt(2)), (1 / np.sqrt(2), 0.85), (0.85, 1.0 + 1e-9)]
BAND_NAMES = ["0-0.25", "0.25-0.5", "0.5-0.707", "0.707-0.85", "0.85-1"]


def region_masks(r_ref, directions, frac=0.10):
    from scipy.spatial import cKDTree
    _, nb = cKDTree(directions).query(directions, 9); lr = np.log(r_ref)
    sharp = np.abs(lr - lr[nb[:, 1:]].mean(1))
    return dict(sharp10=sharp >= np.quantile(sharp, 1 - frac), weak10=r_ref <= np.quantile(r_ref, frac))


def evaluate(r_hat, r_ref, directions, masks=None):
    r_hat, r_ref = np.asarray(r_hat, float), np.asarray(r_ref, float); rel = r_hat / r_ref - 1; a = np.abs(rel)
    masks = masks or region_masks(r_ref, directions); sh = np.abs(np.asarray(directions)[:, 2])
    out = dict(global_p95=float(np.quantile(a, .95)), global_median=float(np.quantile(a, .5)), global_max=float(a.max()),
               signed_p95=float(np.quantile(rel, .95)), signed_p05=float(np.quantile(rel, .05)))
    for reg in ("weak10", "sharp10"):
        m = masks[reg]; out[f"{reg}_signed_p95"] = float(np.quantile(rel[m], .95)); out[f"{reg}_max_overprediction"] = float(rel[m].max())
    for (lo, hi), nm in zip(SHEAR_BANDS, BAND_NAMES):
        m = (sh >= lo) & (sh < hi); out[f"band_{nm}_p95"] = float(np.quantile(a[m], .95))
    out["high_shear_p95"] = out["band_0.85-1_p95"]
    m = sh > 1 / np.sqrt(2) + 1e-9; out["beyond_conventional_envelope_p95"] = float(np.quantile(a[m], .95)) if m.any() else float("nan")
    return out


def coverage_gaps(eval_directions, tested_directions, symmetry=None):
    """angular distance (deg) from every evaluation direction to the nearest tested direction (or its symmetry image)."""
    C = symmetry.covered(tested_directions) if symmetry is not None else np.atleast_2d(tested_directions)
    g = np.degrees(angle_rad(eval_directions, C).min(1))
    return dict(gap_median_deg=float(np.median(g)), gap_p90_deg=float(np.quantile(g, .9)), gap_max_deg=float(g.max()))
