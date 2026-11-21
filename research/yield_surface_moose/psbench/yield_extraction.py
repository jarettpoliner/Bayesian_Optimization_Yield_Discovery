"""Effective-yield definitions evaluated on a stress-controlled proportional loading curve.

All curves are parameterized by the Mandel stress radius r (macroscopic stress
sigma_M = r * d with |d| = 1).
"""

import numpy as np

# Candidate thresholds on the matrix-volume-averaged accumulated equivalent plastic strain.
EQPS_THRESHOLDS = (1e-5, 3e-5, 1e-4, 3e-4, 1e-3)
# Primary effective-yield definition, selected from validation/composite_study.py (mesh mode):
# on S1-S3 over 6 directions its radius differs from h = 1/80 by <= 1.4 % at h = 1/40
# and 1/20, versus 2.4-11 % for first local (quadrature-point) yield.
PRIMARY_EQPS = 1e-4
# Offset for the independent 0.2 % offset diagnostic (equivalent macroscopic plastic strain).
OFFSET = 2e-3


def first_crossing(r, q, level):
    """Radius at which the monotone-in-practice sequence q(r) first reaches `level`,
    by linear interpolation inside the bracketing load step. Returns (r_y, i_lo) or
    (nan, None) if never reached."""
    r, q = np.asarray(r, float), np.asarray(q, float)
    idx = np.nonzero(q >= level)[0]
    if idx.size == 0:
        return float("nan"), None
    k = int(idx[0])
    if k == 0:
        return float(r[0]), 0
    w = (level - q[k - 1]) / (q[k] - q[k - 1])
    return float(r[k - 1] + w * (r[k] - r[k - 1])), k - 1


def macro_plastic_strain_eq(r, eps, compliance):
    """Equivalent macroscopic plastic strain sqrt(2/3 eps_p : eps_p).

    eps: (n, 4) macroscopic total strain tensor components [xx, yy, xy, zz] (tensor shear).
    compliance: macroscopic strain per unit radius in the elastic regime (same layout),
    measured from the elastic probe step. Under proportional stress loading the elastic
    macroscopic strain is exactly r * compliance, so eps_p = eps - r * compliance.
    """
    ep = np.asarray(eps, float) - np.outer(r, compliance)
    sq = ep[:, 0] ** 2 + ep[:, 1] ** 2 + ep[:, 3] ** 2 + 2.0 * ep[:, 2] ** 2
    return np.sqrt(2.0 / 3.0 * sq)


def interp_at(r, y, r_star):
    """Linear interpolation of curve quantity y (n,) or (n, m) at radius r_star."""
    y = np.asarray(y, float)
    if not np.isfinite(r_star):
        return np.full(y.shape[1:], np.nan) if y.ndim > 1 else float("nan")
    if y.ndim == 1:
        return float(np.interp(r_star, r, y))
    return np.array([np.interp(r_star, r, y[:, j]) for j in range(y.shape[1])])
