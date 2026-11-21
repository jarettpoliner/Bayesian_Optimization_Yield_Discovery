"""Analytical yield-radius models used as reproducible oracles (the controller never sees these formulas).

cb2004_mg4li: Cazacu-Barlat 2004 plane-stress yield function for Mg-4%Li as used by Fuhg, Fau, Bouklas & Marino (2023), the benchmark target:
  f = (J2^0)^{3/2} - c J3^0 - sbar^3 with a1 = 0.869, a2 = 3.371, a3 = 3.509, a4 = 1, b1 = -1.591, b2 = 5.414, b3 = 3.957, b4 = 0.259, b5 = b10 = 1,
  c = 2.01, sbar = 300 MPa (homogeneous of degree 3 -> r(d) = sbar / g(d)^{1/3}).
von_mises: plane-stress von Mises with yield stress sigma_y."""
import numpy as np
from .geometry import mandel_to_cauchy

CB2004_MG4LI = dict(a1=0.869, a2=3.371, a3=3.509, a4=1.0, b1=-1.591, b2=5.414, b3=3.957, b4=0.259, b5=1.0, b10=1.0, c=2.01, sbar=300.0)


def _g_cb(S, p):
    x, y, t = S[..., 0], S[..., 1], S[..., 2]
    J2 = p["a1"] / 6 * (x - y) ** 2 + p["a2"] / 6 * y ** 2 + p["a3"] / 6 * x ** 2 + p["a4"] * t ** 2
    J3 = ((p["b1"] + p["b2"]) / 27 * x ** 3 + (p["b3"] + p["b4"]) / 27 * y ** 3 - (1 / 9) * p["b1"] * y * x ** 2 - (1 / 9) * p["b4"] * x * y ** 2
          + (t ** 2 / 3) * (p["b5"] * y + (2 * p["b10"] - p["b5"]) * x))
    return np.maximum(J2, 0) ** 1.5 - p["c"] * J3


def cb2004_radius(D, p=CB2004_MG4LI):
    g = _g_cb(mandel_to_cauchy(D), p)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(g > 0, p["sbar"] / np.cbrt(np.maximum(g, 1e-300)), np.inf)


def von_mises_radius(D, sigma_y=1.0):
    s = mandel_to_cauchy(D); x, y, t = s[..., 0], s[..., 1], s[..., 2]
    return sigma_y / np.sqrt(x * x + y * y - x * y + 3 * t * t)


MODELS = {"cb2004_mg4li": lambda D, **kw: cb2004_radius(D), "von_mises": lambda D, sigma_y=1.0, **kw: von_mises_radius(D, sigma_y)}
