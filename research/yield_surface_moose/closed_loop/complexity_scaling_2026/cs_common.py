"""Complexity-scaling study 2026: shared definitions (candidate library, analytic ground-truth family, complexity descriptors).

SEPARATE from the frozen B1/B2 benchmark.  Imports the frozen V1 package and closed_loop/v2/bomes_v2.py READ-ONLY; nothing frozen is modified.

Candidate library  : 4,096 unoriented loading axes = the strict upper half (z3 > 0) of a 8,192-point Fibonacci sphere in Mandel space
                     z = (s_xx, s_yy, sqrt2 s_xy)/|.|.  Fibonacci z-levels are symmetric about 0 and none is 0 for an even count, so exactly
                     4,096 points have z3 > 0, no axis lies on the rim and no two members are antipodal (checked).  Uniform in area on the
                     projective plane RP^2 (each axis {d, -d} represented once).
Reference set      : 20,000-point Fibonacci sphere (both hemispheres, equal-area weights), disjoint from the library (checked: min angle > 0).
Analytic family    : r(d) = R0 / g(d) with the convex, even, positively homogeneous gauge
                         g(z) = ( q(z)^p + lam * sum_k (w_k |a_k . z|)^p )^(1/p),   q(z) = sqrt(z^T M z)
                     g is a norm (p-norm of norms), so every surface is convex and exactly antipodally symmetric; for p even it is C-infinity.
                     Larger p and lam concentrate curvature at the ridges between the q-dominated and the a_k-dominated regions.
"""
import hashlib, json, os, sys
from pathlib import Path
import numpy as np

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent; CL = HERE.parent
REL = Path.home() / "projects/yield_surface_discovery_release/software/yield_surface_discovery/src"
for p in (str(CL / "v2"), str(CL), str(REL)):
    if p not in sys.path:
        sys.path.insert(0, p)
from yield_surface_discovery.candidates import CandidateLibrary          # noqa: E402  frozen V1
SQ2 = np.sqrt(2.0)
M_J2 = np.array([[1.0, -0.5, 0.0], [-0.5, 1.0, 0.0], [0.0, 0.0, 1.5]])     # plane-stress von Mises in Mandel form (homog_benchmark/analytic_reference.py)
BUDGETS = (8, 16, 33, 64, 128, 256, 512, 1000)


def fibonacci_sphere(n):                                                  # identical to homog_benchmark/analytic_reference.fibonacci_sphere
    i = np.arange(n) + 0.5; z = 1 - 2 * i / n; rr = np.sqrt(1 - z ** 2); ph = np.pi * (1 + 5 ** 0.5) * i
    return np.c_[rr * np.cos(ph), rr * np.sin(ph), z]


def library(n_axes=4096):
    F = fibonacci_sphere(2 * n_axes); A = F[F[:, 2] > 0]; assert len(A) == n_axes
    A = A[:, [1, 2, 0]] if False else A                                   # keep (z1, z2, z3) = (s_xx, s_yy, sqrt2 s_xy) ordering as generated
    lib = CandidateLibrary(A.copy(), np.arange(n_axes), [f"fa{i:04d}" for i in range(n_axes)], ["fibonacci_axis"] * n_axes, f"fibonacci_axes_{n_axes}")
    return lib


def library_hash(lib):
    return hashlib.sha256(np.ascontiguousarray(np.round(lib.directions, 15)).tobytes()).hexdigest()


def references(n=20000):
    return fibonacci_sphere(n)


# ------------------------------------------------------------------ analytic ground-truth family
def _axes(kind, seed=0, k=10):
    if kind == "ico6":                                                    # the 6 axes of the icosahedron (vertices up to sign)
        t = (1 + 5 ** 0.5) / 2; V = np.array([[0, 1, t], [0, -1, t], [1, t, 0], [-1, t, 0], [t, 0, 1], [-t, 0, 1]], float)
        return V / np.linalg.norm(V, axis=1)[:, None], np.ones(6)
    rng = np.random.default_rng(seed); V = rng.normal(size=(k, 3)); V /= np.linalg.norm(V, axis=1)[:, None]
    return V, rng.uniform(0.85, 1.15, k)


M_HILL = np.array([[1.0, -0.35, 0.0], [-0.35, 0.62, 0.0], [0.0, 0.0, 1.9]])   # an anisotropic, positive-definite quadratic (Hill-48 type)
LEVELS = {   # name: (M, p, lam, axes kind, seed, k, R0)   -- R0 fixed below so that the median radius is ~200 MPa
    "A1_J2":         (M_J2,   2,  0.0, None,   0, 0),
    "A2_hill":       (M_HILL, 2,  0.0, None,   0, 0),
    "A3_ico_p6":     (M_J2,   6,  0.6, "ico6", 0, 6),
    "A4_ico_p16":    (M_J2,   16, 1.0, "ico6", 0, 6),
    "A5_rand10_p32": (M_HILL, 32, 1.0, "rand", 11, 10),
    "A6_rand16_p64": (M_HILL, 64, 1.0, "rand", 23, 16),
}
# PREDEFINED complexity classes (frozen in PROTOCOL.md before any comparison), applied identically to analytic and MOOSE surfaces:
# primary descriptor = 99th percentile of the dimensionless maximum principal curvature kappa_1 * r over the equal-area reference sphere
CLASSES = ((3.0, "C1_smooth"), (6.0, "C2_moderate"), (12.0, "C3_localized"), (np.inf, "C4_highly_localized"))


def complexity_class(kappa_p99):
    return next(lab for th, lab in CLASSES if kappa_p99 < th)


def gauge(level, Z):
    M, p, lam, kind, seed, k = LEVELS[level]; Z = np.atleast_2d(Z)
    q = np.sqrt(np.einsum("ij,jk,ik->i", Z, M, Z))
    if lam == 0:
        return q
    A, w = _axes(kind, seed, k); t = np.abs(Z @ A.T) * w[None, :]
    m = np.maximum(q, t.max(1))                                           # scale for a stable p-norm
    return m * ((q / m) ** p + lam * np.sum((t / m[:, None]) ** p, axis=1)) ** (1.0 / p)


_R0 = {}


def radius_fn(level):
    if level not in _R0:
        _R0[level] = 200.0 * float(np.median(gauge(level, fibonacci_sphere(20000))))
    R0 = _R0[level]
    return lambda D: R0 / gauge(level, np.atleast_2d(D))


# ------------------------------------------------------------------ complexity descriptors (computed BEFORE any comparison)
def tangent_basis(D):
    a = np.where(np.abs(D[:, [2]]) < 0.9, np.array([[0, 0, 1.0]]), np.array([[1.0, 0, 0]]))
    e1 = np.cross(D, a); e1 /= np.linalg.norm(e1, axis=1)[:, None]; e2 = np.cross(D, e1); return e1, e2


def curvatures(rfun, D, h=2e-3):
    """principal curvatures of the star-shaped surface X(d) = r(d) d by second-order finite differences in a local chart
    c(u, v) = normalize(d + u e1 + v e2).  Returns (k1 >= k2) with the outward normal (convex => both >= 0)."""
    e1, e2 = tangent_basis(D)
    def X(u, v):
        c = D + u * e1 + v * e2; c /= np.linalg.norm(c, axis=1)[:, None]; return rfun(c)[:, None] * c
    X0 = X(0, 0); Xu = (X(h, 0) - X(-h, 0)) / (2 * h); Xv = (X(0, h) - X(0, -h)) / (2 * h)
    Xuu = (X(h, 0) - 2 * X0 + X(-h, 0)) / h ** 2; Xvv = (X(0, h) - 2 * X0 + X(0, -h)) / h ** 2
    Xuv = (X(h, h) - X(h, -h) - X(-h, h) + X(-h, -h)) / (4 * h ** 2)
    N = np.cross(Xu, Xv); N /= np.linalg.norm(N, axis=1)[:, None]; N *= np.sign(np.sum(N * X0, axis=1))[:, None]
    E, F, G = (Xu * Xu).sum(1), (Xu * Xv).sum(1), (Xv * Xv).sum(1)
    L, Mm, Nn = -(Xuu * N).sum(1), -(Xuv * N).sum(1), -(Xvv * N).sum(1)     # sign: convex surface -> positive curvature w.r.t. outward normal
    den = E * G - F ** 2; K = (L * Nn - Mm ** 2) / den; H = (E * Nn + G * L - 2 * F * Mm) / (2 * den)
    s = np.sqrt(np.maximum(H ** 2 - K, 0)); return H + s, H - s


def descriptors(rfun, n=20000, seed=0):
    D = fibonacci_sphere(n); r = rfun(D); lr = np.log(r); rm = float(np.median(r))
    k1, k2 = curvatures(rfun, D); kn = k1 * r                                # dimensionless max principal curvature (sphere = 1)
    rng = np.random.default_rng(seed); i = rng.integers(0, n, 400000); j = rng.integers(0, n, 400000)
    ang = np.degrees(np.arccos(np.clip(np.abs(np.sum(D[i] * D[j], 1)), -1, 1)))     # AXIAL separation (projective)
    x = lr - lr.mean(); c = x[i] * x[j] / x.var(); bins = np.arange(0, 91, 2.5); corr = [c[(ang >= a) & (ang < a + 2.5)].mean() for a in bins[:-1]]
    below = np.where(np.array(corr) < np.exp(-1))[0]; ell = float(bins[below[0]] if len(below) else 90.0)
    srt = np.sort(kn)[::-1]; top = srt[: n // 10].sum() / srt.sum()
    anti = float(np.max(np.abs(np.log(rfun(-D)) - lr)))
    return dict(anisotropy_max_over_min=float(r.max() / r.min()), std_log_r=float(lr.std()), median_r=rm,
                kappa_max_dimless_p99=float(np.percentile(kn, 99)), kappa_max_dimless_max=float(kn.max()), kappa_cv=float(kn.std() / kn.mean()),
                curvature_top10pct_share=float(top), min_principal_curvature_dimless=float((k2 * r).min()), convex=bool((k2 * r).min() > -1e-3),
                angular_corr_length_deg=ell, antipodal_asym_max=anti, corr_curve=[float(v) for v in corr])


def file_hashes(paths):
    return {str(p): hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in paths}
