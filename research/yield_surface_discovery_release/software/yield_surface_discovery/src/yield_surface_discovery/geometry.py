"""Stress-direction geometry.  Plane stress, Mandel coordinates z = [s_xx, s_yy, sqrt(2) s_xy]; a loading direction is a unit z-vector d and the
yield radius r(d) is the Mandel norm of the stress at first yield along d (so the yield stress is r * d)."""
import numpy as np

SQ2 = np.sqrt(2.0)


def mandel_to_cauchy(z):
    """z = [s_xx, s_yy, sqrt2 s_xy] -> [s_xx, s_yy, s_xy] (tensor shear)."""
    return np.asarray(z, float) * np.array([1.0, 1.0, 1.0 / SQ2])


def cauchy_to_mandel(s):
    return np.asarray(s, float) * np.array([1.0, 1.0, SQ2])


def unit(z, tol=1e-12):
    """normalise to unit length; vectors already unit within tol are returned bit-for-bit unchanged (keeps library directions exact)."""
    z = np.asarray(z, float); n = np.linalg.norm(z, axis=-1, keepdims=True)
    return np.where(np.abs(n - 1.0) <= tol, z, z / n)


def angle_rad(A, B):
    """accurate pairwise angle between unit vectors (atan2 form; well conditioned near 0 and pi)."""
    A = np.atleast_2d(A); B = np.atleast_2d(B)
    cr = np.linalg.norm(np.cross(A[:, None, :], B[None, :, :]), axis=2); return np.arctan2(cr, A @ B.T)


def arccos_angle_rad(A, B):
    """the arccos form used by the frozen angular-maximin implementation (kept for exact reproduction)."""
    return np.arccos(np.clip(np.atleast_2d(A) @ np.atleast_2d(B).T, -1, 1))


def icosphere(level=4):
    """unit icosphere (vertices, faces) by repeated midpoint subdivision; pure numpy."""
    t = (1 + 5 ** 0.5) / 2
    V = [(-1, t, 0), (1, t, 0), (-1, -t, 0), (1, -t, 0), (0, -1, t), (0, 1, t), (0, -1, -t), (0, 1, -t), (t, 0, -1), (t, 0, 1), (-t, 0, -1), (-t, 0, 1)]
    F = [(0, 11, 5), (0, 5, 1), (0, 1, 7), (0, 7, 10), (0, 10, 11), (1, 5, 9), (5, 11, 4), (11, 10, 2), (10, 7, 6), (7, 1, 8),
         (3, 9, 4), (3, 4, 2), (3, 2, 6), (3, 6, 8), (3, 8, 9), (4, 9, 5), (2, 4, 11), (6, 2, 10), (8, 6, 7), (9, 8, 1)]
    V = [np.array(v, float) / np.linalg.norm(v) for v in V]
    for _ in range(level):
        cache, F2 = {}, []

        def mid(a, b):
            k = (min(a, b), max(a, b))
            if k not in cache:
                m = V[a] + V[b]; V.append(m / np.linalg.norm(m)); cache[k] = len(V) - 1
            return cache[k]
        for a, b, c in F:
            ab, bc, ca = mid(a, b), mid(b, c), mid(c, a); F2 += [(a, ab, ca), (b, bc, ab), (c, ca, bc), (ab, bc, ca)]
        F = F2
    return np.array(V), np.array(F, int)


def lonlat_deg(D):
    """(in-plane angle atan2(z2, z1), shear latitude asin(z3)) in degrees; the shear mirror flips the latitude."""
    D = np.atleast_2d(D); return np.degrees(np.arctan2(D[:, 1], D[:, 0])), np.degrees(np.arcsin(np.clip(D[:, 2], -1, 1)))
