#!/usr/bin/env python3
"""Holdout population H1 (independent of everything used to design or evaluate V3): Barlat Yld91 plane-stress yield surfaces.

Family (Barlat, Lege & Brem 1991; Hosford-type isotropic function of a linearly transformed stress), plane stress, sigma_zz = 0:
  material frame   sigma_m = Q(theta) sigma Q(theta)^T                         (random in-plane material orientation theta)
  transformed      s'_xx = [c (s_xx - s_yy) + b s_xx] / 3,  s'_yy = [a s_yy - c (s_xx - s_yy)] / 3,  s'_zz = -(b s_xx + a s_yy) / 3,  s'_xy = h s_xy
  principal        S1,2 = in-plane eigenvalues of [[s'_xx, s'_xy], [s'_xy, s'_yy]],  S3 = s'_zz
  gauge            g(sigma) = ( (|S1 - S2|^m + |S2 - S3|^m + |S3 - S1|^m) / 2 )^(1/m),   r(d) = R0 / g(d)     (d: Mandel direction)
Admissibility: for m >= 1 the Hosford function is a convex symmetric function of the principal values, hence a convex spectral function, and
the map sigma -> s' is linear, so g is convex and positively homogeneous; g is even (antipodal symmetry exact).  a = b = c = h = 1, m = 2 is
von Mises.  Verified numerically below (convexity, symmetry, positivity, curvature) before any controller is run.

Population (fixed here, before any controller run; parameters chosen from physical convention, never from algorithm performance):
  directional complexity  m in {2 (quadratic, Hill-type), 6 (BCC convention), 8 (FCC convention)}
  anisotropy level        log a, log b, log c, log h ~ N(0, s^2) iid,  s in {0.15 (mild), 0.30 (strong)}
  6 levels x 4 realisations = 24.  Realisation i: level (m, s) = LEVELS[i // 4], seed 52000 + i (coefficients and theta ~ U[0, pi)),
  controller start index 2000 + i (fresh starts; never used before).  Replacement rule (geometry only): a realisation failing any check is
  replaced by the next unused seed 52024, 52025, ... of the same level.
R0: median radius over the references = 200.
References: 20,000-point Fibonacci sphere rotated by scipy Rotation.random(random_state=52999) -- a point set never used before; equal-area;
  its minimum axial angle to the 4,096-axis library is checked > 0.
    python surfaces_h1.py  ->  population.json, references_h1.npy, surface_verification.json"""
import hashlib, json, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; CL = HERE.parents[1]; CS = CL / "complexity_scaling_2026"; sys.dont_write_bytecode = True
sys.path.insert(0, str(CS))
import cs_common as C                                                                # noqa: E402  (library, Fibonacci sphere, curvature routine; frozen)
LEVELS = [(2, 0.15), (2, 0.30), (6, 0.15), (6, 0.30), (8, 0.15), (8, 0.30)]; N_REAL = 24; SEED0 = 52000; START0 = 2000; SQ2 = np.sqrt(2.0)


def params(seed, s):
    rng = np.random.default_rng(seed); a, b, c, h = np.exp(s * rng.normal(size=4)); th = float(rng.uniform(0, np.pi)); return dict(a=a, b=b, c=c, h=h, theta=th)


def gauge(Z, p, m):
    Z = np.atleast_2d(Z); sx, sy, txy = Z[:, 0], Z[:, 1], Z[:, 2] / SQ2; t = p["theta"]; co, si = np.cos(t), np.sin(t)
    mx = co * co * sx + si * si * sy + 2 * co * si * txy; my = si * si * sx + co * co * sy - 2 * co * si * txy; mt = (co * co - si * si) * txy + co * si * (sy - sx)
    a, b, c, h = p["a"], p["b"], p["c"], p["h"]
    xx = (c * (mx - my) + b * mx) / 3; yy = (a * my - c * (mx - my)) / 3; zz = -(b * mx + a * my) / 3; xy = h * mt
    cen, rad = (xx + yy) / 2, np.sqrt(((xx - yy) / 2) ** 2 + xy ** 2); S1, S2, S3 = cen + rad, cen - rad, zz
    D = np.c_[np.abs(S1 - S2), np.abs(S2 - S3), np.abs(S3 - S1)]; sc = D.max(1)
    return sc * (np.sum((D / sc[:, None]) ** m, axis=1) / 2) ** (1.0 / m)


class Yld91:
    def __init__(self, rec):
        self.p = {k: rec[k] for k in ("a", "b", "c", "h", "theta")}; self.m = rec["m"]; self.R0 = rec["R0"]

    def radius(self, D):
        return self.R0 / gauge(D, self.p, self.m)


def references():
    from scipy.spatial.transform import Rotation
    return C.fibonacci_sphere(20000) @ Rotation.random(random_state=52999).as_matrix().T


def verify(rec, R, rng):
    F = Yld91(rec); g = lambda Z: gauge(Z, F.p, F.m); r = F.radius(R)
    X = rng.normal(size=(200000, 3)); Y = rng.normal(size=(200000, 3)); conv = float(np.max((g((X + Y) / 2) - (g(X) + g(Y)) / 2) / ((g(X) + g(Y)) / 2)))
    anti = float(np.max(np.abs(g(-R) / g(R) - 1))); k1, k2 = C.curvatures(F.radius, R); kn = k1 * r
    return dict(convexity_max_rel_violation=conv, antipodal_max_rel=anti, positive=bool(np.all(r > 0) and np.all(np.isfinite(r))),
                min_principal_curvature_dimless=float((k2 * r).min()), kappa1_dimless_p99=float(np.percentile(kn, 99)), kappa1_dimless_max=float(kn.max()),
                anisotropy_max_over_min=float(r.max() / r.min()), std_log_r=float(np.log(r).std()),
                ok=bool(conv <= 1e-12 and anti <= 1e-12 and np.all(r > 0) and (k2 * r).min() >= -1e-3))


if __name__ == "__main__":
    R = references(); lib = C.library().directions; assert C.library_hash(C.library()).startswith("8e1654f8")
    minang = float(np.degrees(np.arccos(np.clip(np.abs(R @ lib.T).max(), 0, 1)))); assert minang > 0
    pop, ver, nxt, rng = [], [], N_REAL, np.random.default_rng(1)
    for i in range(N_REAL):
        m, s = LEVELS[i // 4]; seed = SEED0 + i
        while True:
            p = params(seed, s); rec = dict(realisation=i, level=f"m{m}_s{s:.2f}", m=m, s=s, seed=seed, start=START0 + i, **p, R0=1.0)
            rec["R0"] = float(200.0 * np.median(gauge(R, p, m))); v = verify(rec, R, rng)
            ver.append(dict(realisation=i, seed=seed, **v))
            if v["ok"]:
                break
            seed = SEED0 + nxt; nxt += 1
        pop.append(rec)
    np.save(HERE / "references_h1.npy", R)
    json.dump(dict(family="Barlat Yld91 plane stress", levels=LEVELS, realisations=pop, references_sha256=hashlib.sha256(np.ascontiguousarray(R).tobytes()).hexdigest(),
                   references_min_axial_angle_to_library_deg=minang, library_hash=C.library_hash(C.library())), open(HERE / "population.json", "w"), indent=1)
    json.dump(dict(checks=ver, all_accepted_first_try=len(ver) == N_REAL), open(HERE / "surface_verification.json", "w"), indent=1)
    for v, p in zip([x for x in ver if x["ok"]], pop):
        print(p["realisation"], p["level"], f"aniso {v['anisotropy_max_over_min']:.2f} k99 {v['kappa1_dimless_p99']:.2f} kmin {v['min_principal_curvature_dimless']:.3f} conv {v['convexity_max_rel_violation']:.1e}")
    print("references min axial angle to library (deg):", round(minang, 4), "| replacements:", len(ver) - N_REAL)
