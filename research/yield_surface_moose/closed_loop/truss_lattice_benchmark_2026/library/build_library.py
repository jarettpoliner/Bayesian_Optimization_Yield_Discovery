"""Versioned candidate library TLIB-4096 v1 and independent reference design TREF v1 (pure geometry; frozen before any reference MOOSE run).

TLIB-4096 v1 (V1 CandidateLibrary CSV format: id,label,family,z1,z2,z3; unit Mandel direction z = (s_xx, s_yy, sqrt2 s_xy))
  * 4096 axes: hemispherical Fibonacci set on {z3 > 0}: z3_k = 1 - (k + 0.5)/4096, azimuth = k * golden angle, k = 0..4095, then the
    antipode of every axis.  ids 0..4095 = axis representatives (the lowest index, as the frozen AxisLibrary picks), ids 4096..8191 = antipodes.
  * -> 8192 oriented directions, no duplicates, antipode map exact (used only if antipodal symmetry is admissible; it is NOT for TL1, see
    TRUSS_ORACLE_VALIDATION.md, so the campaign candidate set is all 8192 oriented directions with symmetry='none').
  * angular metric: great-circle angle between unit Mandel vectors (frozen geometry.angle_rad); axial angle min(t, pi - t) only in antipodal mode.
  The original 207-axis library is not touched.

TREF v1 (evaluator side only; never read by a controller)
  * AREA set: 400-point full-sphere Fibonacci set rotated by a fixed random rotation (seed 4242) -> area-uniform, so plain means are
    area-weighted; every point must be >= 0.5 deg from every library direction (exact-match reuse impossible); points failing are dropped.
  * SWEEP_EB_PS: great circle cos t (1,1,0)/sqrt2 + sin t (0,0,1), t = phi + 5 k deg, k = 0..71 (the B1/B2 continuity-sweep circle, full 360 deg).
  * SWEEP_BIAX: great circle cos t (1,0,0) + sin t (0,1,0) (sigma_xx - sigma_yy plane), t = phi + 5 k deg, k = 0..71.
    phi = smallest value in {0.50, 0.55, ..., 4.95} deg such that every sweep point is >= 0.25 deg from the library (recorded; 0.5 deg is not
    attainable for 72 points on one circle with a 2.2-deg library spacing). Distinct from every library direction -> no exact-match reuse.
"""
import hashlib, json
from pathlib import Path
import numpy as np, pandas as pd

OUT = Path(__file__).resolve().parent
N_AX = 4096
GOLD = np.pi * (3 - np.sqrt(5))


def library():
    k = np.arange(N_AX); z = 1 - (k + 0.5) / N_AX; rho = np.sqrt(1 - z ** 2); ph = k * GOLD
    A = np.column_stack([rho * np.cos(ph), rho * np.sin(ph), z]); D = np.vstack([A, -A])
    ids = np.arange(2 * N_AX); lab = [f"T{i:04d}+" for i in range(N_AX)] + [f"T{i:04d}-" for i in range(N_AX)]
    fam = ["axis_rep"] * N_AX + ["antipode"] * N_AX
    return pd.DataFrame(dict(id=ids, label=lab, family=fam, z1=D[:, 0], z2=D[:, 1], z3=D[:, 2]))


def ang_deg(P, Q):
    P = np.atleast_2d(P); Q = np.atleast_2d(Q); cr = np.linalg.norm(np.cross(P[:, None, :], Q[None, :, :]), axis=2)
    return np.degrees(np.arctan2(cr, P @ Q.T))


def fib_sphere(n):
    k = np.arange(n); z = 1 - (2 * k + 1) / n; rho = np.sqrt(1 - z ** 2); ph = k * GOLD; return np.column_stack([rho * np.cos(ph), rho * np.sin(ph), z])


def references(D):
    from scipy.spatial.transform import Rotation
    rot = Rotation.random(random_state=4242).as_matrix(); P = fib_sphere(400) @ rot.T
    dmin = ang_deg(P, D).min(1); keep = dmin >= 0.5
    rows = [dict(set="area", z1=p[0], z2=p[1], z3=p[2], min_deg_to_library=float(m)) for p, m in zip(P[keep], dmin[keep])]
    info = dict(area_generated=400, area_dropped=int((~keep).sum()))
    for name, e, s in (("sweep_eb_ps", np.array([1, 1, 0]) / np.sqrt(2), np.array([0, 0, 1.0])), ("sweep_biax", np.array([1.0, 0, 0]), np.array([0, 1.0, 0]))):
        for phi in np.arange(0.5, 5.0, 0.05):
            t = np.radians(phi + 5 * np.arange(72)); Q = np.cos(t)[:, None] * e + np.sin(t)[:, None] * s; m = ang_deg(Q, D).min(1)
            if m.min() >= 0.25:
                break
        else:
            raise SystemExit(f"no admissible phase for {name}")
        info[f"{name}_phase_deg"] = float(phi); info[f"{name}_min_deg_to_library"] = float(m.min())
        rows += [dict(set=name, z1=q[0], z2=q[1], z3=q[2], min_deg_to_library=float(mm), sweep_t_deg=float(phi + 5 * k)) for k, (q, mm) in enumerate(zip(Q, m))]
    R = pd.DataFrame(rows); R.insert(0, "ref_id", [f"TR{i:03d}" for i in range(len(R))]); return R, info


if __name__ == "__main__":
    for f in ("TLIB4096_v1.csv", "TREF_v1.csv"):
        if (OUT / f).exists():
            raise SystemExit(f"{f} exists (frozen); refusing to overwrite")
    L = library(); D = L[["z1", "z2", "z3"]].to_numpy()
    G = ang_deg(D, D); np.fill_diagonal(G, 180.0); nn = G.min(1)
    anti = np.linalg.norm(D[:N_AX] + D[N_AX:], axis=1).max()
    Gax = np.minimum(G[:N_AX, :N_AX], 180 - G[:N_AX, :N_AX]); np.fill_diagonal(Gax, 180); nnax = Gax.min(1)
    R, info = references(D)
    L.to_csv(OUT / "TLIB4096_v1.csv", index=False, float_format="%.17g"); R.to_csv(OUT / "TREF_v1.csv", index=False, float_format="%.17g")
    meta = dict(rule=__doc__, n_directions=len(D), n_axes=N_AX, unit_norm_max_dev=float(np.abs(np.linalg.norm(D, axis=1) - 1).max()),
                min_pairwise_deg=float(nn.min()), nn_deg_mean=float(nn.mean()), nn_deg_max=float(nn.max()), axial_nn_deg_min=float(nnax.min()),
                axial_nn_deg_mean=float(nnax.mean()), antipode_err=float(anti), n_exact_duplicates=int((nn < 1e-9).sum()),
                covering_radius_deg_est=float(ang_deg(fib_sphere(20000), D).min(1).max()), references=dict(n=len(R), per_set=R.set.value_counts().to_dict(), **info,
                min_deg_to_library=float(R.min_deg_to_library.min())),
                sha256={f: hashlib.sha256((OUT / f).read_bytes()).hexdigest() for f in ("TLIB4096_v1.csv", "TREF_v1.csv")})
    (OUT / "TLIB4096_v1_TREF_v1.json").write_text(json.dumps(meta, indent=1)); print(json.dumps({k: v for k, v in meta.items() if k != "rule"}, indent=1))
