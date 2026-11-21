#!/usr/bin/env python3
"""Freeze the independent reference design (same 72 axes for B1 and B2): refs/reference_design.{csv,json,sha256}.
  60 area-stratified axes: 5 isotropic principal-stress regions (homog_benchmark/analytic_reference.principal_region), at least 6 per region,
     remainder proportional to the region's share of S^2 (equal_biaxial 6, tension_compression 20, uniaxial 15, pure_shear 11, biaxial_same_sign 8);
     candidates: 40,000-point Fibonacci set, kept if the axial angle to EVERY one of the 207 library axes is >= 3 deg; within a region,
     deterministic greedy farthest-point (axial angle) selection starting from the lowest-index candidate.
  12-axis continuity sweep: great circle d(t) = cos t (1,1,0)/sqrt2 + sin t (0,0,1) (equal biaxial -> pure shear), t_k = phi + 7.5 k deg,
     k = 0..11, smallest phi in {1.0, 1.25, ..., 7.25} deg with every axis >= 2 deg from all library axes (3 deg is not attainable on this circle
     if the check below says so; the attained minimum is recorded)."""
import hashlib, json, sys
from pathlib import Path
import numpy as np, pandas as pd
FB = Path(__file__).resolve().parent; CL = FB.parents[1]; sys.path.insert(0, str(FB)); sys.path.insert(0, str(CL / "homog_benchmark"))
import designs_common as D                                                          # noqa: E402
import analytic_reference as AR                                                     # noqa: E402  (region definition only)
OUT = FB / "refs"; OUT.mkdir(exist_ok=True)
if (OUT / "reference_design.csv").exists():
    raise SystemExit("frozen")
lib, AX = D.library(); L = AX.directions / np.linalg.norm(AX.directions, axis=1, keepdims=True)
ax_deg = lambda P: np.degrees(np.arccos(np.clip(np.abs(P @ L.T), 0, 1))).min(1)
P = AR.fibonacci_sphere(40000); P = P[P[:, 2] >= 0]; reg = AR.principal_region(P)[0]; dist = ax_deg(P); ok = dist >= 3.0
quota = {"equal_biaxial": 6, "tension_compression": 20, "uniaxial": 15, "pure_shear": 11, "biaxial_same_sign": 8}
rows = []
for rg, q in quota.items():
    idx = np.where(ok & (reg == rg))[0]; C = P[idx]; sel = [0]
    while len(sel) < q:
        dmin = np.degrees(np.arccos(np.clip(np.abs(C @ C[sel].T), 0, 1))).min(1); dmin[sel] = -1; sel.append(int(np.argmax(dmin)))
    for j in sel:
        rows.append(dict(set="area_stratified", region=rg, z1=C[j, 0], z2=C[j, 1], z3=C[j, 2], min_deg_to_library=float(dist[idx[j]])))
e = np.array([1, 1, 0]) / np.sqrt(2); s = np.array([0, 0, 1.0]); best = None
for phi in np.arange(1.0, 7.5, 0.25):
    t = np.radians(phi + 7.5 * np.arange(12)); Q = np.cos(t)[:, None] * e + np.sin(t)[:, None] * s; m = ax_deg(Q).min()
    if m >= 2.0:
        best = (phi, Q); break
assert best is not None, "no admissible continuity-sweep phase"
for k, q in enumerate(best[1]):
    rows.append(dict(set="continuity_sweep", region=AR.principal_region(q[None])[0][0], z1=q[0], z2=q[1], z3=q[2], min_deg_to_library=float(ax_deg(q[None])[0]),
                     sweep_t_deg=best[0] + 7.5 * k))
df = pd.DataFrame(rows); df.insert(0, "ref_id", [f"FBR{i + 1:02d}" for i in range(len(df))]); A = df[["z1", "z2", "z3"]].to_numpy()
G = np.degrees(np.arccos(np.clip(np.abs(A @ A.T), 0, 1))); np.fill_diagonal(G, 180)
info = dict(n=len(df), per_set=df.set.value_counts().to_dict(), per_region_area=df[df.set == "area_stratified"].region.value_counts().to_dict(), sweep_phase_deg=best[0],
            min_deg_to_library_area=float(df[df.set == "area_stratified"].min_deg_to_library.min()), min_deg_to_library_sweep=float(df[df.set == "continuity_sweep"].min_deg_to_library.min()),
            min_pairwise_axial_deg=float(G.min()), region_share_of_sphere=pd.Series(AR.principal_region(AR.fibonacci_sphere(20000))[0]).value_counts(normalize=True).to_dict())
df.to_csv(OUT / "reference_design.csv", index=False, float_format="%.17g"); json.dump(dict(rule=__doc__, **info), open(OUT / "reference_design.json", "w"), indent=1)
(OUT / "reference_design.sha256").write_text(hashlib.sha256((OUT / "reference_design.csv").read_bytes()).hexdigest() + "  reference_design.csv\n"); print(json.dumps(info, indent=1))
