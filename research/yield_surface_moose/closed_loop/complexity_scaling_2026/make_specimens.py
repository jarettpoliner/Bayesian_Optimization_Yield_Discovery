#!/usr/bin/env python3
"""Physical pilot geometries P2-P6 (PROTOCOL.md), meshed with the existing psbench.microstructure.build_mesh (gmsh 4.4, periodic node
matching verified by build_mesh).  Writes ONLY complexity_scaling_2026/specimens/<P>/mesh_h<nel>/.  psbench is imported read-only."""
import json, math, sys
from pathlib import Path
import numpy as np
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent; PS = HERE.parents[1]; sys.path.insert(0, str(PS))
from psbench.microstructure import Inclusion, build_mesh, place_inclusions, _fits          # noqa: E402  read-only
GAP, MARGIN = 0.03, 0.03


def regular(n=3, vf=0.15):
    r = math.sqrt(vf / (n * n * math.pi)); c = (np.arange(n) + 0.5) / n
    return [Inclusion(float(x), float(y), r, r, 0.0) for x in c for y in c]


def multiscale(seed=5501, vf=0.18):
    rng = np.random.default_rng(seed); placed, area, tries = [], 0.0, 0
    while tries < 200000 and area < vf - 0.005:
        tries += 1; a = float(rng.uniform(0.03, 0.11)); ar = float(rng.uniform(1.0, 4.0)); b = max(a / ar, 0.03)
        inc = Inclusion(float(rng.uniform(0, 1)), float(rng.uniform(0, 1)), a, b, float(rng.uniform(0, math.pi)))
        if area + inc.area > vf + 1e-12:
            continue
        if _fits(inc, placed, 1.0, GAP, MARGIN):
            placed.append(inc); area += inc.area
    return placed


SPECS = {
    "P2_single": dict(kind="single", incl=lambda: [Inclusion(0.5, 0.5, math.sqrt(0.15 / math.pi), math.sqrt(0.15 / math.pi), 0.0)]),
    "P3_regular": dict(kind="regular 3x3", incl=lambda: regular()),
    "P4_clustered_s4401": dict(kind="clustered", incl=lambda: place_inclusions("clustered", dict(vf_target=0.15, n_clusters=3, cluster_std=0.10, r_median=0.05, r_sigma=0.35,
                                                                                                   r_min=0.03, r_max=0.10, gap=GAP, margin=MARGIN), 4401)),
    "P5_multiscale_s5501": dict(kind="multiscale", incl=lambda: multiscale(5501)),
    "P6_aligned_hiAR_s6601": dict(kind="aligned AR6", incl=lambda: place_inclusions("aligned", dict(vf_target=0.15, a=0.18, b=0.03, theta_deg=0.0, gap=GAP, margin=MARGIN), 6601)),
}

if __name__ == "__main__":
    which = sys.argv[1:] or list(SPECS)
    for name in which:
        nels = (40, 80) if name.startswith("P6") else (40,)
        inc = SPECS[name]["incl"]()
        for nel in nels:
            d = HERE / "specimens" / name / f"mesh_h{nel}"
            if (d / "rve.msh").exists():
                continue
            build_mesh(d, inc, L=1.0, h=1.0 / nel)
            (d / "inclusions.json").write_text(json.dumps([i.__dict__ for i in inc], indent=1))
        st = json.loads((HERE / "specimens" / name / "mesh_h40/mesh_stats.json").read_text())
        print(name, len(inc), "inclusions", st)
