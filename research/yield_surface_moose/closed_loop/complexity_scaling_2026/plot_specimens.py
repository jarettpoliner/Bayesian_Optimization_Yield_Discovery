"""figure of the predefined physical pilot geometries P2-P6 (actual meshes) -> analysis/pilot_microstructures.{png,pdf,svg}"""
import sys, json
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection, LineCollection
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE.parents[1] / "presentation_draft/visualization/scripts"))
from common import read_msh
NAMES = {"P2_single": "P2 single circle", "P3_regular": "P3 regular 3×3", "P4_clustered_s4401": "P4 clustered (seed 4401)",
         "P5_multiscale_s5501": "P5 multiscale (seed 5501)", "P6_aligned_hiAR_s6601": "P6 aligned AR 6 (seed 6601)"}
fig, axs = plt.subplots(1, 5, figsize=(19, 4.3))
for ax, (k, t) in zip(axs, NAMES.items()):
    N, Q, ph = read_msh(HERE / "specimens" / k / "mesh_h40/rve.msh"); st = json.load(open(HERE / "specimens" / k / "mesh_h40/mesh_stats.json"))
    ax.add_collection(PolyCollection(N[Q], facecolors=np.where(ph[:, None] == 2, [[0.61, 0.29, 0.10, 1]], [[0.84, 0.86, 0.89, 1]]), edgecolors="#9a9a9a", linewidths=0.12))
    E = {}
    for e, q in enumerate(Q):
        for a, b in zip(q, np.roll(q, -1)):
            E.setdefault(tuple(sorted((a, b))), []).append(ph[e])
    ax.add_collection(LineCollection([N[list(kk)] for kk, v in E.items() if len(v) == 2 and v[0] != v[1]], colors="k", linewidths=0.8))
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([]); ax.set_title(t, fontsize=12)
    ax.set_xlabel(f"{st['n_quad']:,} QUAD4 · f = {st['vf_mesh']:.3f}", fontsize=10.5)
fig.suptitle("New physical pilot geometries (B1 phase properties, frozen oracle r1, h = 1/40; predefined in PROTOCOL.md)", fontsize=13.5)
(HERE / "analysis").mkdir(exist_ok=True)
for e in ("png", "pdf", "svg"):
    fig.savefig(HERE / f"analysis/pilot_microstructures.{e}", bbox_inches="tight", dpi=250)
