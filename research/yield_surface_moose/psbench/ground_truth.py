"""HIDDEN ground-truth specimen definitions (simulator side only).

Nothing in this module may be passed to the BO controller. The BO side only ever
sees `specimen_id` strings and the observation dictionaries returned by
`psbench.observe`.

Units: MPa for moduli and stresses; RVE edge length L = 1 (results are scale-free).
"""

import json
from pathlib import Path

from .microstructure import build_mesh, place_inclusions

ROOT = Path(__file__).resolve().parents[1]
SPECIMEN_DIR = ROOT / "specimens"

# Common placement controls: minimum inclusion-inclusion gap and boundary margin,
# both ~1 element (h = 1/40) so every matrix ligament is resolved by >= 1 element.
_PLACE = dict(gap=0.03, margin=0.03)

SPECIMENS = {
    # Homogeneous J2 material used only for analytical validation.
    "homog": dict(
        family="none",
        seed=0,
        geometry={},
        matrix=dict(E=70000.0, nu=0.33, sy=250.0, H=1000.0),
        inclusion=None,
    ),
    # A. Random particulate: equal circles, uniform positions (Al / SiC-like).
    "S1": dict(
        family="particulate",
        seed=1101,
        geometry=dict(vf_target=0.15, radius=0.08, **_PLACE),
        matrix=dict(E=70000.0, nu=0.33, sy=250.0, H=1500.0),
        inclusion=dict(E=400000.0, nu=0.20),
    ),
    # B. Aligned fibers: aspect-3 ellipses at a common 30 deg orientation.
    "S2": dict(
        family="aligned",
        seed=2202,
        geometry=dict(vf_target=0.18, a=0.12, b=0.04, theta_deg=30.0, **_PLACE),
        matrix=dict(E=72000.0, nu=0.34, sy=200.0, H=800.0),
        inclusion=dict(E=230000.0, nu=0.22),
    ),
    # C. Clustered: lognormal radii, inclusions drawn around 3 cluster centers.
    "S3": dict(
        family="clustered",
        seed=3303,
        geometry=dict(vf_target=0.15, n_clusters=3, cluster_std=0.10, r_median=0.05,
                      r_sigma=0.35, r_min=0.03, r_max=0.10, **_PLACE),
        matrix=dict(E=65000.0, nu=0.32, sy=300.0, H=2500.0),
        inclusion=dict(E=300000.0, nu=0.25),
    ),
}


def mesh_dir(specimen_id, nel=40):
    return SPECIMEN_DIR / specimen_id / f"mesh_h{nel}"


def build_specimen(specimen_id, nel=40):
    """Generate (or reuse) the fixed mesh of a specimen at resolution h = 1/nel.

    The inclusion layout depends only on the specimen seed, never on nel, so meshes
    at different resolutions describe the same geometry (used for mesh convergence).
    """
    spec = SPECIMENS[specimen_id]
    d = mesh_dir(specimen_id, nel)
    msh = d / "rve.msh"
    if not msh.exists():
        incl = [] if spec["family"] == "none" else place_inclusions(
            spec["family"], spec["geometry"], spec["seed"])
        build_mesh(d, incl, L=1.0, h=1.0 / nel)
        (SPECIMEN_DIR / specimen_id / "truth.json").write_text(json.dumps(spec, indent=2))
    stats = json.loads((d / "mesh_stats.json").read_text())
    return msh, stats, spec
