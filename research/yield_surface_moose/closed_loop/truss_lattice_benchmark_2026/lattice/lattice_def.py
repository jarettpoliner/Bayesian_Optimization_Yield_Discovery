"""FROZEN truss-lattice definition TL1 (written 2026-11-20 before ANY MOOSE run of this workstream; never edited afterwards).

Architecture ("TL1: reinforced triangulated frame + kagome infill")
  * Fine triangular lattice of bar length a = 1 mm, decomposed (exactly) into
      - a coarse triangulated FRAME (spacing 2a; every other lattice line), and
      - a KAGOME INFILL (the remaining lattice lines).
  * Primitive periodic cell (rhombus): L1 = (2a, 0), L2 = (a, sqrt3 a); 4 nodes, 12 bars:
      N0 (0,0)            frame vertex (coordination 6, all frame bars)
      N1 (a,0), N2 (a/2, sqrt3 a/2), N3 (3a/2, sqrt3 a/2)   kagome nodes = mid-points of frame edges (2 frame + 4 kagome bars)
  * Bar families (6 families, 3 orientations):  FH = frame 0 deg, FI60 / FI120 = frame +-60 deg, K0 / K60 / K120 = kagome.
  * Two elastoplastic material groups: FRAME and INFILL.  Areas differ by family (x-reinforced frame) -> direction-dependent stiffness.
  * Mechanical stability: fully triangulated (every face a triangle), periodic static indeterminacy 12 - (2*4 - 2) - 3 = 3.
    Verified (not assumed) by the rank of the periodic equilibrium matrix in validation/independent_truss.py.

Rationale (physical, a priori; NOT tuned for any algorithm): a stiff, strong, x-reinforced frame carrying load through a statically indeterminate
network together with a weaker infill gives (i) anisotropic stiffness, (ii) load paths that switch between frame and infill with loading
direction, (iii) redistribution after first yield.  Parameters are round engineering values (steel-like frame, softer infill) fixed here once.

Mechanics: small-strain axial bars, J2 (= uniaxial) plasticity with LINEAR isotropic hardening (MOOSE PlasticTruss); no buckling, no damage,
no bending (a pin-jointed truss carries axial force only).
Units: mm, N, MPa.
"""
import json, hashlib
import numpy as np

A_BAR = 1.0                                   # bar length [mm]
S3 = np.sqrt(3.0)
L1 = np.array([2 * A_BAR, 0.0])
L2 = np.array([A_BAR, S3 * A_BAR])
NODES = np.array([[0.0, 0.0], [A_BAR, 0.0], [A_BAR / 2, S3 * A_BAR / 2], [3 * A_BAR / 2, S3 * A_BAR / 2]])
NODE_ROLE = ["frame_vertex", "kagome", "kagome", "kagome"]

MATERIALS = {   # E [MPa], yield stress [MPa], linear hardening modulus H (d sigma_y / d eps_p) [MPa]
    "FRAME":  dict(E=200000.0, sigma_y=300.0, H=2000.0),
    "INFILL": dict(E=200000.0, sigma_y=200.0, H=2000.0),
}
FAMILY = {      # family -> (material group, area mm^2, orientation deg)
    "FH":    ("FRAME", 0.05, 0),
    "FI60":  ("FRAME", 0.03, 60),
    "FI120": ("FRAME", 0.03, 120),
    "K0":    ("INFILL", 0.02, 0),
    "K60":   ("INFILL", 0.02, 60),
    "K120":  ("INFILL", 0.02, 120),
}
BLOCK_ID = {"FH": 1, "FI60": 2, "FI120": 3, "K0": 4, "K60": 5, "K120": 6}


def cell_area():
    return float(abs(L1[0] * L2[1] - L1[1] * L2[0]))


def bars():
    """every bond of the fine triangular lattice owned by the cell: from node i in directions 0, 60, 120 deg.
    returns list of dict(i, j, shift=(n1, n2)) meaning bar from X_i to X_j + n1 L1 + n2 L2, plus family."""
    out = []
    for i, Xi in enumerate(NODES):
        for th in (0, 60, 120):
            e = np.array([np.cos(np.radians(th)), np.sin(np.radians(th))]); P = Xi + A_BAR * e; hit = None
            for j, Xj in enumerate(NODES):
                for n1 in (-1, 0, 1):
                    for n2 in (-1, 0, 1):
                        if np.linalg.norm(Xj + n1 * L1 + n2 * L2 - P) < 1e-9:
                            hit = (j, n1, n2)
            assert hit is not None
            j, n1, n2 = hit
            frame = NODE_ROLE[i] == "frame_vertex" or NODE_ROLE[j] == "frame_vertex"
            fam = ("FH" if th == 0 else f"FI{th}") if frame else f"K{th}"
            out.append(dict(i=i, j=j, shift=(n1, n2), family=fam, theta_deg=th))
    return out


def summary():
    B = bars(); A = cell_area()
    rho = sum(FAMILY[b["family"]][1] * A_BAR for b in B) / A
    return dict(name="TL1", a=A_BAR, L1=L1.tolist(), L2=L2.tolist(), cell_area=A, nodes=NODES.tolist(), node_role=NODE_ROLE,
                materials=MATERIALS, families={k: dict(material=v[0], area=v[1], theta_deg=v[2]) for k, v in FAMILY.items()},
                bars=[dict(b, shift=list(b["shift"])) for b in B], n_bars=len(B), relative_density=rho,
                static_indeterminacy=len(B) - (2 * len(NODES) - 2) - 3)


if __name__ == "__main__":
    s = summary(); txt = json.dumps(s, indent=1)
    open(__file__.replace("lattice_def.py", "TL1_frozen.json"), "w").write(txt)
    print(txt); print("sha256(TL1_frozen.json) =", hashlib.sha256(txt.encode()).hexdigest())
