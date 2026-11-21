"""Frozen campaign design rules for the final B1/B2 benchmark (controller-side; blind: imports only the frozen BO/V2 code).
  BO starts / axial-maximin starts : K = (0, 1, 2, 3, 5, 6, 7, 8, 9)  -- the 9 starts with DISTINCT 8-axis initial designs (start 4 duplicates 2)
  first axis of start k            : axis of V1 candidate 0 (k = 0) or of default_rng(7000 + k).choice(414)   (V1/V2 rule)
  axial maximin                    : bomes_v2.axial_maximin over the 207 axis representatives, 33 + 12 spares
  random unique axes, seeds 1..10  : default_rng(9000 + k).choice(207, 45, replace=False); first 33 successful (12 spares)
Identical to homog_benchmark/nonadaptive.py (verified in the Stage-1 audit)."""
import sys
from pathlib import Path
import numpy as np
CL = Path(__file__).resolve().parents[2]
for p in (CL / "v2", CL):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))
from bomes_v2 import AxisLibrary, axial_maximin                                     # noqa: E402
from yield_surface_discovery.io import load_config, build_library                  # noqa: E402

BUDGET, SPARE = 33, 12
STARTS = (0, 1, 2, 3, 5, 6, 7, 8, 9)
RANDOM_SEEDS = tuple(range(1, 11))
VARIANTS = ("B1", "B2")


def library():
    lib = build_library(load_config(CL / "configs/S2_campaign.yaml")); return lib, AxisLibrary.from_library(lib)


def first_axis(AX, n_lib, start):
    return int(AX.class_of[0 if start == 0 else int(np.random.default_rng(7000 + start).choice(n_lib))])


def nonadaptive_designs(lib, AX):
    n = BUDGET + SPARE; out = {}
    for s in STARTS:
        out[f"axmm_s{s}"] = [int(c) for c in axial_maximin(AX.directions, first_axis(AX, len(lib.directions), s), n)]
    for k in RANDOM_SEEDS:
        out[f"random_{k}"] = [int(c) for c in np.random.default_rng(9000 + k).choice(len(AX.directions), n, replace=False)]
    return out
