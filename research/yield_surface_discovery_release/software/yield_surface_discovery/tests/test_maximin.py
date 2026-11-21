import numpy as np
from _helpers import library, registered
from yield_surface_discovery.initialization import angular_maximin, first_candidate
from yield_surface_discovery.geometry import angle_rad


def test_legacy_maximin_reproduces_registered_initial_designs():
    lib = library(); reg = registered()
    for s in range(5):
        f = first_candidate(s, len(lib), "registered_legacy")
        assert angular_maximin(lib.directions, f, 8, ties="legacy") == reg[f"C_entropy_sym_s{s}"]["trajectory_original_ids"][:8]


def test_deterministic_maximin_properties():
    lib = library(dedup=True); D = lib.directions; sel = angular_maximin(D, 0, 8)
    assert len(set(sel)) == 8 and sel[0] == 0
    for k in range(1, 8):                                              # each pick maximises the min angle to the previous picks
        d = angle_rad(D, D[sel[:k]]).min(1); d[sel[:k]] = -1
        assert abs(d[sel[k]] - d.max()) < 1e-9 and sel[k] == int(np.where(d >= d.max() - 1e-9)[0].min())


def test_first_candidate_rules():
    assert first_candidate(0, 414) == 0 and first_candidate(3, 414) == first_candidate(3, 414)
    assert first_candidate(1, 416, "registered_legacy") not in (0, 18, 9, 80)
