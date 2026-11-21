import json
import numpy as np
from _helpers import DATA, library
from yield_surface_discovery import DirectGP, SHEAR_MIRROR, NO_SYMMETRY


def test_matches_registered_reconstructor():
    """predictions equal those of the registered implementation (blind_lib.fit_blind) recorded in data/reproduction/gp_reference.json"""
    ref = json.loads((DATA / "reproduction" / "gp_reference.json").read_text()); P2 = library().directions
    Q = library("arbitrary_direction_reference").directions[:200]
    for c in ref["cases"]:
        gp = DirectGP(SHEAR_MIRROR if c["symmetry"] else NO_SYMMETRY).fit(P2[c["candidate_ids"]], np.array(c["radii"]))
        assert np.allclose(gp.predict_radius(Q), c["pred_radius"], rtol=1e-6, atol=0)
        assert np.allclose(np.sqrt(gp.predict_log(Q)[1]), c["pred_log_sd"], rtol=1e-5, atol=1e-9)
        assert abs(gp.mean_ - c["log_mean"]) < 1e-12


def test_frozen_hyperparameter_schedule():
    P2 = library().directions; rng = np.random.default_rng(0); idx = rng.choice(len(P2), 12, replace=False); r = 200 * np.exp(rng.normal(0, .2, 12))
    hp = DirectGP().fit(P2[idx], r).hyperparameters()
    assert not hp["optimized"] and hp["length_scale"] == 0.5 and hp["white_noise"] == 1e-3 and hp["n_training_points"] == 12
    hp = DirectGP(SHEAR_MIRROR).fit(P2[idx], r).hyperparameters()            # mirror pseudo-observations count towards the 20-point switch
    assert hp["n_training_points"] > 12 and hp["n_real_observations"] == 12
    idx = rng.choice(len(P2), 22, replace=False); assert DirectGP().fit(P2[idx], 200 * np.exp(rng.normal(0, .2, 22))).hyperparameters()["optimized"]


def test_symmetric_fit_is_mirror_symmetric_and_log_space():
    P2 = library().directions; rng = np.random.default_rng(1); idx = rng.choice(len(P2), 10, replace=False); r = 300 * np.exp(rng.normal(0, .3, 10))
    gp = DirectGP(SHEAR_MIRROR).fit(P2[idx], r); Q = library("arbitrary_direction_reference").directions[4:104]
    assert np.allclose(gp.predict_radius(Q), gp.predict_radius(Q * [1, 1, -1]), rtol=2e-3)
    assert np.all(gp.predict_radius(Q) > 0)
