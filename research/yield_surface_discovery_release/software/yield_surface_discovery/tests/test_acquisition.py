import numpy as np
from _helpers import library
from yield_surface_discovery import DirectGP
from yield_surface_discovery.acquisition import predictive_variance_scores, select


def test_scores_are_gp_predictive_variance():
    lib = library(); rng = np.random.default_rng(2); idx = rng.choice(len(lib), 9, replace=False)
    gp = DirectGP().fit(lib.directions[idx], 100 * np.exp(rng.normal(0, .2, 9))); mask = np.zeros(len(lib), bool); mask[idx] = True
    sc, mean, var = predictive_variance_scores(gp, lib.directions, mask)
    assert np.all(np.isneginf(sc[idx])) and np.allclose(sc[~mask], var[~mask]) and np.allclose(var, gp.predict_log(lib.directions)[1])
    j, tied = select(sc); assert j in tied and sc[j] >= sc[~mask].max() - 1e-12


def test_deterministic_tie_rule():
    s = np.array([0.5, 1.0, 1.0 - 1e-13, 0.9, 1.0, -np.inf])
    j, tied = select(s); assert j == 1 and list(tied) == [1, 2, 4]
    j, tied = select(np.array([0.2, 1.0 - 1e-6, 1.0])); assert j == 2 and list(tied) == [2]
    assert select(np.array([1.0, 1.0]), ties="legacy")[0] == 0
