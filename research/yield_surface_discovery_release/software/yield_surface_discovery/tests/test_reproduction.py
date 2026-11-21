"""Reproduction of the REGISTERED results (same-surrogate Arms C = physics-informed and A = generic blind) and of the minimal example."""
import json, sys
import numpy as np, pytest
from _helpers import ROOT, DATA, library, registered, eval_set, pinned_environment
from yield_surface_discovery import AdaptiveController, ControllerSettings, SHEAR_MIRROR, NO_SYMMETRY, CachedOracle, DirectGP
from yield_surface_discovery.oracle import analytical_oracle
from yield_surface_discovery.analytical import cb2004_radius
from yield_surface_discovery.acquisition import predictive_variance_scores, select
from yield_surface_discovery.metrics import evaluate

ARMS = {"C_entropy_sym": SHEAR_MIRROR, "A_entropy_nosym": NO_SYMMETRY}


def replay_registered(arm, start, max_steps=20):
    """feed the registered sequence and radii; at every acquisition step the registered choice must lie in the package's tie set."""
    lib = library(); ref = registered()[f"{arm}_s{start}"]; seq, r = ref["trajectory_original_ids"], np.array(ref["radii"])
    assert np.allclose(cb2004_radius(lib.directions[seq]), r, rtol=1e-14)
    for k in range(8, min(28, 8 + max_steps)):
        gp = DirectGP(ARMS[arm]).fit(lib.directions[seq[:k]], r[:k]); mask = np.zeros(len(lib), bool); mask[seq[:k]] = True
        sc, _, _ = predictive_variance_scores(gp, lib.directions, mask); _, tied = select(sc, rel_tol=1e-8)
        assert seq[k] in set(tied.tolist()), f"{arm} s{start}: registered test {k + 1} not among the package's maximisers"
    return min(20, max_steps)


@pytest.mark.parametrize("arm", list(ARMS))
@pytest.mark.parametrize("start", range(5))
def test_replay_registered_choices(arm, start):
    replay_registered(arm, start)


@pytest.mark.parametrize("arm", list(ARMS))
def test_registered_metrics_from_registered_points(arm):
    lib = library(); D, rref, M = eval_set()
    for s in range(5):
        ref = registered()[f"{arm}_s{s}"]; seq, r = ref["trajectory_original_ids"], np.array(ref["radii"])
        for n, v in ref["continuous_abs_p95"].items():
            n = int(n); m = evaluate(DirectGP(ARMS[arm]).fit(lib.directions[seq[:n]], r[:n]).predict_radius(D), rref, D, M)
            assert abs(m["global_p95"] - v) < 1e-9, (arm, s, n, m["global_p95"], v)


@pytest.mark.parametrize("arm", list(ARMS))
def test_bit_exact_legacy_rerun(arm):
    """re-running the controller in legacy mode reproduces the registered 28-test sequences bit-for-bit -- ONLY in the pinned environment
    (legacy mode lets floating-point round-off break exact ties, as the registered runs did)."""
    if not pinned_environment():
        pytest.xfail("numpy/scipy/scikit-learn differ from requirements-lock.txt: round-off ties may resolve differently in legacy mode")
    lib = library(); st = ControllerSettings(start_rule="registered_legacy", maximin_ties="legacy", acquisition_ties="legacy")
    for s in range(5):
        R = AdaptiveController(lib, CachedOracle(analytical_oracle("cb2004_mg4li")), ARMS[arm], st).run(s)
        assert R.sequence == registered()[f"{arm}_s{s}"]["trajectory_original_ids"]


def test_minimal_example_reproduces_expected(tmp_path):
    sys.path.insert(0, str(ROOT / "scripts")); import run_minimal_example as ME
    ok, got = ME.run(tmp_path / "minimal", check=True); assert ok, got
    assert (tmp_path / "minimal" / "run_start0.png").exists() and (tmp_path / "minimal" / "tests.csv").exists()
