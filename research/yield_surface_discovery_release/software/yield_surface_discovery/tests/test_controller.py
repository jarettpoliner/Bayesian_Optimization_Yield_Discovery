import numpy as np
from _helpers import library
from yield_surface_discovery import AdaptiveController, ControllerSettings, SHEAR_MIRROR, NO_SYMMETRY, YieldOracle, CachedOracle
from yield_surface_discovery.oracle import AnalyticalOracle
from yield_surface_discovery.analytical import von_mises_radius


class Recorder(YieldOracle):
    def __init__(self, fail_first=0):
        self.inner = AnalyticalOracle(lambda D: von_mises_radius(D, 250.0)); self.log = []; self.fail_first = fail_first

    def evaluate(self, direction):
        self.log.append(np.asarray(direction).copy()); o = self.inner.evaluate(direction)
        if self.fail_first > 0:
            self.fail_first -= 1; o.success, o.status = False, "solver_failed"
        return o


def test_budget_counts_real_tests_only_and_no_oracle_before_design():
    lib = library(dedup=True); rec = Recorder(); c = AdaptiveController(lib, rec, SHEAR_MIRROR, ControllerSettings(budget=14))
    design = c.initial_design(0); assert len(rec.log) == 0                       # the initial design is geometry only
    R = c.run(0); assert R.sequence[:8] == design and len(R.sequence) == 14 and len(rec.log) == 14 and len(set(R.sequence)) == 14
    assert all(s["test_number"] >= 9 for s in R.steps)


def test_failed_tests_are_skipped():
    lib = library(dedup=True); rec = Recorder(fail_first=2); R = AdaptiveController(lib, rec, NO_SYMMETRY, ControllerSettings(budget=10)).run(0)
    assert len(R.sequence) == 10 and len(rec.log) == 12 and sum("failed" in s.get("event", "") for s in R.steps) == 2


def test_deterministic_mode_is_repeatable_and_cached_oracle_counts():
    lib = library(dedup=True); runs = []
    for _ in range(2):
        o = CachedOracle(AnalyticalOracle(lambda D: von_mises_radius(D, 1.0))); runs.append(AdaptiveController(lib, o, SHEAR_MIRROR, ControllerSettings(budget=12)).run(2))
        assert o.calls == 12
    assert runs[0].sequence == runs[1].sequence
