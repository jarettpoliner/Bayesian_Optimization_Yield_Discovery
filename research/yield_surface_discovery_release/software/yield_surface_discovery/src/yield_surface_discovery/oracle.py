"""Solver-independent oracle interface.  The adaptive controller calls ONLY YieldOracle.evaluate(direction).

An oracle executes one proportional loading test along a requested unit Mandel direction and reports first yield.  Implementations:
  AnalyticalOracle (reproducible benchmark), solvers.external.ExternalCommandOracle (any command-line solver),
  solvers.moose.MooseYieldOracle (MOOSE adapter; optional, MOOSE not required), or your own subclass (e.g. laboratory hardware)."""
from abc import ABC, abstractmethod
from dataclasses import dataclass, field, asdict
import numpy as np
from .geometry import mandel_to_cauchy, unit


@dataclass
class YieldObservation:
    requested_direction: np.ndarray            # unit Mandel direction that was requested
    achieved_direction: np.ndarray             # unit Mandel direction of the measured stress at yield
    stress_at_yield: np.ndarray                # Cauchy (s_xx, s_yy, s_xy) at first yield
    yield_radius: float                        # Mandel norm of the yield stress
    success: bool = True
    status: str = "ok"
    metadata: dict = field(default_factory=dict)

    def to_dict(self):
        d = asdict(self)
        for k in ("requested_direction", "achieved_direction", "stress_at_yield"):
            d[k] = [float(v) for v in np.asarray(d[k]).ravel()]
        return d


class YieldOracle(ABC):
    @abstractmethod
    def evaluate(self, direction) -> YieldObservation:
        """run one loading test along `direction` (unit Mandel 3-vector) and return the first-yield observation."""

    def evaluate_many(self, directions):
        """run several independent tests (the initial design).  Default: one evaluate() per direction; oracles may run them as a batch."""
        return [self.evaluate(d) for d in np.atleast_2d(directions)]


class AnalyticalOracle(YieldOracle):
    """exact yield radius from a closed-form model r(d) (the formula stays inside the oracle)."""
    def __init__(self, radius_fn, name="analytical"):
        self.radius_fn, self.name = radius_fn, name

    def evaluate(self, direction):
        return self.evaluate_many(np.atleast_2d(direction))[0]

    def evaluate_many(self, directions):
        D = unit(np.atleast_2d(np.asarray(directions, float))); R = np.asarray(self.radius_fn(D), float); out = []
        for d, r in zip(D, R):
            ok = bool(np.isfinite(r) and r > 0)
            out.append(YieldObservation(d, d.copy(), mandel_to_cauchy(r * d), float(r), ok, "ok" if ok else "no_yield", dict(oracle=self.name)))
        return out


class CachedOracle(YieldOracle):
    """runs each physical test at most once (repeat requests return the stored observation); counts real experiments."""
    def __init__(self, oracle, decimals=12):
        self.oracle, self.decimals, self.cache, self.calls = oracle, decimals, {}, 0

    def _key(self, d):
        return tuple(np.round(np.asarray(d, float), self.decimals))

    def evaluate(self, direction):
        k = self._key(direction)
        if k not in self.cache:
            self.cache[k] = self.oracle.evaluate(direction); self.calls += 1
        return self.cache[k]

    def evaluate_many(self, directions):
        D = np.atleast_2d(directions); new = [d for d in D if self._key(d) not in self.cache]
        if new:
            for d, o in zip(new, self.oracle.evaluate_many(np.array(new))):
                self.cache[self._key(d)] = o; self.calls += 1
        return [self.cache[self._key(d)] for d in D]


def analytical_oracle(model, **params):
    from .analytical import MODELS
    if model not in MODELS:
        raise ValueError(f"unknown analytical model {model!r}; available: {sorted(MODELS)}")
    return AnalyticalOracle(lambda D: MODELS[model](D, **params), name=model)
