"""Adapter between a BO controller and the MOOSE yield experiment.

The BO controller sees a black box f(d) -> observation on the unit sphere of Mandel
stress directions d = sigma_M / |sigma_M|, sigma_M = (sigma_xx, sigma_yy, sqrt2 tau_xy).
It never receives specimen parameters; it only names the specimen by its id.

NOTE: the existing BO implementation was not found on this machine, so its exact
calling convention (array shapes, direction ordering, whether it consumes r or
log r, failure handling) is not yet verified. This adapter exposes the conventions
stated in the project brief; `MooseYieldOracle.__call__` is the single place to
adapt once that code is available. The BO mathematics is not touched here.
"""

import numpy as np

from .observe import observe


class MooseYieldOracle:
    """Callable black box for one fixed specimen.

    oracle(d)          -> yield radius r(d) (float); NaN if the trial was rejected
    oracle.log_radius  -> log r(d), the GP target in the log-yield-radius representation
    oracle.history     -> list of (direction, observation) pairs, in call order
    """

    def __init__(self, specimen_id, reject="nan", **oracle_options):
        """reject: 'nan' returns NaN for rejected trials, 'raise' raises RuntimeError."""
        self.specimen_id = specimen_id
        self.reject = reject
        self.options = oracle_options
        self.history = []

    def observe(self, direction):
        d = np.asarray(direction, float).reshape(3)
        n = np.linalg.norm(d)
        if not np.isfinite(n) or n == 0:
            raise ValueError("direction must be a nonzero 3-vector")
        obs = observe(d / n, self.specimen_id, **self.options)
        self.history.append((d / n, obs))
        if not obs["accepted"] and self.reject == "raise":
            raise RuntimeError(f"trial {obs['trial_id']} rejected")
        return obs

    def __call__(self, direction):
        obs = self.observe(direction)
        return obs["yield_radius"] if obs["accepted"] else float("nan")

    def log_radius(self, direction):
        obs = self.observe(direction)
        return obs["log_yield_radius"] if obs["accepted"] else float("nan")

    def training_data(self):
        """(D, log_r) arrays of accepted observations, for a GP update."""
        acc = [(d, o["log_yield_radius"]) for d, o in self.history if o["accepted"]]
        if not acc:
            return np.zeros((0, 3)), np.zeros(0)
        D, y = zip(*acc)
        return np.array(D), np.array(y)
