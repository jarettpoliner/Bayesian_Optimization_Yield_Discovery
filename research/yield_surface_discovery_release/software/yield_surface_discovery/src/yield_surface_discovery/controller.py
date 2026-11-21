"""Adaptive controller (frozen algorithm).  It sees only candidate directions and YieldOracle observations.

  1. 8 initial tests by angular maximin over the candidate library (geometry only; no oracle values used);
  2. query the oracle once per selected test;
  3. fit the direct Matern-3/2 GP on log yield radius (+ shear-mirror pseudo-observations if the physics-informed symmetry is enabled);
  4. score every untested candidate by the GP predictive variance; select the maximum (deterministic tie handling);
  5. repeat 2-4 until the budget of REAL tests is reached.  Pseudo-observations never count as tests."""
from dataclasses import dataclass, asdict, field
import numpy as np, pandas as pd
from .gp import DirectGP
from .initialization import angular_maximin, first_candidate, LEGACY_EXCLUDE
from .acquisition import predictive_variance_scores, select
from .geometry import angle_rad
from .symmetry import NO_SYMMETRY


@dataclass
class ControllerSettings:
    n_initial: int = 8
    budget: int = 28
    start_rule: str = "seeded"                 # 'seeded' | 'registered_legacy'
    seed_base: int = 7000
    first_index: int = 0
    legacy_exclude: tuple = LEGACY_EXCLUDE
    maximin_ties: str = "deterministic"        # 'deterministic' | 'legacy'
    acquisition: str = "predictive_variance"
    acquisition_ties: str = "deterministic"    # 'deterministic' | 'legacy'
    tie_rel: float = 1e-10
    tie_abs: float = 1e-12
    maximin_tol_rad: float = 1e-9
    fit_direction: str = "requested"           # 'requested' (frozen) | 'achieved' (use the measured stress direction)
    achieved_warn_rad: float = 1e-6
    record_scores: bool = False


@dataclass
class RunResult:
    start: int
    sequence: list
    labels: list
    families: list
    observations: list
    steps: list
    symmetry: str
    settings: dict
    scores: list = field(default_factory=list)

    def fit(self, controller, n_tests):
        """surface GP refitted on the first n_tests real tests of this run (same settings as the controller)."""
        return controller.fit(self.observations[:n_tests], self.sequence[:n_tests])

    def directions(self, which="requested"):
        return np.array([o.requested_direction if which == "requested" else o.achieved_direction for o in self.observations])

    def radii(self):
        return np.array([o.yield_radius for o in self.observations])

    def tests_frame(self):
        rows = []
        for k, (i, lab, fam, o) in enumerate(zip(self.sequence, self.labels, self.families, self.observations)):
            rows.append(dict(start=self.start, test_number=k + 1, phase="angular-maximin initial" if k < self.settings["n_initial"] else "predictive uncertainty",
                             candidate_index=i, label=lab, family=fam, req_z1=o.requested_direction[0], req_z2=o.requested_direction[1], req_z3=o.requested_direction[2],
                             ach_z1=o.achieved_direction[0], ach_z2=o.achieved_direction[1], ach_z3=o.achieved_direction[2], s_xx=o.stress_at_yield[0],
                             s_yy=o.stress_at_yield[1], s_xy=o.stress_at_yield[2], yield_radius=o.yield_radius, status=o.status))
        return pd.DataFrame(rows)


class AdaptiveController:
    def __init__(self, library, oracle, symmetry=NO_SYMMETRY, settings=None, gp_kwargs=None):
        self.library, self.oracle, self.symmetry = library, oracle, symmetry
        self.settings = settings or ControllerSettings(); self.gp_kwargs = gp_kwargs or {}

    def new_gp(self):
        return DirectGP(symmetry=self.symmetry, **self.gp_kwargs)

    def initial_design(self, start, unavailable=()):
        s = self.settings; D = self.library.directions
        f = first_candidate(start, len(D), s.start_rule, s.seed_base, s.first_index, s.legacy_exclude)
        return angular_maximin(D, f, s.n_initial, s.maximin_ties, s.maximin_tol_rad, unavailable)

    def _observe(self, idx):
        return self._accept(self.oracle.evaluate(self.library.directions[idx]))

    def _accept(self, o):
        if o.success and np.isfinite(o.yield_radius) and o.yield_radius > 0:
            dev = float(angle_rad(o.requested_direction, o.achieved_direction)[0, 0])
            if dev > self.settings.achieved_warn_rad:
                o.metadata["achieved_vs_requested_rad"] = dev
            return o
        return None

    def fit(self, observations, indices=None):
        """fit the surface GP.  fit_direction='requested' uses the exact library directions of the tested candidates (frozen behaviour);
        'achieved' uses the measured stress directions reported by the oracle."""
        if self.settings.fit_direction == "requested" and indices is not None:
            D = self.library.directions[np.asarray(indices)]
        else:
            D = np.array([o.requested_direction if self.settings.fit_direction == "requested" else o.achieved_direction for o in observations])
        return self.new_gp().fit(D, np.array([o.yield_radius for o in observations]))

    def run(self, start=0):
        s = self.settings; D = self.library.directions; n = len(D)
        if s.acquisition != "predictive_variance":
            raise ValueError("the frozen method uses acquisition 'predictive_variance' only")
        unavailable, seq, obs, steps, scores = set(), [], [], [], []
        design = self.initial_design(start)                     # geometry only: no oracle call before this line returns
        pending = list(design)
        while pending:                                            # the initial tests are independent: run them as a batch
            for idx, o in zip(pending, self.oracle.evaluate_many(D[pending])):
                o = self._accept(o)
                if o is None:
                    unavailable.add(idx); steps.append(dict(test_number=len(seq) + 1, event="failed initial test (candidate marked unavailable)", candidate_index=idx))
                else:
                    seq.append(idx); obs.append(o)
            if len(seq) >= s.n_initial:
                break
            seeds = seq if seq else [min(i for i in range(n) if i not in unavailable)]      # extend the design from the successful tests
            ext = angular_maximin(D, seeds, s.n_initial if seq else s.n_initial, s.maximin_ties, s.maximin_tol_rad, unavailable)
            pending = [i for i in ext if i not in seq]
        while len(seq) < s.budget:
            gp = self.fit(obs, seq); mask = np.zeros(n, bool); mask[seq] = True; mask[list(unavailable)] = True
            sc, mean, var = predictive_variance_scores(gp, D, mask); j, tied = select(sc, s.acquisition_ties, s.tie_rel, s.tie_abs)
            hp = gp.hyperparameters()
            step = dict(test_number=len(seq) + 1, candidate_index=j, label=self.library.labels[j], score=float(sc[j]), gp_mean_log_radius=float(mean[j]),
                        gp_predictive_variance=float(var[j]), tie_set_size=int(len(tied)), tie_break_invoked=bool(len(tied) > 1),
                        n_training_points=hp["n_training_points"], hyperparameters_optimized=hp["optimized"], amplitude=hp.get("amplitude"),
                        length_scale=hp.get("length_scale"), white_noise=hp.get("white_noise"))
            if s.record_scores:
                scores.append(sc.copy())
            o = self._observe(j)
            if o is None:
                unavailable.add(j); step["event"] = "failed test (candidate marked unavailable)"; steps.append(step); continue
            seq.append(j); obs.append(o); steps.append(step)
        return RunResult(start, seq, [self.library.labels[i] for i in seq], [self.library.families[i] for i in seq], obs, steps, self.symmetry.name,
                         asdict(s), scores)
