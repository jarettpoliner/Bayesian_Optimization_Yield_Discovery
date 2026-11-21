"""BO+MES Version 2: symmetry-aware maximum-entropy sampling (separate from, and importing read-only from, the frozen V1 package).

Unchanged from V1
  * surrogate: Matern-3/2 GP (sklearn ConstantKernel * Matern(nu=1.5) + WhiteKernel) on y = log r, constant mean = median(y),
    same hyper-parameter bounds, same fixed values below the ML-II threshold, same jitter, same dedup of identical stress points;
  * acquisition: maximum entropy sampling, alpha(d) = H[f(d) | D] = 0.5 log(2 pi e sigma_n^2(d)), i.e. d_next = argmax sigma_n^2(d);
  * tie handling: V1 acquisition.select (deterministic, rel 1e-10, abs 1e-12, lowest index wins).

Changed (each is a separate switch so its effect can be isolated)
  1. symmetry="antipodal": candidates are unoriented loading AXES {d, -d}; one representative per class (lowest V1 index); a tested axis is
     excluded as a whole; the GP receives the antipodal image of every real observation as a pseudo-observation (V1's own Symmetry mechanism,
     operation -I, pseudo-noise 1e-4) so the posterior is centrally symmetric.  The initial design is angular maximin on the AXIAL distance
     min(theta, pi - theta).  Only valid where r(d) = r(-d) has been verified.
  2. schedule: ML-II threshold and refit cadence counted in REAL unique observations (V1 counts training points incl. pseudo-observations).
  3. acquisition variance: LATENT posterior variance of f (WhiteKernel excluded).  V1 scores latent + white noise; the white term is the same
     constant for every candidate, so the argmax is identical (verified in audit_v1.py) -- this only makes the entropy definition explicit.
"""
from dataclasses import dataclass, field, asdict
import sys, warnings
from pathlib import Path
import numpy as np
from scipy.linalg import solve_triangular

BO_PKG = Path.home() / "projects/yield_surface_discovery_release/software/yield_surface_discovery"
if str(BO_PKG / "src") not in sys.path:
    sys.path.insert(0, str(BO_PKG / "src"))
from yield_surface_discovery.symmetry import Symmetry, NO_SYMMETRY            # noqa: E402  (read-only import of frozen V1 code)
from yield_surface_discovery.acquisition import select                         # noqa: E402
from yield_surface_discovery.geometry import angle_rad                         # noqa: E402

ANTIPODAL = Symmetry("antipodal", (-np.eye(3),))


# ------------------------------------------------------------------ candidate axes
@dataclass
class AxisLibrary:
    """antipodal equivalence classes of a V1 CandidateLibrary.  rep[c] = V1 index of the representative of class c; members[c] = V1 indices."""
    directions: np.ndarray            # (n_classes, 3) representative directions
    rep: np.ndarray
    members: list
    class_of: np.ndarray              # V1 index -> class id
    antipode_of: np.ndarray           # V1 index -> V1 index of -d (-1 if absent)
    antipode_err: np.ndarray          # |d_i + d_antipode| (numerical check)
    labels: list

    @classmethod
    def from_library(cls, lib, tol=1e-9):
        D = lib.directions; n = len(D)
        dist = np.linalg.norm(D[:, None, :] + D[None, :, :], axis=2)
        j = dist.argmin(1); err = dist[np.arange(n), j]; anti = np.where(err < tol, j, -1)
        if np.any(anti == np.arange(n)):
            raise ValueError("a direction is its own antipode (zero vector?)")
        if np.any((anti >= 0) & (anti[np.maximum(anti, 0)] != np.arange(n))):
            raise ValueError("antipode map is not an involution")
        class_of = -np.ones(n, int); rep, members = [], []
        for i in range(n):
            if class_of[i] >= 0:
                continue
            m = [i] + ([int(anti[i])] if anti[i] >= 0 else []); c = len(rep)
            class_of[m] = c; rep.append(i); members.append(m)
        rep = np.array(rep)
        return cls(D[rep], rep, members, class_of, anti, err, [lib.labels[i] for i in rep])

    def __len__(self):
        return len(self.rep)


def axial_angle(A, B):
    t = angle_rad(A, B); return np.minimum(t, np.pi - t)


def axial_maximin(Dax, seeds, n_total=8, tol=1e-9, unavailable=()):
    """V1 angular_maximin (deterministic ties, lowest index) with the axial distance min(theta, pi-theta)."""
    sel = [int(seeds)] if np.isscalar(seeds) else [int(i) for i in seeds]; bad = set(int(u) for u in unavailable)
    while len(sel) < n_total:
        dist = axial_angle(Dax, Dax[sel]).min(1); dist[sel] = -np.inf; dist[list(bad)] = -np.inf
        m = dist.max(); sel.append(int(np.where(dist >= m - tol)[0].min()))
    return sel


# ------------------------------------------------------------------ GP
@dataclass
class Schedule:
    """ML-II hyper-parameter schedule, counted in REAL unique observations.
    start: first n at which ML-II is run (None = never).  every: refit cadence after start (1 = every observation; k = at start, start+k, ...;
    between refits the last fitted hyper-parameters are held fixed)."""
    name: str = "A20"
    start: int = 20
    every: int = 1

    def refit_now(self, n):
        return self.start is not None and n >= self.start and (n - self.start) % self.every == 0

    def active(self, n):
        return self.start is not None and n >= self.start


SCHEDULES = {"A20": Schedule("A20", 20, 1), "B8": Schedule("B8", 8, 1), "C12": Schedule("C12", 12, 1), "D8p4": Schedule("D8p4", 8, 4),
             "FIXED": Schedule("FIXED", None, 1)}


class GPv2:
    """Matern-3/2 GP on log r with an optional symmetry (V1 Symmetry object) and a real-observation schedule."""

    def __init__(self, symmetry=NO_SYMMETRY, schedule=SCHEDULES["A20"], nu=1.5, length_scale=0.5, length_scale_bounds=(0.05, 5.0),
                 amplitude_bounds=(1e-5, 1e1), white_noise=1e-3, white_noise_bounds=(1e-6, 1e-1), n_restarts=2, random_state=0, jitter=1e-8,
                 dedup_decimals=9, held=None):
        self.symmetry, self.schedule, self.nu = symmetry, schedule, nu
        self.length_scale, self.ls_b, self.amp_b, self.white_noise, self.wn_b = length_scale, tuple(length_scale_bounds), tuple(amplitude_bounds), white_noise, tuple(white_noise_bounds)
        self.n_restarts, self.random_state, self.jitter, self.dedup_decimals = n_restarts, random_state, jitter, dedup_decimals
        self.held = held                                     # hyper-parameters (amp, ls, wn) to hold between periodic refits

    def fit(self, directions, radii):
        from sklearn.gaussian_process import GaussianProcessRegressor
        from sklearn.gaussian_process.kernels import Matern, WhiteKernel, ConstantKernel
        D = np.atleast_2d(np.asarray(directions, float)); r = np.asarray(radii, float)
        if np.any(~np.isfinite(r)) or np.any(r <= 0):
            raise ValueError("yield radii must be finite and positive")
        P = D * r[:, None]; _, ui = np.unique(np.round(P, self.dedup_decimals), axis=0, return_index=True); ui = np.sort(ui)
        self.d_obs_, self.r_obs_ = D[ui], r[ui]; y = np.log(self.r_obs_); self.mean_ = float(np.median(y)); yc = y - self.mean_
        X, Y, alpha = self.symmetry.augment(self.d_obs_, yc); self.X_train_, self.y_train_ = X, Y
        a = self.jitter + (0 if alpha is None else alpha); n = len(self.d_obs_); amp0 = max(float(np.var(Y)), 1e-4)
        if self.schedule.refit_now(n):
            k = ConstantKernel(amp0, self.amp_b) * Matern(self.length_scale, self.ls_b, nu=self.nu) + WhiteKernel(self.white_noise, self.wn_b)
            gp = GaussianProcessRegressor(k, alpha=a, n_restarts_optimizer=self.n_restarts, random_state=self.random_state); self.mode_ = "ml2"
        elif self.schedule.active(n) and self.held is not None:
            amp, ls, wn = self.held
            k = ConstantKernel(amp, "fixed") * Matern(ls, "fixed", nu=self.nu) + WhiteKernel(wn, "fixed")
            gp = GaussianProcessRegressor(k, alpha=a, optimizer=None); self.mode_ = "held"
        else:
            k = ConstantKernel(amp0, "fixed") * Matern(self.length_scale, "fixed", nu=self.nu) + WhiteKernel(self.white_noise, "fixed")
            gp = GaussianProcessRegressor(k, alpha=a, optimizer=None); self.mode_ = "fixed"
        with warnings.catch_warnings():
            warnings.simplefilter("ignore"); gp.fit(X, Y)
        self.gp_ = gp; return self

    # -- hyper-parameters
    def theta(self):
        p = self.gp_.kernel_.get_params()
        return float(p["k1__k1__constant_value"]), float(p["k1__k2__length_scale"]), float(p["k2__noise_level"])

    def hyperparameters(self):
        amp, ls, wn = self.theta()
        at_b = [nm for nm, v, (lo, hi) in (("amplitude", amp, self.amp_b), ("length_scale", ls, self.ls_b), ("white_noise", wn, self.wn_b))
                if self.mode_ == "ml2" and (v <= lo * 1.0001 or v >= hi * 0.9999)]
        return dict(mode=self.mode_, n_real=int(len(self.d_obs_)), n_training_points=int(len(self.X_train_)), amplitude=amp, length_scale=ls,
                    white_noise=wn, at_bounds=at_b, log_marginal_likelihood=float(self.gp_.log_marginal_likelihood_value_))

    # -- prediction
    def predict_log(self, directions):
        """mean of log r and LATENT variance of f (no white noise)."""
        X = np.atleast_2d(directions); g = self.gp_; k1 = g.kernel_.k1
        Ks = k1(X, g.X_train_); mu = Ks @ g.alpha_
        v = solve_triangular(g.L_, Ks.T, lower=True, check_finite=False)
        var = np.maximum(k1.diag(X) - np.einsum("ij,ij->j", v, v), 0.0)
        return self.mean_ + mu, var

    def predict_total_var(self, directions):
        return self.predict_log(directions)[1] + self.theta()[2]

    def predict_radius(self, directions):
        return np.exp(self.predict_log(directions)[0])


def entropy(var_latent):
    return 0.5 * np.log(2 * np.pi * np.e * np.maximum(var_latent, 1e-300))


# ------------------------------------------------------------------ controller
@dataclass
class V2Settings:
    symmetry: str = "antipodal"          # 'antipodal' | 'none'
    schedule: str = "A20"
    n_initial: int = 8
    budget: int = 33                     # REAL tests
    first_index: int = 0                 # V1 index of the first initial test (start 0); start s >= 1 -> default_rng(7000 + s) as in V1
    seed_base: int = 7000
    tie_rel: float = 1e-10
    tie_abs: float = 1e-12
    maximin_tol_rad: float = 1e-9
    n_restarts: int = 2
    record_scores: bool = True


class V2Controller:
    """sees only the candidate directions and oracle observations (same YieldOracle interface as V1)."""

    def __init__(self, lib, oracle, settings=None):
        self.lib, self.oracle, self.s = lib, oracle, settings or V2Settings()
        self.sym = ANTIPODAL if self.s.symmetry == "antipodal" else NO_SYMMETRY
        if self.s.symmetry == "antipodal":
            self.ax = AxisLibrary.from_library(lib)
            self.C = self.ax.directions; self.cidx = self.ax.rep          # candidate set = axis representatives
        else:
            self.ax = None; self.C = lib.directions; self.cidx = np.arange(len(lib.directions))
        self.held = None

    def first(self, start):
        n1 = len(self.lib.directions)
        f = self.s.first_index if start == 0 else int(np.random.default_rng(self.s.seed_base + start).choice(n1))
        return int(self.ax.class_of[f]) if self.ax is not None else f

    def initial_design(self, start, unavailable=()):
        f = self.first(start)
        if self.ax is not None:
            return axial_maximin(self.C, f, self.s.n_initial, self.s.maximin_tol_rad, unavailable)
        from yield_surface_discovery.initialization import angular_maximin
        return angular_maximin(self.C, f, self.s.n_initial, "deterministic", self.s.maximin_tol_rad, unavailable)

    def make_gp(self, schedule=None):
        return GPv2(self.sym, SCHEDULES[schedule or self.s.schedule], n_restarts=self.s.n_restarts, held=self.held)

    def fit(self, obs_dirs, radii, schedule=None):
        gp = self.make_gp(schedule).fit(obs_dirs, radii)
        if gp.mode_ == "ml2":
            self.held = gp.theta()
        return gp

    def _ok(self, o):
        return o is not None and o.success and np.isfinite(o.yield_radius) and o.yield_radius > 0

    def run(self, start=0):
        s = self.s; seq, obs, steps, scores, unavailable = [], [], [], [], set()
        pending = self.initial_design(start)
        while pending:
            for c, o in zip(pending, self.oracle.evaluate_many(self.C[pending])):
                if self._ok(o):
                    seq.append(c); obs.append(o)
                else:
                    unavailable.add(c); steps.append(dict(test_number=len(seq) + 1, event="failed initial test", candidate=c))
            if len(seq) >= s.n_initial:
                break
            ext = (axial_maximin(self.C, seq, s.n_initial, s.maximin_tol_rad, unavailable) if self.ax is not None else
                   __import__("yield_surface_discovery.initialization", fromlist=["x"]).angular_maximin(self.C, seq, s.n_initial, "deterministic", s.maximin_tol_rad, unavailable))
            pending = [i for i in ext if i not in seq]
        while len(seq) < s.budget:
            gp = self.fit(self.C[seq], [o.yield_radius for o in obs])
            mean, var = gp.predict_log(self.C); sc = var.copy(); sc[seq] = -np.inf; sc[list(unavailable)] = -np.inf
            j, tied = select(sc, "deterministic", s.tie_rel, s.tie_abs); hp = gp.hyperparameters()
            step = dict(test_number=len(seq) + 1, candidate=int(j), v1_index=int(self.cidx[j]), label=self.lib.labels[self.cidx[j]],
                        latent_var=float(var[j]), entropy=float(entropy(var[j])), gp_mean_log_radius=float(mean[j]), tie_set_size=int(len(tied)),
                        tie_set=[int(t) for t in tied], **{k: v for k, v in hp.items()})
            if s.record_scores:
                scores.append(sc.copy())
            o = self.oracle.evaluate(self.C[j])
            if not self._ok(o):
                unavailable.add(j); step["event"] = "failed test (axis marked unavailable)"; steps.append(step); continue
            seq.append(int(j)); obs.append(o); steps.append(step)
        return dict(start=start, sequence=seq, v1_indices=[int(self.cidx[c]) for c in seq], labels=[self.lib.labels[self.cidx[c]] for c in seq],
                    observations=obs, steps=steps, scores=scores, settings=asdict(s))
