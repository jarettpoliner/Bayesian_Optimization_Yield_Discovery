"""BO V3: leave-one-out-residual-weighted variance (LRV) acquisition.  Separate, versioned module (v3.0); the frozen V2 controller, GP,
kernel, schedule, symmetry, library handling and tie rule are imported READ-ONLY from closed_loop/v2/bomes_v2.py and the V1 package.

Exactly ONE thing differs from frozen V2 BO+MES: the score that is maximised.
  V2:  s(x) = sigma_n^2(x)                                   (latent posterior variance; maximum entropy sampling)
  V3:  s(x) = sigma_n^2(x) * rho_n(x)                         (LRV)

Specification (fixed before any test; see PROTOCOL.md section 2):
  GP fit       identical to V2Controller.fit (A20 schedule, antipodal pseudo-observations, held hyper-parameters, bounds, restarts); sigma_n^2 is
               the frozen GPv2.predict_log latent variance.  Only observations already made are used.
  LOO          for each REAL observation i, the closed-form block leave-out of B_i = {x_i, -x_i (its antipodal pseudo-observation, if present)} from
               the fitted GP at its current hyper-parameters:
                 Q = K^-1 (K = training covariance incl. white noise, jitter and pseudo-noise, from the fitted Cholesky factor),
                 e_B = (Q_BB)^-1 (Q y)_B       (LOO residuals y_B - mu_{-B}(x_B)),     C_B = (Q_BB)^-1   (LOO predictive covariance),
                 u_i = e_B[real] / sqrt(C_B[real, real])    (standardised LOO residual of observation i).
  Interpolation rho_n(x) = sum_i c_i(x) u_i^2 / sum_i c_i(x),  c_i(x) = max( m(x, x_i), m(x, -x_i) ),  m = Matern(nu = 1.5, length 0.5).
               The interpolation length is V2's own frozen pre-threshold length scale (GPv2 default length_scale = 0.5), NOT the current ML-II
               length scale: ML-II saturates at its upper bound 5.0 in most steps, where m >= 0.85 between ANY two axes, so the weights would be
               nearly uniform and V3 would collapse to V2 by construction.  (Fixed a priori from the kernel and the bound; no V3 result used.)
  Normalisation rho_n <- rho_n / mean_i(u_i^2)  (score scale only; argmax unaffected).
  Safeguards   (numerical only) if mean(u^2) is not finite or <= 1e-300, or any u_i is not finite, rho_n := 1 (V2 behaviour, flagged);
               LOO variances floored at 1e-300; rho_n floored at 1e-12; tested / unavailable axes get -inf as in V2.
  Ties         frozen V1 acquisition.select (deterministic, rel 1e-10, abs 1e-12).
  Reduction    if all u_i are equal (or uniform_risk=True), rho_n == 1 exactly and the sequence equals V2 BO+MES bitwise.
No free parameter is introduced and none may be changed after the protocol freeze."""
import sys
from dataclasses import asdict
from pathlib import Path
import numpy as np
from scipy.linalg import cho_solve
from sklearn.gaussian_process.kernels import Matern

CL = Path(__file__).resolve().parents[1]
for p in (str(CL / "v2"),):
    if p not in sys.path:
        sys.path.insert(0, p)
import bomes_v2 as B                                                                  # noqa: E402  frozen
from bomes_v2 import axial_maximin, entropy, GPv2                                    # noqa: E402  frozen
from yield_surface_discovery.acquisition import select                              # noqa: E402  frozen V1
VERSION = "v3.0-LRV"
RISK_LENGTH = GPv2().length_scale                                                    # 0.5, the frozen V2 pre-threshold length scale
assert RISK_LENGTH == 0.5


def loo_standardised(gp):
    """standardised block-LOO residuals u_i of the REAL observations of a fitted frozen GPv2 (pseudo-image removed together)."""
    g = gp.gp_; X = gp.X_train_; n = len(gp.d_obs_); N = len(X)
    Q = cho_solve((g.L_, True), np.eye(N)); a = g.alpha_.ravel(); u = np.empty(n)
    for i in range(n):
        im = [j for j in range(n, N) if np.linalg.norm(X[j] + X[i]) < 1e-9]; Bk = [i] + im[:1]
        QBB = Q[np.ix_(Bk, Bk)]; CB = np.linalg.inv(QBB); eB = CB @ a[Bk]
        u[i] = eB[0] / np.sqrt(max(CB[0, 0], 1e-300))
    return u


def risk(gp, cand, u):
    """rho_n over candidate axes (normalised to mean 1 over the observations' u^2) and a flag if the V2 fallback was used."""
    u2 = u ** 2; m2 = float(np.mean(u2)) if len(u2) else np.nan
    if not np.all(np.isfinite(u2)) or not np.isfinite(m2) or m2 <= 1e-300:
        return np.ones(len(cand)), True
    k = Matern(RISK_LENGTH, "fixed", nu=1.5); D = gp.d_obs_
    c = np.maximum(k(cand, D), k(cand, -D)); rho = (c @ u2) / c.sum(1) / m2
    return np.maximum(rho, 1e-12), False


class V3Controller(B.V2Controller):
    """frozen V2Controller with the LRV score.  Sees only the candidate axes and its own past observations."""

    def __init__(self, lib, oracle, settings=None, uniform_risk=False):
        super().__init__(lib, oracle, settings); self.uniform_risk = bool(uniform_risk)

    def run(self, start=0):
        s = self.s; seq, obs, steps, scores, unavailable = [], [], [], [], set()
        pending = self.initial_design(start)
        while pending:                                                               # identical to V2Controller.run
            for c, o in zip(pending, self.oracle.evaluate_many(self.C[pending])):
                if self._ok(o):
                    seq.append(c); obs.append(o)
                else:
                    unavailable.add(c); steps.append(dict(test_number=len(seq) + 1, event="failed initial test", candidate=c))
            if len(seq) >= s.n_initial:
                break
            ext = axial_maximin(self.C, seq, s.n_initial, s.maximin_tol_rad, unavailable)
            pending = [i for i in ext if i not in seq]
        while len(seq) < s.budget:
            gp = self.fit(self.C[seq], [o.yield_radius for o in obs])
            mean, var = gp.predict_log(self.C)
            if self.uniform_risk:
                u = np.ones(len(seq)); rho, fb = np.ones(len(self.C)), False
            else:
                u = loo_standardised(gp); rho, fb = risk(gp, self.C, u)
            sc = var * rho; sc[seq] = -np.inf; sc[list(unavailable)] = -np.inf
            v2 = var.copy(); v2[seq] = -np.inf; v2[list(unavailable)] = -np.inf
            j, tied = select(sc, "deterministic", s.tie_rel, s.tie_abs); jm, _ = select(v2, "deterministic", s.tie_rel, s.tie_abs); hp = gp.hyperparameters()
            step = dict(test_number=len(seq) + 1, candidate=int(j), mes_candidate=int(jm), differs_from_mes=bool(j != jm), latent_var=float(var[j]),
                        rho=float(rho[j]), rho_min=float(rho[np.isfinite(sc)].min()), rho_max=float(rho[np.isfinite(sc)].max()), u_abs_max=float(np.abs(u).max()),
                        u_rms=float(np.sqrt(np.mean(u ** 2))), fallback=bool(fb), score=float(sc[j]), entropy=float(entropy(var[j])),
                        gp_mean_log_radius=float(mean[j]), tie_set_size=int(len(tied)), **{k: v for k, v in hp.items()})
            if s.record_scores:
                scores.append(sc.copy())
            o = self.oracle.evaluate(self.C[j])
            if not self._ok(o):
                unavailable.add(j); step["event"] = "failed test (axis marked unavailable)"; steps.append(step); continue
            seq.append(int(j)); obs.append(o); steps.append(step)
        return dict(start=start, sequence=seq, observations=obs, steps=steps, scores=scores, settings=asdict(s), version=VERSION, uniform_risk=self.uniform_risk)
