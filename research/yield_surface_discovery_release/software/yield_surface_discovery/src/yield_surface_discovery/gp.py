"""Direct Matern-3/2 Gaussian process on log yield radius (the complete surface model; no analytical template, no residual).

y(d) = log r(d),  mean m = median of observed y (a generic constant),  r_hat(d) = exp(m + mu_GP(d)).
Kernel: ConstantKernel(amp) x Matern(nu = 1.5) + WhiteKernel.  Frozen hyper-parameter schedule:
  training set < 20 points (pseudo-observations count): amp = max(var(y), 1e-4), length scale 0.5 and white noise 1e-3 FIXED (no optimiser);
  training set >= 20 points: ML-II over amp [1e-5, 10], length scale [0.05, 5], white noise [1e-6, 0.1] (2 restarts, random_state 0).
Jitter 1e-8 (+ pseudo_noise^2 on symmetry pseudo-observations).  Observations at identical stress points (rounded to 1e-9) are merged.
This reproduces the registered implementation exactly (see tests/test_reproduction.py)."""
import warnings
import numpy as np
from .symmetry import NO_SYMMETRY


class DirectGP:
    def __init__(self, symmetry=NO_SYMMETRY, nu=1.5, length_scale=0.5, length_scale_bounds=(0.05, 5.0), amplitude_bounds=(1e-5, 1e1),
                 white_noise=1e-3, white_noise_bounds=(1e-6, 1e-1), ml2_min_points=20, n_restarts=2, random_state=0, jitter=1e-8,
                 optimize_hyperparameters=True, dedup_decimals=9):
        self.symmetry = symmetry; self.nu = nu; self.length_scale = length_scale; self.length_scale_bounds = tuple(length_scale_bounds)
        self.amplitude_bounds = tuple(amplitude_bounds); self.white_noise = white_noise; self.white_noise_bounds = tuple(white_noise_bounds)
        self.ml2_min_points = ml2_min_points; self.n_restarts = n_restarts; self.random_state = random_state; self.jitter = jitter
        self.optimize_hyperparameters = optimize_hyperparameters; self.dedup_decimals = dedup_decimals; self.gp_ = None

    def fit(self, directions, radii):
        from sklearn.gaussian_process import GaussianProcessRegressor
        from sklearn.gaussian_process.kernels import Matern, WhiteKernel, ConstantKernel
        D = np.atleast_2d(np.asarray(directions, float)); r = np.asarray(radii, float)
        if np.any(~np.isfinite(r)) or np.any(r <= 0):
            raise ValueError("yield radii must be finite and positive")
        P = D * r[:, None]; _, ui = np.unique(np.round(P, self.dedup_decimals), axis=0, return_index=True); ui = np.sort(ui)
        self.d_obs_, self.r_obs_ = D[ui], r[ui]; y = np.log(self.r_obs_); self.mean_ = float(np.median(y)); yc = y - self.mean_
        X, Y, alpha = self.symmetry.augment(self.d_obs_, yc); self.X_train_, self.y_train_ = X, Y
        amp = max(float(np.var(Y)), 1e-4); a = self.jitter + (0 if alpha is None else alpha)
        if (not self.optimize_hyperparameters) or len(X) < self.ml2_min_points:
            k = ConstantKernel(amp, "fixed") * Matern(self.length_scale, "fixed", nu=self.nu) + WhiteKernel(self.white_noise, "fixed")
            gp = GaussianProcessRegressor(k, alpha=a, optimizer=None); self.optimized_ = False
        else:
            k = ConstantKernel(amp, self.amplitude_bounds) * Matern(self.length_scale, self.length_scale_bounds, nu=self.nu) + WhiteKernel(self.white_noise, self.white_noise_bounds)
            gp = GaussianProcessRegressor(k, alpha=a, n_restarts_optimizer=self.n_restarts, random_state=self.random_state); self.optimized_ = True
        with warnings.catch_warnings():
            warnings.simplefilter("ignore"); gp.fit(X, Y)
        self.gp_ = gp; return self

    def predict_log(self, directions, return_var=True):
        """mean of log r (including the constant m) and predictive variance (latent + white noise) at unit directions."""
        mu, sd = self.gp_.predict(np.atleast_2d(directions), return_std=True)
        return (self.mean_ + mu, sd ** 2) if return_var else self.mean_ + mu

    def predict_radius(self, directions):
        return np.exp(self.mean_ + self.gp_.predict(np.atleast_2d(directions)))

    def hyperparameters(self):
        p = self.gp_.kernel_.get_params(); out = dict(optimized=self.optimized_, n_training_points=int(len(self.X_train_)), n_real_observations=int(len(self.d_obs_)))
        for k, v in p.items():
            if k.endswith("constant_value"):
                out["amplitude"] = float(v)
            elif k.endswith("length_scale") and not k.endswith("bounds"):
                out["length_scale"] = float(np.atleast_1d(v)[0])
            elif k.endswith("noise_level") and not k.endswith("bounds"):
                out["white_noise"] = float(v)
        out["kernel"] = str(self.gp_.kernel_); return out
