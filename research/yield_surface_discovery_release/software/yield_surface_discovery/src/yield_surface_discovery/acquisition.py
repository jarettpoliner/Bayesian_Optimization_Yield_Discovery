"""Acquisition: maximum predictive-uncertainty (entropy) sampling from the SAME GP that reconstructs the surface.
score(c) = predictive variance of log r at candidate c (monotone in the Gaussian predictive entropy).  This is NOT max-value entropy search."""
import numpy as np


def predictive_variance_scores(gp, directions, unavailable_mask):
    mean, var = gp.predict_log(directions); sc = var.copy(); sc[unavailable_mask] = -np.inf; return sc, mean, var


def select(scores, ties="deterministic", rel_tol=1e-10, abs_tol=1e-12):
    """returns (chosen index, tie set).  Tie set = {c: score(c) >= max - (|max| rel_tol + abs_tol)}.
    ties='deterministic': lowest index in the tie set (recommended; reproducible across platforms).
    ties='legacy': np.argmax (round-off decides among tied candidates; used only to reproduce the registered runs)."""
    m = np.max(scores)
    if not np.isfinite(m):
        raise RuntimeError("no available candidate left")
    tied = np.where(scores >= m - (abs(m) * rel_tol + abs_tol))[0]
    j = int(np.argmax(scores)) if ties == "legacy" else int(tied.min())
    return j, tied
