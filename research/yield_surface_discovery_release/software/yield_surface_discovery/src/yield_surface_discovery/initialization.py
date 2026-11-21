"""Initial design: angular maximin over the candidate library (geometry only; never uses yield values)."""
import numpy as np
from .geometry import angle_rad, arccos_angle_rad

LEGACY_EXCLUDE = (0, 18, 9, 80)       # candidate ids excluded from the random first draw in the registered runs (Hill-calibration-era rule)


def first_candidate(start, n_candidates, rule="seeded", seed_base=7000, first_index=0, legacy_exclude=LEGACY_EXCLUDE):
    """start 0 -> first_index; start s >= 1 -> numpy default_rng(seed_base + s) draw.
    rule 'seeded': draw over all candidates.  rule 'registered_legacy': draw over candidates minus legacy_exclude (reproduces the registered runs)."""
    if start == 0:
        return int(first_index)
    rng = np.random.default_rng(seed_base + start)
    if rule == "seeded":
        return int(rng.choice(n_candidates))
    if rule == "registered_legacy":
        return int(rng.choice(np.setdiff1d(np.arange(n_candidates), list(legacy_exclude))))
    raise ValueError(f"unknown start rule {rule!r}")


def angular_maximin(D, seeds, n_total=8, ties="deterministic", tol=1e-9, unavailable=()):
    """greedy max-min angle design continuing from `seeds` (an index or a list of already selected indices).
    ties='legacy': exactly the registered implementation (np.argmax on arccos angles);
    ties='deterministic': atan2 angles, candidates within tol (rad) of the best are tied and the lowest index wins."""
    D = np.asarray(D); sel = [int(seeds)] if np.isscalar(seeds) else [int(i) for i in seeds]; bad = set(int(u) for u in unavailable)
    while len(sel) < n_total:
        if ties == "legacy":
            rem = np.setdiff1d(np.arange(len(D)), sel + sorted(bad))
            sel.append(int(rem[np.argmax(arccos_angle_rad(D[rem], D[sel]).min(1))]))
        else:
            dist = angle_rad(D, D[sel]).min(1); dist[sel] = -np.inf; dist[list(bad)] = -np.inf
            m = dist.max(); sel.append(int(np.where(dist >= m - tol)[0].min()))
    return sel
