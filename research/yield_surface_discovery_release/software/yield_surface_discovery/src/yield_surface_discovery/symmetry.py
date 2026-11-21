"""Known material symmetry used as an OPTIONAL physics-informed prior.

Only one symmetry is implemented: the shear mirror s_xy -> -s_xy (Mandel z3 -> -z3), exact for orthotropic materials whose axes are the
test axes.  It enters the GP as noise-1e-4 pseudo-observations at the mirror image of every real observation (frozen rule: a mirror is added iff
it lies > 1e-6 rad from every real observation).  Without it the method is the generic blind baseline."""
from dataclasses import dataclass, field
import numpy as np

MIRROR_MATRIX = np.diag([1.0, 1.0, -1.0])


@dataclass(frozen=True)
class Symmetry:
    name: str
    operations: tuple = field(default_factory=tuple)      # extra 3x3 maps (identity excluded)
    pseudo_noise: float = 1e-4
    duplicate_angle_rad: float = 1e-6

    @property
    def enabled(self):
        return len(self.operations) > 0

    def images(self, D):
        D = np.atleast_2d(D); return [D @ np.asarray(g).T for g in self.operations]

    def augment(self, d_u, yc):
        """frozen augmentation (identical to the registered implementation): returns X, y, alpha_extra (None if no symmetry)."""
        if not self.enabled:
            return d_u, yc, None
        px = np.concatenate(self.images(d_u)); keep = np.arccos(np.clip(px @ d_u.T, -1, 1)).min(1) > self.duplicate_angle_rad
        X = np.vstack([d_u, px[keep]]); y = np.r_[yc, np.tile(yc, len(self.operations))[keep]]
        alpha = np.r_[np.zeros(len(d_u)), np.full(int(keep.sum()), self.pseudo_noise ** 2)]
        return X, y, alpha

    def covered(self, D):
        """real directions plus their symmetry images."""
        D = np.atleast_2d(D); return np.vstack([D] + self.images(D)) if self.enabled else D


NO_SYMMETRY = Symmetry("none", ())
SHEAR_MIRROR = Symmetry("shear_mirror", (MIRROR_MATRIX,))


def symmetry_from_name(name_or_bool):
    if name_or_bool in (True, "shear_mirror", "shear-mirror", "mirror"):
        return SHEAR_MIRROR
    if name_or_bool in (False, None, "none", "off"):
        return NO_SYMMETRY
    raise ValueError(f"unknown symmetry {name_or_bool!r}; only 'shear_mirror' or 'none' are supported")


def equivalence_classes(D, symmetry=SHEAR_MIRROR, tol=1e-9):
    """class id per direction: d ~ d' iff |d - d'| < tol or |g d - d'| < tol for a symmetry operation g."""
    D = np.atleast_2d(D); cls = -np.ones(len(D), int); c = 0
    for i in range(len(D)):
        if cls[i] >= 0:
            continue
        cls[i] = c; imgs = [D[i]] + [np.asarray(g) @ D[i] for g in symmetry.operations]
        for j in range(i + 1, len(D)):
            if cls[j] < 0 and any(np.linalg.norm(im - D[j]) < tol for im in imgs):
                cls[j] = c
        c += 1
    return cls


def redundant_count(D, symmetry=SHEAR_MIRROR, tol=1e-9):
    """number of real tests whose information is already implied (duplicate or symmetry image of another selected test)."""
    return int(len(np.atleast_2d(D)) - len(set(equivalence_classes(D, symmetry, tol).tolist())))
