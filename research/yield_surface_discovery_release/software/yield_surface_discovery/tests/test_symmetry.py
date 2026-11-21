import numpy as np
from _helpers import library
from yield_surface_discovery.symmetry import SHEAR_MIRROR, NO_SYMMETRY, equivalence_classes, redundant_count, MIRROR_MATRIX


def test_mirror_and_augmentation_rule():
    assert np.array_equal(MIRROR_MATRIX, np.diag([1.0, 1.0, -1.0]))
    d = np.array([[1.0, 0, 0], [0, 0.6, 0.8]]); X, y, a = SHEAR_MIRROR.augment(d, np.array([0.1, -0.1]))
    assert len(X) == 3 and np.allclose(X[2], [0, 0.6, -0.8]) and y[2] == -0.1 and a[2] == 1e-8      # self-mirror direction adds nothing
    X, y, a = NO_SYMMETRY.augment(d, np.array([0.1, -0.1])); assert a is None and len(X) == 2


def test_expanded_library_deduplication_and_classes():
    lib = library(dedup=True); assert len(lib) == 414
    assert sorted((r[1], r[3]) for r in lib.removed_duplicates) == [("shear+160", "shear-70"), ("shear-160", "shear+70")]
    cls = equivalence_classes(lib.directions); assert len(set(cls.tolist())) == 243


def test_redundancy_count():
    d = np.array([[0, 0, 1.0], [0, 0, -1.0], [1.0, 0, 0]]); assert redundant_count(d) == 1 and redundant_count(d, NO_SYMMETRY) == 0
