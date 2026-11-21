"""yield_surface_discovery: adaptive selection of physical loading tests to reconstruct an unknown plane-stress yield surface.

Frozen practical method: 8 angular-maximin initial tests -> direct Matern-3/2 GP on log yield radius -> maximum predictive-uncertainty
sampling from that same GP -> repeat; optional physics-informed augmentation with the known shear-mirror symmetry (s_xy -> -s_xy).
The adaptive controller only talks to a YieldOracle (analytical model, external solver, MOOSE, or a laboratory)."""
__version__ = "1.0.0"
from .geometry import mandel_to_cauchy, cauchy_to_mandel, unit, angle_rad
from .symmetry import Symmetry, NO_SYMMETRY, SHEAR_MIRROR, symmetry_from_name
from .candidates import CandidateLibrary
from .gp import DirectGP
from .oracle import YieldObservation, YieldOracle, AnalyticalOracle, CachedOracle
from .controller import AdaptiveController, ControllerSettings, RunResult
