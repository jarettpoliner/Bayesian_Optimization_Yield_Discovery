"""Yield Surface BO — a portable MOOSE-connected Bayesian experimental-design framework with stable uncertainty-driven acquisition (V2 BO+MES,
default) and an independently validated analytical error-weighted acquisition method (V3).  Public API:

    from yield_surface_bo import run_campaign, load_config, METHODS, check_compat, load_library, library_report, fit_gp, predict
    res = run_campaign(load_config("examples/analytic_v3.yaml"))

Methods: v2_mes (validated; default for live MOOSE), v3_weighted (research; analytic H1 advantage, physical superiority not established),
v3_signed (experimental), angular_maximin and random (baselines).  All selection mathematics is the frozen research code (see _frozen.py)."""
import os as _os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):       # the frozen studies ran with single-threaded BLAS (reproducibility)
    _os.environ.setdefault(_v, "1")
__version__ = "1.0.0"
from .core import METHODS, DEFAULT_LIVE_METHOD, BUILTIN_LIBRARIES, check_compat, load_library, library_report, StopCampaign   # noqa: F401
from .campaign import run_campaign, load_config, fit_gp, predict, evaluate_errors                                          # noqa: F401
