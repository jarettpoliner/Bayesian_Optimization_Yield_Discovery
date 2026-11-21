"""BO-facing observation interface.

This is the only module the BO side should import. It exposes measured, macroscopic
quantities of a virtual experiment and nothing about the specimen definition
(phase properties, microstructure, mesh) or the oracle internals.
"""

import json
import math
import time

from .oracle import RUN_DIR, run_moose_trial

# Measured quantities the BO side may see. Everything else in the oracle record
# (phase properties, geometry, mesh, paths, step-planning radii) stays simulator side.
PUBLIC_KEYS = (
    "trial_id",
    "yield_radius",
    "log_yield_radius",
    "stress_mandel",
    "macro_strain",
    "mean_matrix_eqps",
    "sigma_zz",
    "direction_error",
    "converged",
    "accepted",
    "runtime_seconds",
)

LOG = RUN_DIR / "observations.jsonl"


def observe(direction, specimen_id, **oracle_options):
    """Run one yield experiment along the Mandel stress direction `direction`
    (sigma_xx, sigma_yy, sqrt(2) tau_xy; normalized internally) on the fixed specimen
    `specimen_id`, and return only public measured quantities.

    `accepted` is False when the run did not converge or failed the realized-direction
    or plane-stress checks; such observations should not be used to update the GP.
    """
    rec = run_moose_trial(direction, specimen_id, **oracle_options)
    r = rec["yield_radius"]
    obs = dict(
        trial_id=f"{specimen_id}/{rec['tag']}",
        yield_radius=r,
        log_yield_radius=math.log(r) if r and math.isfinite(r) and r > 0 else float("nan"),
        stress_mandel=rec["stress_mandel"],
        macro_strain=rec["macro_strain"],
        mean_matrix_eqps=rec["mean_matrix_eqps"],
        sigma_zz=rec["sigma_zz_qp_absmax"],
        direction_error=rec["direction_error"],
        converged=rec["converged"],
        accepted=rec["accepted"],
        runtime_seconds=rec["runtime_seconds"],
    )
    assert set(obs) == set(PUBLIC_KEYS)
    LOG.parent.mkdir(parents=True, exist_ok=True)
    with open(LOG, "a") as f:
        f.write(json.dumps(dict(time=time.time(), specimen_id=specimen_id,
                                direction=list(map(float, direction)), **obs), default=float) + "\n")
    return obs
