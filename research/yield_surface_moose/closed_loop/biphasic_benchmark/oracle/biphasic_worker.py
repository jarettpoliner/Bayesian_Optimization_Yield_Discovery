"""Simulator-side worker (moose env, own process): run ONE two-phase trial and print the public observation (same pattern as
closed_loop/moose_worker.py, which is not modified).   python biphasic_worker.py --variant B1 --direction z1 z2 z3 --tag t"""
import argparse, json, math, sys, time, traceback
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
MARK = "@@OBSERVATION@@"
ap = argparse.ArgumentParser(); ap.add_argument("--variant", required=True); ap.add_argument("--direction", nargs=3, type=float, required=True)
ap.add_argument("--tag", required=True); a = ap.parse_args(); t0 = time.perf_counter()
try:
    from biphasic_oracle import observe
    obs = observe(a.direction, a.variant, a.tag); obs["error"] = None
except Exception as e:
    obs = dict(trial_id=f"{a.variant}/{a.tag}", yield_radius=None, accepted=False, converged=False, runtime_seconds=time.perf_counter() - t0,
               error=f"{type(e).__name__}: {e}", traceback=traceback.format_exc())
clean = lambda v: (None if isinstance(v, float) and math.isnan(v) else v)
print(MARK + json.dumps({k: clean(v) for k, v in obs.items()}, default=float), flush=True)
