"""Worker r2 = biphasic_worker_r1.py plus --nel and --dr-rel pass-through to biphasic_oracle_r1.run_trial (unchanged oracle r1).
python biphasic_worker_r2.py --variant B1 --direction z1 z2 z3 --tag t [--pin-mode corners] [--nel 40] [--dr-rel 0.05]"""
import argparse, json, math, sys, time, traceback
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
MARK = "@@OBSERVATION@@"
ap = argparse.ArgumentParser(); ap.add_argument("--variant", required=True); ap.add_argument("--direction", nargs=3, type=float, required=True)
ap.add_argument("--tag", required=True); ap.add_argument("--pin-mode", default="corners", choices=["corners", "origin"])
ap.add_argument("--nel", type=int, default=40); ap.add_argument("--dr-rel", type=float, default=0.05); a = ap.parse_args(); t0 = time.perf_counter()
try:
    from biphasic_oracle_r1 import run_trial
    rec = run_trial(a.direction, a.variant, a.tag, nel=a.nel, dr_rel=a.dr_rel, pin_mode=a.pin_mode); r = rec["yield_radius"]
    obs = dict(trial_id=f"{a.variant}/{a.tag}", yield_radius=r, log_yield_radius=math.log(r) if r and math.isfinite(r) and r > 0 else float("nan"),
               stress_mandel=rec["stress_mandel"], mean_rve_eqps=rec["at_yield"]["eqps_rve_avg"], direction_error=rec["direction_error"], converged=rec["converged"],
               accepted=rec["accepted"], runtime_seconds=rec["runtime_seconds"], error=None)
except Exception as e:
    obs = dict(trial_id=f"{a.variant}/{a.tag}", yield_radius=None, accepted=False, converged=False, runtime_seconds=time.perf_counter() - t0,
               error=f"{type(e).__name__}: {e}", traceback=traceback.format_exc())
clean = lambda v: (None if isinstance(v, float) and math.isnan(v) else v)
print(MARK + json.dumps({k: clean(v) for k, v in obs.items()}, default=float), flush=True)
