"""MOOSE pilot worker (moose env): the FROZEN oracle biphasic_oracle_r1.run_trial, unmodified, with three module-level names rebound at
run time: mesh_file -> complexity_scaling_2026/specimens/<config>/mesh_h<nel>/rve.msh, RUN_ROOT -> complexity_scaling_2026/runs, and the
material table entry "CSB1" = the frozen B1 phase properties (materials_frozen.json).  Same solver, criteria, tolerances and Mandel convention.
python pilot_worker.py --config P2_single --direction z1 z2 z3 --tag t --nel 40"""
import argparse, json, math, sys, time, traceback
from pathlib import Path
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent; OR = HERE.parent / "biphasic_benchmark/oracle"; sys.path.insert(0, str(OR))
MARK = "@@OBSERVATION@@"
ap = argparse.ArgumentParser(); ap.add_argument("--config", required=True); ap.add_argument("--direction", nargs=3, type=float, required=True)
ap.add_argument("--tag", required=True); ap.add_argument("--nel", type=int, default=40); a = ap.parse_args(); t0 = time.perf_counter()
try:
    import biphasic_oracle_r1 as O
    O.MATERIALS = dict(O.MATERIALS, CSB1=O.MATERIALS["B1"])
    O.mesh_file = lambda nel: HERE / "specimens" / a.config / f"mesh_h{nel}" / "rve.msh"
    O.RUN_ROOT = HERE / "runs" / a.config
    rec = O.run_trial(a.direction, "CSB1", a.tag, nel=a.nel, dr_rel=0.05, pin_mode="corners"); r = rec["yield_radius"]
    obs = dict(trial_id=f"{a.config}/{a.tag}", yield_radius=r, stress_mandel=rec["stress_mandel"], accepted=rec["accepted"], converged=rec["converged"],
               flags=rec["flags"], runtime_seconds=rec["runtime_seconds"], first_local_yield=rec["first_local_yield_radius_probe"], at_yield=rec["at_yield"], error=None)
except Exception as e:
    obs = dict(trial_id=f"{a.config}/{a.tag}", yield_radius=None, accepted=False, converged=False, runtime_seconds=time.perf_counter() - t0,
               error=f"{type(e).__name__}: {e}", traceback=traceback.format_exc())
clean = lambda v: (None if isinstance(v, float) and math.isnan(v) else v)
print(MARK + json.dumps({k: clean(v) for k, v in obs.items()}, default=float), flush=True)
