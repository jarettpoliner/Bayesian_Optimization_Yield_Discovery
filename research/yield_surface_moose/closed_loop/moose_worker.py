"""Simulator-side worker: run ONE MOOSE yield experiment and print the public observation as JSON.

Runs in the `moose` environment, in its own process, so the BO controller process never imports psbench
(and therefore never loads psbench.ground_truth or the specimen truth files).  The only inputs are the specimen id,
the Mandel loading direction and a run tag; the only output is the dict returned by psbench.observe.observe
(plus an `error` field if the trial raised).

    python moose_worker.py --specimen S2 --direction z1 z2 z3 --tag bo_xxx_t01
"""
import argparse
import json
import math
import sys
import time
import traceback
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))          # plane_stress/ (for psbench)

MARK = "@@OBSERVATION@@"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--specimen", required=True)
    ap.add_argument("--direction", nargs=3, type=float, required=True)
    ap.add_argument("--tag", required=True)
    a = ap.parse_args()
    t0 = time.perf_counter()
    try:
        from psbench.observe import observe
        obs = observe(a.direction, a.specimen, tag=a.tag)
        obs["error"] = None
    except Exception as e:                                             # a failed trial must be reported, never swallowed
        obs = dict(trial_id=f"{a.specimen}/{a.tag}", yield_radius=float("nan"), log_yield_radius=float("nan"),
                   accepted=False, converged=False, runtime_seconds=time.perf_counter() - t0,
                   error=f"{type(e).__name__}: {e}", traceback=traceback.format_exc())
    def clean(v):
        if isinstance(v, float) and not math.isfinite(v):
            return None if math.isnan(v) else str(v)
        if isinstance(v, (list, tuple)):
            return [clean(x) for x in v]
        return v
    print(MARK + json.dumps({k: clean(v) for k, v in obs.items()}, default=float), flush=True)


if __name__ == "__main__":
    main()
