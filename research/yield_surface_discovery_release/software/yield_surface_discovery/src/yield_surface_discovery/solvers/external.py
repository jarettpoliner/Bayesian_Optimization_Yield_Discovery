"""Generic command-line solver adapter.

For each requested test the oracle creates a run directory, formats `command` (a list of strings) with the placeholders
  {z1} {z2} {z3}            requested unit Mandel direction
  {s_xx} {s_yy} {s_xy}      the same direction as Cauchy components (tensor shear)
  {run_dir} {run_id}        the run directory / identifier
runs it, keeps stdout/stderr, and reads a JSON result (default {run_dir}/result.json) with keys:
  yield_radius (Mandel norm)  OR  stress_at_yield [s_xx, s_yy, s_xy];   optional: achieved_direction [z1, z2, z3], status, metadata.
Non-zero exit codes, time-outs and missing/invalid results return success = False (the controller then marks the candidate unavailable)."""
import json, subprocess, time
from pathlib import Path
import numpy as np
from ..oracle import YieldOracle, YieldObservation
from ..geometry import mandel_to_cauchy, cauchy_to_mandel, unit


class ExternalCommandOracle(YieldOracle):
    def __init__(self, command, work_dir, timeout_s=3600, result_file="result.json", env=None):
        self.command, self.work_dir, self.timeout_s, self.result_file, self.env = list(command), Path(work_dir), timeout_s, result_file, env
        self.work_dir.mkdir(parents=True, exist_ok=True); self.counter = 0

    def evaluate(self, direction):
        d = unit(np.asarray(direction, float)); s = mandel_to_cauchy(d); self.counter += 1
        run_id = f"test_{self.counter:04d}"; run_dir = self.work_dir / run_id; run_dir.mkdir(parents=True, exist_ok=True)
        fmt = dict(z1=repr(float(d[0])), z2=repr(float(d[1])), z3=repr(float(d[2])), s_xx=repr(float(s[0])), s_yy=repr(float(s[1])), s_xy=repr(float(s[2])),
                   run_dir=str(run_dir), run_id=run_id)
        cmd = [c.format(**fmt) for c in self.command]; meta = dict(command=cmd, run_dir=str(run_dir)); t0 = time.time()
        try:
            p = subprocess.run(cmd, cwd=run_dir, capture_output=True, text=True, timeout=self.timeout_s, env=self.env)
            (run_dir / "stdout.log").write_text(p.stdout); (run_dir / "stderr.log").write_text(p.stderr); meta.update(returncode=p.returncode, wall_s=time.time() - t0)
            if p.returncode != 0:
                return YieldObservation(d, d, np.full(3, np.nan), float("nan"), False, "solver_failed", meta)
            res = json.loads((run_dir / self.result_file).read_text())
        except subprocess.TimeoutExpired:
            return YieldObservation(d, d, np.full(3, np.nan), float("nan"), False, "timeout", meta)
        except (OSError, ValueError) as e:
            meta["error"] = repr(e); return YieldObservation(d, d, np.full(3, np.nan), float("nan"), False, "invalid_result", meta)
        if "stress_at_yield" in res:
            S = np.asarray(res["stress_at_yield"], float); z = cauchy_to_mandel(S); r = float(np.linalg.norm(z)); ach = z / r
        else:
            r = float(res["yield_radius"]); ach = unit(np.asarray(res.get("achieved_direction", d), float)); S = mandel_to_cauchy(r * ach)
        status = res.get("status", "ok"); meta.update(res.get("metadata", {}))
        return YieldObservation(d, ach, S, r, bool(status == "ok" and np.isfinite(r) and r > 0), status, meta)
