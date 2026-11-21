"""MOOSE adapter (OPTIONAL; MOOSE is never imported or required by the package).

STATUS: implemented and unit-tested against a fake executable that mimics MOOSE's CSV output (tests/test_moose_adapter.py).
It has NOT been executed against a real MOOSE installation.  Object names in the example input template must be checked against the installed
MOOSE revision (e.g. `moose-opt --registry`).

One virtual experiment = homogeneous single element, stress-controlled proportional loading along the requested direction
(applied Cauchy stress = lambda * mandel_to_cauchy(d), so lambda is the Mandel radius), monotonic in lambda:
  1. coarse ramp: n_ramp equal increments up to lambda_cap (rig capacity; independent of the hidden material); if no yield,
     escalate the capacity (x escalation_factor) up to capacity_escalations times, else status 'no_yield_within_capacity';
  2. first-yield observable: volume-averaged equivalent plastic strain column `eqps` > eqps_threshold (default 1e-8);
  3. bisection between the last elastic and first plastic levels until (hi - lo)/hi <= bisection_rel_tol (default 1e-4), each probe a short run
     (lo, mid);  yield radius = (lo + hi)/2;
  4. checks: achieved stress direction (at the last elastic level) within direction_tol_rad of the request and out-of-plane stress below
     out_of_plane_tol * |stress|; otherwise status 'direction_check_failed'.
Template placeholders (string.Template, ${name}; unknown ${...} expressions are left for MOOSE):
  ${sxx} ${syy} ${sxy} ${t_levels} ${t_end} ${run_id} plus every key of `material_parameters` (kept on the oracle side only).
Expected CSV (MOOSE `[Outputs] csv = true`, file <run_id>.csv): columns time, eqps, sxx, syy, sxy, szz, sxz, syz."""
import string, subprocess, time
from pathlib import Path
import numpy as np, pandas as pd
from ..oracle import YieldOracle, YieldObservation
from ..geometry import mandel_to_cauchy, cauchy_to_mandel, unit, angle_rad


class MooseYieldOracle(YieldOracle):
    def __init__(self, executable, input_template, work_dir, material_parameters=None, lambda_cap=2000.0, n_ramp=20, bisection_rel_tol=1e-4,
                 eqps_threshold=1e-8, capacity_escalations=1, escalation_factor=4.0, launcher=(), extra_args=(), timeout_s=600,
                 direction_tol_rad=1e-6, out_of_plane_tol=1e-8, keep_runs=True):
        self.executable = str(executable); self.template = string.Template(Path(input_template).read_text()); self.work_dir = Path(work_dir)
        self.material_parameters = dict(material_parameters or {}); self.lambda_cap = float(lambda_cap); self.n_ramp = int(n_ramp)
        self.tol = float(bisection_rel_tol); self.thr = float(eqps_threshold); self.n_esc = int(capacity_escalations); self.f_esc = float(escalation_factor)
        self.launcher, self.extra_args, self.timeout_s = list(launcher), list(extra_args), timeout_s
        self.direction_tol_rad, self.out_of_plane_tol, self.keep_runs = direction_tol_rad, out_of_plane_tol, keep_runs
        self.work_dir.mkdir(parents=True, exist_ok=True); self.counter = 0

    # -------------------------------------------------------------- one MOOSE execution
    def _run(self, test_dir, tag, s_hat, levels):
        run_id = f"{tag}"; run_dir = test_dir / tag; run_dir.mkdir(parents=True, exist_ok=True)
        text = self.template.safe_substitute(sxx=repr(float(s_hat[0])), syy=repr(float(s_hat[1])), sxy=repr(float(s_hat[2])),
                                             t_levels=" ".join(repr(float(v)) for v in levels), t_end=repr(float(levels[-1])), run_id=run_id,
                                             **{k: repr(v) if isinstance(v, float) else str(v) for k, v in self.material_parameters.items()})
        inp = run_dir / "input.i"; inp.write_text(text)
        cmd = self.launcher + [self.executable, "-i", str(inp)] + self.extra_args; t0 = time.time()
        p = subprocess.run(cmd, cwd=run_dir, capture_output=True, text=True, timeout=self.timeout_s)
        (run_dir / "stdout.log").write_text(p.stdout); (run_dir / "stderr.log").write_text(p.stderr)
        info = dict(tag=tag, levels=[float(v) for v in levels], returncode=p.returncode, wall_s=time.time() - t0, run_dir=str(run_dir))
        csv = run_dir / f"{run_id}.csv"
        if p.returncode != 0 or not csv.exists():
            return None, info
        df = pd.read_csv(csv); missing = {"time", "eqps", "sxx", "syy", "sxy"} - set(df.columns)
        if missing:
            info["error"] = f"CSV lacks columns {sorted(missing)}"; return None, info
        return df, info

    @staticmethod
    def _row(df, t):
        return df.iloc[int(np.argmin(np.abs(df["time"].to_numpy() - t)))]

    def evaluate(self, direction):
        d = unit(np.asarray(direction, float)); s_hat = mandel_to_cauchy(d); self.counter += 1
        test_dir = self.work_dir / f"test_{self.counter:04d}"; runs = []; meta = dict(test_dir=str(test_dir), runs=runs)
        fail = lambda status: YieldObservation(d, d, np.full(3, np.nan), float("nan"), False, status, meta)
        cap, lo, hi, df_lo = self.lambda_cap, None, None, None
        for esc in range(self.n_esc + 1):
            levels = np.linspace(cap / self.n_ramp, cap, self.n_ramp); df, info = self._run(test_dir, f"ramp{esc}", s_hat, levels); runs.append(info)
            if df is None:
                return fail("solver_failed")
            df = df[df["time"] > 0]; plastic = df["eqps"].to_numpy() > self.thr
            if plastic.any():
                k = int(np.argmax(plastic)); hi = float(df["time"].iloc[k]); lo = float(df["time"].iloc[k - 1]) if k > 0 else 0.0
                df_lo = df.iloc[[k - 1]] if k > 0 else None; df_hi = df.iloc[[k]]; break
            cap *= self.f_esc
        if hi is None:
            return fail("no_yield_within_capacity")
        n_probe = 0
        while (hi - lo) / hi > self.tol:
            mid = 0.5 * (lo + hi); levels = [lo, mid] if lo > 0 else [mid]; n_probe += 1
            df, info = self._run(test_dir, f"probe{n_probe:02d}", s_hat, levels); runs.append(info)
            if df is None:
                return fail("solver_failed")
            row = self._row(df, mid)
            if float(row["eqps"]) > self.thr:
                hi, df_hi = mid, df[np.isclose(df["time"], mid)]
            else:
                lo, df_lo = mid, df[np.isclose(df["time"], mid)]
        r = 0.5 * (lo + hi); ref = df_lo if df_lo is not None and len(df_lo) else df_hi
        sig = ref[["sxx", "syy", "sxy"]].to_numpy(float)[0]; z = cauchy_to_mandel(sig); ach = z / np.linalg.norm(z)
        oop = [float(ref[c].iloc[0]) for c in ("szz", "sxz", "syz") if c in ref.columns]
        dev = float(angle_rad(d, ach)[0, 0]); meta.update(lambda_lo=lo, lambda_hi=hi, n_runs=len(runs), achieved_vs_requested_rad=dev,
                                                          max_out_of_plane=float(max(np.abs(oop))) if oop else None)
        stress = mandel_to_cauchy(r * ach)
        if dev > self.direction_tol_rad or (oop and max(np.abs(oop)) > self.out_of_plane_tol * np.linalg.norm(sig)):
            return YieldObservation(d, ach, stress, r, False, "direction_check_failed", meta)
        return YieldObservation(d, ach, stress, r, True, "ok", meta)
