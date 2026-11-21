"""Simulation oracle (simulator side): run_moose_trial(direction, specimen_id).

Each trial starts from the virgin specimen (a fresh MOOSE process) and loads it proportionally in macroscopic stress space along the requested Mandel
direction, sigma_M = r * d, with periodic BCs and exact stress control through the
homogenization constraints. It runs three MOOSE solves:

  1. probe : one elastic step at r = 1 MPa -> directional elastic compliance and
             first-local-yield radius r1 = sy / max_qp(vm_matrix) (r = 1 MPa).
  2. ramp  : 1 MPa, r1 (elastic jump to the onset of local yield, where eqps(r) has
             its kink), then geometric steps r1 * (1 + dr_rel)^k until the
             matrix-averaged eqps exceeds `stop_eqps` (Terminator); a checkpoint is
             written at every step.
  3. refine: restart from the ramp checkpoint at the lower end of the bracket of the
             primary criterion (same state, same load path) and take `n_refine`
             sub-steps inside the bracket; the primary yield radius is interpolated
             inside the refined bracket.

The matrix yield stress is used only to plan load steps (step 2), never to compute
the yield radius, and nothing here is returned to the BO side except through
psbench.observe.
"""

import json
import math
import os
import shutil
import subprocess
import time
import uuid
from pathlib import Path

import numpy as np
import pandas as pd

from . import yield_extraction as yx
from .ground_truth import ROOT, build_specimen
from .moose_input import write_input

MOOSE_EXE = os.environ.get("PSBENCH_MOOSE_EXE",
                           str(Path.home() / "projects/moose/modules/combined/combined-opt"))
CONDA_ENV = os.environ.get("PSBENCH_CONDA_ENV", "moose")
RUN_DIR = ROOT / "runs"
SQRT2 = math.sqrt(2.0)

# Acceptance tolerances for a trial
DIRECTION_TOL_DEG = 0.5  # realized vs requested macroscopic stress direction
# max_qp |sigma_zz| relative to the stress radius; the nonlinear tolerance (1e-9 abs on
# residuals of order 1e1-1e2) bounds the out-of-plane residual far below this
SIGMA_ZZ_REL_TOL = 1e-6


def _command(input_file, nprocs):
    exe = [MOOSE_EXE, "-i", input_file, "--color", "off"]
    if nprocs > 1:
        exe = ["mpiexec", "-n", str(nprocs)] + exe
    if Path(os.environ.get("CONDA_PREFIX", "")).name != CONDA_ENV:
        mm = shutil.which("micromamba") or os.environ.get("MAMBA_EXE", "micromamba")
        exe = [mm, "run", "-n", CONDA_ENV] + exe
    return exe


def _run(workdir, name, mesh, spec, d, times, stop_eqps, nprocs, **kw):
    inp = f"{name}.i"
    write_input(workdir / inp, str(mesh), spec, d, times, stop_eqps, **kw)
    t0 = time.perf_counter()
    res = subprocess.run(_command(inp, nprocs), cwd=workdir, capture_output=True, text=True)
    wall = time.perf_counter() - t0
    log = res.stdout + res.stderr
    (workdir / f"{name}.log").write_text(log)
    csv = workdir / f"{name}_out.csv"
    df = pd.read_csv(csv) if csv.exists() else None
    info = dict(
        returncode=res.returncode,
        wall_s=wall,
        failed_solves=log.count("Solve Did NOT Converge") + log.count("Solve Failed"),
        terminated=("Terminator" in log and "eqps_m_avg" in log) or (
            df is not None and len(df) and df["eqps_m_avg"].iloc[-1] > stop_eqps),
    )
    return df, info


def _curve(df):
    """Extract arrays from a postprocessor CSV (drop the t = 0 initial row)."""
    df = df[df["time"] > 0].reset_index(drop=True)
    r = df["time"].to_numpy()
    sig = df[["s_xx", "s_yy", "s_xy"]].to_numpy()
    sig_M = np.c_[sig[:, 0], sig[:, 1], SQRT2 * sig[:, 2]]
    eps = df[["e_xx", "e_yy", "e_xy", "e_zz"]].to_numpy()
    return df, r, sig_M, eps


def run_moose_trial(direction, specimen_id, nel=40, dr_rel=0.05, n_refine=5, refine=True,
                    stop_eqps=4e-3, nprocs=1, tag=None, field_output=True, exodus=False,
                    max_steps=40, keep_checkpoints=False):
    """Run one virtual yield experiment. Returns the full (simulator-side) record."""
    t_start = time.perf_counter()
    d = np.asarray(direction, float)
    d = d / np.linalg.norm(d)
    mesh, mstats, spec = build_specimen(specimen_id, nel)
    tag = tag or uuid.uuid4().hex[:8]
    workdir = RUN_DIR / specimen_id / f"h{nel}" / tag
    workdir.mkdir(parents=True, exist_ok=True)
    mesh = mesh.resolve()
    runs = {}

    # 1. elastic probe
    dfp, runs["probe"] = _run(workdir, "probe", mesh, spec, d, [0.0, 1.0], stop_eqps, nprocs)
    if dfp is None or runs["probe"]["returncode"] != 0:
        raise RuntimeError(f"probe failed, see {workdir}/probe.log")
    p = dfp[dfp["time"] == 1.0].iloc[0]
    if p["eqps_m_max"] > 0:
        raise RuntimeError("elastic probe at r = 1 MPa is already plastic")
    compliance = p[["e_xx", "e_yy", "e_xy", "e_zz"]].to_numpy(float)
    r1 = spec["matrix"]["sy"] / p["vm_m_max"]  # first local (quadrature-point) yield radius

    # 2. ramp
    ramp = [0.0, 1.0] + [r1 * (1.0 + dr_rel) ** k for k in range(0, max_steps)]
    dfr, runs["ramp"] = _run(workdir, "ramp", mesh, spec, d, ramp, stop_eqps, nprocs,
                             field_output=field_output and not refine,
                             exodus=exodus and not refine, checkpoint=refine)
    if dfr is None:
        raise RuntimeError(f"ramp failed, see {workdir}/ramp.log")
    dfr, r, sig_M, eps = _curve(dfr)
    eqps = dfr["eqps_m_avg"].to_numpy()

    thresholds = {f"{t:g}": yx.first_crossing(r, eqps, t)[0] for t in yx.EQPS_THRESHOLDS}
    ep_macro = yx.macro_plastic_strain_eq(r, eps, compliance)
    r_offset = yx.first_crossing(r, ep_macro, yx.OFFSET)[0]
    r_y, k_lo = yx.first_crossing(r, eqps, yx.PRIMARY_EQPS)

    # 3. refinement of the primary bracket
    r_y_coarse = r_y
    df_y, r_c, sig_c, eps_c = dfr, r, sig_M, eps
    if refine and k_lo is not None and k_lo >= 1:
        r_lo, r_hi = r[k_lo], r[k_lo + 1]
        seq = [0.0] + [r_lo + (r_hi - r_lo) * j / n_refine for j in range(1, n_refine + 1)]
        # CSV row k (t > 0) is time step k + 1, which is checkpoint file k + 1
        cp = workdir / "ramp_cp_cp" / f"{k_lo + 1:04d}"
        dff, runs["refine"] = _run(workdir, "refine", mesh, spec, d, seq, yx.PRIMARY_EQPS * 1.0000001,
                                   nprocs, field_output=field_output, exodus=exodus,
                                   restart_from=str(cp))
        if dff is not None and runs["refine"]["returncode"] == 0:
            dff = dff[dff["time"] > 0]
            df_y, r_c, sig_c, eps_c = _curve(pd.concat([dfr.iloc[: k_lo + 1], dff], ignore_index=True))
            r_y, _ = yx.first_crossing(r_c, df_y["eqps_m_avg"].to_numpy(), yx.PRIMARY_EQPS)
    if not keep_checkpoints:
        shutil.rmtree(workdir / "ramp_cp_cp", ignore_errors=True)

    # Quantities at the primary yield point (interpolated in the refined bracket)
    sig_y = yx.interp_at(r_c, sig_c, r_y)
    eps_y = yx.interp_at(r_c, eps_c, r_y)
    dir_real = sig_y / np.linalg.norm(sig_y) if np.all(np.isfinite(sig_y)) else np.full(3, np.nan)
    dir_err = float(np.degrees(np.arccos(np.clip(dir_real @ d, -1.0, 1.0))))
    # plane-stress diagnostics over every recorded step of every run
    allpp = pd.concat([x for x in (dfr, df_y) if x is not None])
    szz_qp = float(np.max(np.abs(allpp[["szz_qp_max", "szz_qp_min"]].to_numpy())))
    sxz_qp = float(np.max(np.abs(allpp[["sxz_qp_max", "sxz_qp_min", "syz_qp_max", "syz_qp_min"]].to_numpy())))
    szz_rel = szz_qp / max(r.max(), 1.0)
    converged = all(v["returncode"] == 0 for v in runs.values()) and np.isfinite(r_y)
    flags = []
    if not np.isfinite(r_y):
        flags.append("yield_not_reached")
    if dir_err > DIRECTION_TOL_DEG:
        flags.append("direction_error")
    if szz_rel > SIGMA_ZZ_REL_TOL:
        flags.append("plane_stress_violation")
    if any(v["failed_solves"] for v in runs.values()):
        flags.append("cutback_occurred")

    def at(col):
        return yx.interp_at(r_c, df_y[col].to_numpy(), r_y) if col in df_y else None

    rec = dict(
        specimen_id=specimen_id,
        tag=tag,
        nel=nel,
        direction=d.tolist(),
        yield_radius=r_y,
        yield_radius_coarse=r_y_coarse,
        first_local_yield_radius=float(r1),
        threshold_radii=thresholds,
        offset_radius=r_offset,
        stress_mandel=sig_y.tolist(),
        macro_strain=[float(eps_y[0]), float(eps_y[1]), float(SQRT2 * eps_y[2])]
        if np.all(np.isfinite(eps_y)) else None,
        macro_strain_tensor=dict(xx=float(eps_y[0]), yy=float(eps_y[1]), xy=float(eps_y[2]), zz=float(eps_y[3])),
        sigma_zz_avg=at("s_zz"),
        sigma_zz_qp_absmax=szz_qp,
        sigma_zz_rel=szz_rel,
        sigma_out_of_plane_shear_qp_absmax=sxz_qp,
        mean_matrix_eqps=at("eqps_m_avg"),
        max_matrix_eqps=at("eqps_m_max"),
        max_vm=at("vm_max"),
        max_vm_inclusion=at("vm_i_max"),
        mean_vm_inclusion=at("vm_i_avg"),
        elastic_compliance=compliance.tolist(),
        direction_error=dir_err,
        converged=bool(converged),
        flags=flags,
        accepted=bool(converged and not ({"direction_error", "plane_stress_violation"} & set(flags))),
        nonlinear_iterations=int(dfr["nl_its"].sum() + (df_y["nl_its"].sum() if df_y is not dfr else 0)),
        n_steps=int(len(dfr) + (len(df_y) if df_y is not dfr else 0)),
        runs=runs,
        runtime_seconds=time.perf_counter() - t_start,
        mesh=mstats,
        workdir=str(workdir),
    )
    (workdir / "record.json").write_text(json.dumps(rec, indent=2, default=float))
    return rec
