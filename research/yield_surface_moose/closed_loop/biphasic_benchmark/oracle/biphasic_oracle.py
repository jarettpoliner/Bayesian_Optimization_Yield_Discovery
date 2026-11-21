"""Two-plastic-phase simulation oracle (SIMULATOR SIDE).  Separate from psbench (which is imported read-only for nothing but paths):
same three-solve procedure as psbench.oracle.run_moose_trial (elastic probe at r = 1 MPa, geometric ramp with checkpoints, checkpoint-restart
refinement inside the bracket), same formulation, mesh, solver settings, acceptance tolerances and Mandel convention, with these changes:
  * input writer: proposed_inputs/moose_input_biphasic.py (both phases J2 + linear isotropic hardening);
  * step planning: r1 = min over plastic phases of sy_p / max_qp vm_p at the elastic probe (exact first local yield of the first phase);
  * Terminator on the RVE volume-average eqps (eqps_rve_avg > stop_eqps);
  * PRIMARY effective yield: first crossing of eqps_rve_avg = 1e-4 (refined); SECONDARY: eqps_m_avg = 1e-4 (refined separately; this is the
    psbench criterion, used for regression against stored S2 records);
  * extended record: per-step phase-resolved history, volume-consistency check, first local yield per phase, phase shares at yield.
Each trial: fresh MOOSE processes in runs/<variant>/h<nel>/<tag>/, record.json written there."""
import json, math, os, shutil, subprocess, sys, time
from pathlib import Path
import numpy as np, pandas as pd

HERE = Path(__file__).resolve().parent; BB = HERE.parent; PS = BB.parents[1]
sys.path.insert(0, str(BB / "proposed_inputs"))
import moose_input_biphasic as MB                                                   # noqa: E402

MOOSE_EXE = str(Path.home() / "projects/moose/modules/combined/combined-opt")
RUN_ROOT = BB / "runs"; SQRT2 = math.sqrt(2.0)
DIRECTION_TOL_DEG = 0.5; SIGMA_ZZ_REL_TOL = 1e-6                                    # identical to psbench.oracle
CRITERIA = (("eqps_rve_avg", 1e-4), ("eqps_m_avg", 1e-4))                            # primary first
MATERIALS = json.loads((HERE / "materials_frozen.json").read_text())


def mesh_file(nel):
    return PS / f"specimens/S2/mesh_h{nel}/rve.msh"


def _run(workdir, name, mesh, spec, d, times, stop_value, stop_pp, **kw):
    MB.write_input(workdir / f"{name}.i", str(mesh), spec, d, times, stop_value, stop_pp=stop_pp, **kw)
    t0 = time.perf_counter()
    res = subprocess.run([MOOSE_EXE, "-i", f"{name}.i", "--color", "off"], cwd=workdir, capture_output=True, text=True)
    log = res.stdout + res.stderr; (workdir / f"{name}.log").write_text(log)
    csv = workdir / f"{name}_out.csv"; df = pd.read_csv(csv) if csv.exists() else None
    return df, dict(returncode=res.returncode, wall_s=time.perf_counter() - t0, failed_solves=log.count("Solve Did NOT Converge") + log.count("Solve Failed"))


def first_crossing(r, q, level):
    idx = np.nonzero(np.asarray(q) >= level)[0]
    if idx.size == 0:
        return float("nan"), None
    k = int(idx[0])
    if k == 0:
        return float(r[0]), 0
    return float(r[k - 1] + (level - q[k - 1]) / (q[k] - q[k - 1]) * (r[k] - r[k - 1])), k - 1


def interp(r, y, rs):
    return float(np.interp(rs, r, y)) if np.isfinite(rs) else float("nan")


def run_trial(direction, variant, tag, nel=40, dr_rel=0.05, n_refine=5, stop_eqps=4e-3, max_steps=40, field_output=True):
    t_start = time.perf_counter(); d = np.asarray(direction, float); d = d / np.linalg.norm(d)
    spec = MATERIALS[variant]; mesh = mesh_file(nel).resolve(); wd = RUN_ROOT / variant / f"h{nel}" / tag
    if wd.exists():
        raise RuntimeError(f"{wd} exists; refusing to overwrite")
    wd.mkdir(parents=True); runs = {}
    dfp, runs["probe"] = _run(wd, "probe", mesh, spec, d, [0.0, 1.0], stop_eqps, "eqps_rve_avg")
    if dfp is None or runs["probe"]["returncode"] != 0:
        raise RuntimeError(f"probe failed, see {wd}/probe.log")
    p = dfp[dfp["time"] == 1.0].iloc[0]
    if p["eqps_m_max"] > 0 or p["eqps_i_max"] > 0:
        raise RuntimeError("elastic probe at r = 1 MPa is already plastic")
    compliance = p[["e_xx", "e_yy", "e_xy", "e_zz"]].to_numpy(float)
    r1_m = spec["matrix"]["sy"] / p["vm_m_max"]; r1_i = spec["inclusion"]["sy"] / p["vm_i_max"]; r1 = min(r1_m, r1_i)
    ramp = [0.0, 1.0] + [r1 * (1.0 + dr_rel) ** k for k in range(max_steps)]
    dfr, runs["ramp"] = _run(wd, "ramp", mesh, spec, d, ramp, stop_eqps, "eqps_rve_avg", checkpoint=True)
    if dfr is None:
        raise RuntimeError(f"ramp failed, see {wd}/ramp.log")
    dfr = dfr[dfr["time"] > 0].reset_index(drop=True); r = dfr["time"].to_numpy()
    out = {}
    for j, (pp, lev) in enumerate(CRITERIA):
        r_c, k_lo = first_crossing(r, dfr[pp].to_numpy(), lev); cur = dfr; name = f"refine_{pp}"
        if k_lo is not None and k_lo >= 1:
            r_lo, r_hi = r[k_lo], r[k_lo + 1]
            seq = [0.0] + [r_lo + (r_hi - r_lo) * i / n_refine for i in range(1, n_refine + 1)]
            dff, runs[name] = _run(wd, name, mesh, spec, d, seq, lev * 1.0000001, pp, field_output=(field_output and j == 0),
                                   restart_from=str(wd / "ramp_cp_cp" / f"{k_lo + 1:04d}"))
            if dff is not None and runs[name]["returncode"] == 0:
                cur = pd.concat([dfr.iloc[: k_lo + 1], dff[dff["time"] > 0]], ignore_index=True)
                r_ref, _ = first_crossing(cur["time"].to_numpy(), cur[pp].to_numpy(), lev)
            else:
                r_ref = float("nan")
        else:
            r_ref = r_c
        out[pp] = dict(level=lev, r_coarse=r_c, r_refined=r_ref, curve=cur)
    shutil.rmtree(wd / "ramp_cp_cp", ignore_errors=True)
    prim = out[CRITERIA[0][0]]; cr = prim["curve"]; rc = cr["time"].to_numpy(); ry = prim["r_refined"]
    sig = np.c_[cr.s_xx, cr.s_yy, SQRT2 * cr.s_xy]; sig_y = np.array([interp(rc, sig[:, k], ry) for k in range(3)])
    dir_real = sig_y / np.linalg.norm(sig_y) if np.all(np.isfinite(sig_y)) else np.full(3, np.nan)
    dir_err = float(np.degrees(np.arccos(np.clip(dir_real @ d, -1, 1))))
    allpp = pd.concat([dfr] + [o["curve"] for o in out.values()])
    szz = float(np.max(np.abs(allpp[["szz_qp_max", "szz_qp_min"]].to_numpy()))); szz_rel = szz / max(r.max(), 1.0)
    vm_, vi_ = float(p.vol_m), float(p.vol_i); fm, fi = vm_ / (vm_ + vi_), vi_ / (vm_ + vi_)
    cons = np.abs(allpp.eqps_rve_avg - (fm * allpp.eqps_m_avg + fi * allpp.eqps_i_avg)) / np.maximum(allpp.eqps_rve_avg.abs(), 1e-300)
    cons_abs = np.abs(allpp.eqps_rve_avg - (fm * allpp.eqps_m_avg + fi * allpp.eqps_i_avg))
    first_plastic = {ph: (float(r[np.argmax(dfr[f"eqps_{t}_max"].to_numpy() > 0)]) if np.any(dfr[f"eqps_{t}_max"].to_numpy() > 0) else float("nan"))
                     for ph, t in (("matrix", "m"), ("inclusion", "i"))}
    converged = all(v["returncode"] == 0 for v in runs.values()) and np.isfinite(ry)
    flags = [f for f, c in (("yield_not_reached", not np.isfinite(ry)), ("direction_error", dir_err > DIRECTION_TOL_DEG),
                            ("plane_stress_violation", szz_rel > SIGMA_ZZ_REL_TOL), ("cutback_occurred", any(v["failed_solves"] for v in runs.values()))) if c]
    at = lambda col: interp(rc, cr[col].to_numpy(), ry)
    hist_cols = ["time", "s_xx", "s_yy", "s_xy", "e_xx", "e_yy", "e_xy", "e_zz", "eqps_rve_avg", "eqps_m_avg", "eqps_i_avg", "eqps_m_max", "eqps_i_max",
                 "pvf_m", "pvf_i", "wp_rve_avg", "wp_m_avg", "wp_i_avg", "vm_m_max", "vm_i_max", "vm_m_avg", "vm_i_avg", "nl_its"]
    dfr[hist_cols].to_csv(wd / "ramp_history.csv", index=False)
    rec = dict(variant=variant, tag=tag, nel=nel, direction=d.tolist(), yield_radius=ry,
               criteria={pp: dict(level=o["level"], r_coarse=o["r_coarse"], r_refined=o["r_refined"]) for pp, o in out.items()},
               first_local_yield_radius_probe=dict(matrix=float(r1_m), inclusion=float(r1_i)), first_plastic_ramp_step_radius=first_plastic, r1_planning=float(r1),
               stress_mandel=sig_y.tolist(), direction_error=dir_err, sigma_zz_qp_absmax=szz, sigma_zz_rel=szz_rel,
               at_yield=dict(eqps_rve_avg=at("eqps_rve_avg"), eqps_m_avg=at("eqps_m_avg"), eqps_i_avg=at("eqps_i_avg"), pvf_m=at("pvf_m"), pvf_i=at("pvf_i"),
                             wp_rve_avg=at("wp_rve_avg"), e_xx=at("e_xx"), e_yy=at("e_yy"), e_xy=at("e_xy"), e_zz=at("e_zz")),
               end_of_ramp=dict(r=float(r[-1]), **{c: float(dfr[c].iloc[-1]) for c in ("eqps_rve_avg", "eqps_m_avg", "eqps_i_avg", "pvf_m", "pvf_i")}),
               volumes=dict(matrix=vm_, inclusion=vi_, total=vm_ + vi_, f_m=fm, f_i=fi),
               volume_consistency=dict(max_abs=float(cons_abs.max()), max_rel_where_positive=float(cons[allpp.eqps_rve_avg > 0].max()) if (allpp.eqps_rve_avg > 0).any() else 0.0),
               elastic_compliance=compliance.tolist(), converged=bool(converged), flags=flags,
               accepted=bool(converged and not ({"direction_error", "plane_stress_violation"} & set(flags))),
               nonlinear_iterations=int(allpp.nl_its.sum()), n_steps=int(len(allpp)), runs=runs, runtime_seconds=time.perf_counter() - t_start, workdir=str(wd))
    (wd / "record.json").write_text(json.dumps(rec, indent=2, default=float))
    return rec


PUBLIC_KEYS = ("trial_id", "yield_radius", "log_yield_radius", "stress_mandel", "mean_rve_eqps", "direction_error", "converged", "accepted", "runtime_seconds")


def observe(direction, variant, tag):
    rec = run_trial(direction, variant, tag); r = rec["yield_radius"]
    return dict(trial_id=f"{variant}/{tag}", yield_radius=r, log_yield_radius=math.log(r) if r and math.isfinite(r) and r > 0 else float("nan"),
                stress_mandel=rec["stress_mandel"], mean_rve_eqps=rec["at_yield"]["eqps_rve_avg"], direction_error=rec["direction_error"],
                converged=rec["converged"], accepted=rec["accepted"], runtime_seconds=rec["runtime_seconds"])
