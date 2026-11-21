#!/usr/bin/env python3
"""Deterministic interface tests for any oracle defined in a config (run BEFORE a closed loop on a new solver, e.g. MOOSE).
Evaluates uniaxial tension, uniaxial compression, equibiaxial tension, pure shear and combined normal + shear; reports radius, achieved-vs-
requested direction, status and run metadata; optionally compares with an analytical cross-check (e.g. --reference von_mises --sigma-y 250).
    python scripts/check_solver_interface.py --config configs/example_moose.yaml --reference von_mises --sigma-y 250 --tol 1e-3"""
import os, sys
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT / "src"))
import argparse
import numpy as np, pandas as pd
from yield_surface_discovery.io import load_config, build_oracle
from yield_surface_discovery.geometry import cauchy_to_mandel, unit, angle_rad
from yield_surface_discovery.analytical import MODELS

CASES = {"uniaxial tension (x)": [1, 0, 0], "uniaxial compression (x)": [-1, 0, 0], "equibiaxial tension": [1, 1, 0], "pure shear (s_xy)": [0, 0, 1],
         "combined normal + shear": [0.6, 0, 0.8]}            # Cauchy components (s_xx, s_yy, s_xy) -> unit Mandel direction

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True); ap.add_argument("--reference", choices=sorted(MODELS)); ap.add_argument("--sigma-y", type=float, default=1.0)
    ap.add_argument("--tol", type=float, default=1e-3, help="relative tolerance for the analytical cross-check"); ap.add_argument("--output")
    a = ap.parse_args(); cfg = load_config(a.config); oracle = build_oracle(cfg); rows, ok = [], True
    for name, s in CASES.items():
        d = unit(cauchy_to_mandel(np.array(s, float))); o = oracle.evaluate(d)
        r = dict(case=name, z1=d[0], z2=d[1], z3=d[2], status=o.status, success=o.success, yield_radius=o.yield_radius, s_xx=o.stress_at_yield[0], s_yy=o.stress_at_yield[1],
                 s_xy=o.stress_at_yield[2], achieved_vs_requested_deg=float(np.degrees(angle_rad(d, o.achieved_direction)[0, 0])) if o.success else np.nan,
                 n_runs=o.metadata.get("n_runs"), run_dir=o.metadata.get("test_dir", o.metadata.get("run_dir")))
        if a.reference and o.success:
            ref = float(MODELS[a.reference](d[None], sigma_y=a.sigma_y)[0]); r.update(reference_radius=ref, rel_error=o.yield_radius / ref - 1); ok &= abs(r["rel_error"]) <= a.tol
        ok &= bool(o.success); rows.append(r)
    df = pd.DataFrame(rows); out = Path(a.output or Path(cfg["output"]["directory"]) / "interface_tests.csv"); out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False); pd.set_option("display.width", 200); print(df.drop(columns=["run_dir"]).to_string(index=False)); print("written:", out)
    print("INTERFACE CHECK", "PASSED" if ok else "FAILED"); sys.exit(0 if ok else 1)
