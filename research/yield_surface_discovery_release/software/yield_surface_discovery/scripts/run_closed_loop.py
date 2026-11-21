#!/usr/bin/env python3
"""Run the frozen controller against ANY oracle defined in a config (analytical, external_command, moose).
    python scripts/run_closed_loop.py --config configs/example_solver.yaml
Writes the selected test sequence, every observation (requested/achieved direction, stress at yield, radius, status, metadata), the GP
reconstruction (surface mesh .obj, s_xy = 0 section with uncertainty), hyper-parameters, and metrics when an analytical reference exists."""
import os, sys
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")                     # single-threaded BLAS: bit-reproducible results
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT / "src"))
import argparse
from yield_surface_discovery.workflows import run_config, plot_run

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True); ap.add_argument("--output"); ap.add_argument("--start", type=int, default=None)
    a = ap.parse_args()
    res = run_config(a.config, a.output, None if a.start is None else [a.start], export_mesh=True)
    for s in sorted(set(res["tests"].start)):
        print("figure:", plot_run(res["output"], s, title=f"closed loop, start {s}"))
    print(res["tests"][["start", "test_number", "label", "family", "yield_radius", "status"]].to_string(index=False))
    print("outputs written to", res["output"])
