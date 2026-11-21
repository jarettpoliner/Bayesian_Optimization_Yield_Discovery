#!/usr/bin/env python3
"""Run the frozen controller against an ANALYTICAL oracle and score it on held-out directions.
    python scripts/run_benchmark.py --config configs/benchmark_physics_informed.yaml
    python scripts/run_benchmark.py --config configs/benchmark_blind.yaml
Outputs (config output.directory or --output): tests.csv, steps.csv, per_run.csv, summary.csv, observations_start*.json, surface/section files, run_info.json."""
import os, sys
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")                     # single-threaded BLAS: bit-reproducible results
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT / "src"))
import argparse
from yield_surface_discovery.workflows import run_config

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True); ap.add_argument("--output"); ap.add_argument("--starts", type=int, nargs="*"); ap.add_argument("--no-mesh", action="store_true")
    a = ap.parse_args()
    res = run_config(a.config, a.output, a.starts, export_mesh=False if a.no_mesh else None)
    print("outputs written to", res["output"])
