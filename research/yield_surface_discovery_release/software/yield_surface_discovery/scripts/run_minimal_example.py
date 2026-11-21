#!/usr/bin/env python3
"""Minimal bundled analytical example (about 10 s): physics-informed controller, expanded library, 12 real tests, Mg-4%Li analytical oracle.
Writes the selected tests, observations, GP reconstruction, a metric and one figure, then (with --check) compares with the stored expected
result data/minimal_example/expected.json.
    python scripts/run_minimal_example.py --output ./outputs/minimal_example --check"""
import os, sys
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")                     # single-threaded BLAS: bit-reproducible results
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT / "src"))
import argparse, json
import numpy as np, pandas as pd
from yield_surface_discovery.workflows import run_config, plot_run

CONFIG = ROOT / "data" / "minimal_example" / "config.yaml"; EXPECTED = ROOT / "data" / "minimal_example" / "expected.json"


def run(output, check=True, write_expected=False):
    res = run_config(CONFIG, output, quiet=True); t = res["tests"]; pr = res["per_run"]
    fig = plot_run(res["output"], 0, title="minimal example: physics-informed controller, 12 real tests (Mg-4%Li analytical oracle)")
    got = dict(labels=t.label.tolist(), global_p95=float(pr[pr.n_real_tests == 12].global_p95.iloc[0]),
               weak10_max_overprediction=float(pr[pr.n_real_tests == 12].weak10_max_overprediction.iloc[0]))
    print("selected tests:", " ".join(got["labels"])); print(f"global p95 relative error at 12 tests: {100 * got['global_p95']:.4f} %"); print("figure:", fig)
    if write_expected:
        EXPECTED.write_text(json.dumps(dict(**got, tolerance_rel=1e-8, note="continuous GP reconstruction scored on the 3183 held-out directions"), indent=1))
    if check:
        exp = json.loads(EXPECTED.read_text()); ok = got["labels"] == exp["labels"] and abs(got["global_p95"] / exp["global_p95"] - 1) <= exp["tolerance_rel"]
        print("CHECK vs expected:", "PASS" if ok else f"FAIL (expected {exp['labels']} / {exp['global_p95']})"); return ok, got
    return True, got


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--output", default=str(ROOT / "outputs" / "minimal_example")); ap.add_argument("--check", action="store_true")
    ap.add_argument("--write-expected", action="store_true", help=argparse.SUPPRESS)
    a = ap.parse_args(); ok, _ = run(a.output, a.check, a.write_expected); sys.exit(0 if ok else 1)
