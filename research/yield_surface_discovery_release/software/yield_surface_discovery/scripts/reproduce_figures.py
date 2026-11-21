#!/usr/bin/env python3
"""Rebuild every release figure and table from the bundled frozen results (data/frozen_results); nothing is re-run.
    python scripts/reproduce_figures.py --results data/frozen_results --output ./figures --tables ./tables"""
import os, sys
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")                     # single-threaded BLAS: bit-reproducible results
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT / "src"))
import argparse
from yield_surface_discovery import figures, tables

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", default=str(ROOT / "data" / "frozen_results")); ap.add_argument("--output", default=str(ROOT / "outputs" / "figures"))
    ap.add_argument("--tables", default=None); ap.add_argument("--only", nargs="*", help="figure ids, e.g. fig1 fig3 figS4")
    a = ap.parse_args()
    for p in figures.make_all(a.results, a.output, a.only):
        print("figure:", p)
    if a.tables:
        for p in tables.make_all(a.results, a.tables):
            print("table:", p)
