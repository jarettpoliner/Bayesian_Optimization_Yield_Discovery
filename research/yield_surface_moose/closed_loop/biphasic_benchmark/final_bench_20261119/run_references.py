#!/usr/bin/env python3
"""Stage 2 (EVALUATOR side): run the frozen 72 reference axes for one variant via biphasic_worker_r2 (validated configuration).
Never read by any controller or sampler.  Output: refs/<variant>/trials.jsonl (all records incl. failures).  python run_references.py B1 [parallel]"""
import hashlib, json, sys
from pathlib import Path
import pandas as pd
FB = Path(__file__).resolve().parent; sys.path.insert(0, str(FB))
from bench_oracle import BiphasicOracle                                             # noqa: E402
v = sys.argv[1]; par = int(sys.argv[2]) if len(sys.argv) > 2 else 11
des = pd.read_csv(FB / "refs/reference_design.csv")
assert hashlib.sha256((FB / "refs/reference_design.csv").read_bytes()).hexdigest() == (FB / "refs/reference_design.sha256").read_text().split()[0]
out = FB / "refs" / v; out.mkdir(parents=True, exist_ok=False)
o = BiphasicOracle(v, f"FBREF_{v}", out / "trials.jsonl", (), par); obs = o.evaluate_many(des[["z1", "z2", "z3"]].to_numpy())
print(f"done {v}: {sum(x.success for x in obs)}/{len(obs)} accepted", flush=True)
