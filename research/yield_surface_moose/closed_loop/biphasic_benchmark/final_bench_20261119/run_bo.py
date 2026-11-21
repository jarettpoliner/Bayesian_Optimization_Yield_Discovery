#!/usr/bin/env python3
"""Stage 3/4: one V2-A20 BO+MES campaign (FROZEN v2/bomes_v2.V2Controller, V2Settings(schedule='A20', symmetry='antipodal', budget=33)).
Antipodal symmetry ENABLED (validated in S1b-min): 207 unoriented axes, pseudo-observations not counted as experiments.  Shear mirror off.
Exact (variant, direction) reuse only from this benchmark's pooled non-adaptive trials of the SAME variant.  BLIND (no material / reference access).
    python run_bo.py B1 0            ->  campaigns/B1/bo_s0/
    python run_bo.py B1 0 --replay   ->  campaigns/B1/bo_s0_replay/ (ReplayOracle on bo_s0/trials.jsonl; no MOOSE)"""
import argparse, hashlib, json, sys
from pathlib import Path
FB = Path(__file__).resolve().parent; CL = FB.parents[1]; sys.path.insert(0, str(FB)); sys.path.insert(0, str(CL / "v2")); sys.path.insert(0, str(CL))
import designs_common as D                                                          # noqa: E402
from bench_oracle import BiphasicOracle                                             # noqa: E402
from bomes_v2 import V2Controller, V2Settings                                       # noqa: E402
import run_v2_moose as RV                                                           # noqa: E402
from run_bo_moose import ReplayOracle, bo_source_hashes                             # noqa: E402
ap = argparse.ArgumentParser(); ap.add_argument("variant"); ap.add_argument("start", type=int); ap.add_argument("--replay", action="store_true"); a = ap.parse_args()
assert all(v["matches_release"] for v in bo_source_hashes().values())
lib, AX = D.library(); base = FB / "campaigns" / a.variant / f"bo_s{a.start}"
if a.replay:
    out = Path(str(base) + "_replay"); out.mkdir(parents=True, exist_ok=False); oracle = ReplayOracle(base / "trials.jsonl")
else:
    out = base; out.mkdir(parents=True, exist_ok=False)
    oracle = BiphasicOracle(a.variant, f"FB_{a.variant}_bo_s{a.start}", out / "trials.jsonl", [FB / "campaigns" / a.variant / "nonadaptive/trials.jsonl"], 8)
res = V2Controller(lib, oracle, V2Settings(schedule="A20", symmetry="antipodal", budget=D.BUDGET)).run(a.start)
ns = __import__("argparse").Namespace(replay=a.replay, specimen=a.variant, mode="v2", schedule="A20", start=a.start, budget=D.BUDGET, out=str(out))
RV.write_outputs(out, lib, res, "v2", ns, oracle)
info = json.load(open(out / "run_info.json")); sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
info["final_bench_sources"] = {f: sha(FB / f) for f in ("bench_oracle.py", "designs_common.py", "run_bo.py")}; info["worker_r2"] = sha(FB.parent / "oracle/biphasic_worker_r2.py")
(out / "run_info.json").write_text(json.dumps(info, indent=1, default=str)); print(f"done {a.variant} s{a.start}: {len(res['sequence'])} tests", flush=True)
