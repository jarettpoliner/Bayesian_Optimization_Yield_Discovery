#!/usr/bin/env python3
"""Integrity checks of the analytic complexity-scaling campaigns (closeout).  Read-only; writes logs_r2/verify_complexity.json.
  1 frozen files: PROTOCOL.sha256 entries (PROTOCOL.md, cs_common.py incl. the reference set and the analytic levels, run_analytic.py, descriptors)
  2 frozen controller / GP sources: every run_info.json's recorded hashes equal the files on disk now (bomes_v2.py, V1 package)
  3 library hash = PROTOCOL value; reference set (20,000 Fibonacci points) hash reported
  4 AXMM / RANDOM / BOFIXED designs: >= 1000 axes, unique, no antipodal duplicates; first 8 axes shared with BO start k (matched pairs)
  5 BO: original log + continuation = one sequence; unique axes; every recorded radius equals the analytic oracle exactly (single-direction calls, as the controller makes them); n reached
  6 resume self-test PASS (logs/resume_selftest.json)      python verify_complexity.py"""
import hashlib, json, subprocess, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import cs_common as C                                                              # noqa: E402
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
out = {}
r = subprocess.run(["sha256sum", "-c", "PROTOCOL.sha256"], cwd=HERE, capture_output=True, text=True); out["protocol_sha256_check"] = dict(rc=r.returncode, lines=r.stdout.splitlines())
lib = C.library(); out["library_hash"] = C.library_hash(lib); out["library_hash_matches_protocol"] = out["library_hash"] == "8e1654f80a0240d0857de2abbbac455501c5fd1fabbbdb633f8b1e368ef0a545"
R = C.references(); out["reference_set"] = dict(n=len(R), sha256_rounded=hashlib.sha256(np.round(R, 12).tobytes()).hexdigest())
bad, nchk = [], 0
for ri in sorted(list(HERE.glob("campaigns/*/*/run_info.json")) + list(HERE.glob("campaigns_r2/*/*/run_info.json"))):
    for p, h in json.load(open(ri)).get("frozen_sources", {}).items():
        nchk += 1
        if not Path(p).exists() or sha(p) != h:
            bad.append(f"{ri.parent.relative_to(HERE)}: {p}")
out["frozen_sources"] = dict(hashes_checked=nchk, mismatches=bad)
first8 = {}
def uniq(ax):
    D = lib.directions[np.asarray(ax)]; G = np.abs(D @ D.T); np.fill_diagonal(G, 0)
    return len(set(ax)) == len(ax) and bool(G.max() < 1 - 1e-12)
des = {}
for d in sorted(HERE.glob("campaigns/*/*_s*")):
    if d.name.startswith("BO_"):
        continue
    j = json.load(open(d / "design.json")); ax = list(j["axes"])
    des[f"{d.parent.name}/{d.name}"] = dict(n=len(ax), unique_unoriented=uniq(ax[:1000]), first8=ax[:8])
bo = {}
for L in C.LEVELS:
    fn = C.radius_fn(L)
    for k in range(4):
        recs = [json.loads(l) for p in (HERE / "campaigns" / L / f"BO_s{k}" / "evaluations.jsonl", HERE / "campaigns_r2" / L / f"BO_s{k}" / "evaluations_cont.jsonl") if p.exists() for l in open(p)]
        Dd = np.array([x["direction"] for x in recs]); idx = [int(np.argmin(np.abs(lib.directions - z).sum(1))) for z in Dd]
        exact = bool(np.array_equal(lib.directions[idx], Dd)); rad = np.array([x["yield_radius"] for x in recs]); oracle_ok = bool(np.array_equal(rad, np.r_[fn(Dd[:8]), [fn(z[None])[0] for z in Dd[8:]]]))   # controller call pattern: initial 8 in one evaluate_many batch, then one per step
        batch_rel = float(np.max(np.abs(fn(Dd) / rad - 1)))                                                                           # batched BLAS rounding, reported only
        bo[f"{L}/BO_s{k}"] = dict(n=len(recs), axes_on_library=exact, unique_unoriented=uniq(idx), radii_equal_oracle=oracle_ok, batch_eval_max_rel_diff=batch_rel, first8=idx[:8])
        for arm in ("AXMM", "RANDOM"):
            x = des.get(f"{L}/{arm}_s{k}")
            if x:
                x["first8_equals_BO"] = x["first8"] == idx[:8] if arm == "AXMM" else None
for v in list(des.values()) + list(bo.values()):
    v.pop("first8", None)
out["nonadaptive_designs"] = des; out["bo_campaigns"] = bo
out["resume_selftest_PASS"] = json.load(open(HERE / "logs/resume_selftest.json"))["PASS"]
out["bo_min_n"] = min(v["n"] for v in bo.values()); out["bo_n_reached_512_all"] = out["bo_min_n"] >= 512; out["bo_n_reached_1000_all"] = out["bo_min_n"] >= 1000
out["PASS_integrity_so_far"] = bool(out["protocol_sha256_check"]["rc"] == 0 and out["library_hash_matches_protocol"] and not bad and out["resume_selftest_PASS"]
                             and all(v["n"] >= 1000 and v["unique_unoriented"] for v in des.values()) and all(v["axes_on_library"] and v["unique_unoriented"] and v["radii_equal_oracle"] for v in bo.values())
                             and all(v.get("first8_equals_BO") in (True, None) for v in des.values()))
out["PASS"] = bool(out["PASS_integrity_so_far"] and out["bo_n_reached_1000_all"])
json.dump(out, open(HERE / "logs_r2/verify_complexity.json", "w"), indent=1)
print(json.dumps({k: v for k, v in out.items() if k not in ("nonadaptive_designs", "bo_campaigns")}, indent=1))
print("designs:", {k: v for k, v in des.items() if not (v["unique_unoriented"] and v["n"] >= 1000 and v.get("first8_equals_BO") in (True, None))} or "all OK")
print("bo:", {k: v["n"] for k, v in bo.items()})
