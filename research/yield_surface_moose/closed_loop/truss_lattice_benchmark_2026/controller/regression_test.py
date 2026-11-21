"""TV2C regression test: the wrapper must reproduce the ORIGINAL frozen controller decisions bitwise on a known recorded campaign.
Test: B1 BO start 0 of the final B1/B2 benchmark (release data), replayed (ReplayOracle, recorded MOOSE observations, no MOOSE run):
  R1  wrapper path with the original 207-axis library object                 -> sequence, step records and all 26x207 acquisition scores identical
  R2  original library written to CSV in the TLIB format and re-loaded through the wrapper's loader -> identical directions (bitwise) and R1 again
  R3  wrapper baselines reproduce the recorded axial-maximin / random designs of B1 (nonadaptive/designs.json) on the original library."""
import json, sys, tempfile
from pathlib import Path
import numpy as np, pandas as pd
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT / "controller"))
import truss_v2_compat as W
CL = ROOT.parent; REL = CL.parent / "release_B1_B2_BO_benchmark_20261119"
sys.path.insert(0, str(CL / "biphasic_benchmark/final_bench_20261119"))   # original (non-release) copy; resolves configs/ correctly
import designs_common as DC
from run_bo_moose import ReplayOracle

hashes = W.check_frozen(); camp = REL / "data/campaigns/B1/bo_s0"
tests = pd.read_csv(camp / "tests.csv", float_precision="round_trip"); steps = pd.read_csv(camp / "steps.csv", float_precision="round_trip"); scores = np.load(camp / "acquisition_scores.npy")
lib, AX = DC.library(); out = {}


def replay(lib_):
    res = W.V2Controller(lib_, ReplayOracle(camp / "trials.jsonl"), W.V2Settings(schedule="A20", symmetry="antipodal", budget=33)).run(0)
    seq_ok = res["v1_indices"] == tests.v1_index.tolist()
    st = pd.DataFrame(res["steps"]); cols = ["candidate", "v1_index", "latent_var", "entropy", "gp_mean_log_radius", "amplitude", "length_scale", "white_noise", "log_marginal_likelihood"]
    steps_ok = all(np.array_equal(st[c].to_numpy(), steps[c].to_numpy()) for c in cols)
    sc_ok = np.array_equal(np.array(res["scores"]), scores)
    r_ok = np.array_equal(np.array([o.yield_radius for o in res["observations"]]), tests.yield_radius.to_numpy())
    return dict(sequence=seq_ok, steps_bitwise=steps_ok, scores_bitwise=sc_ok, radii=r_ok, n_tests=len(res["sequence"]))


out["R1"] = replay(lib)
with tempfile.TemporaryDirectory() as td:
    f = Path(td) / "lib207.csv"
    pd.DataFrame(dict(id=lib.ids, label=lib.labels, family=lib.families, z1=lib.directions[:, 0], z2=lib.directions[:, 1], z3=lib.directions[:, 2])).to_csv(f, index=False, float_format="%.17g")
    lib2 = W.load_library(f)
    out["R2"] = dict(directions_bitwise=bool(np.array_equal(lib2.directions, lib.directions)), **replay(lib2))
fd = REL / "data/campaigns/B1/nonadaptive_designs.json"; d = json.load(open(fd))["designs"] if fd.exists() else None
if d is None:
    out["R3"] = "skipped: designs.json not in release data"
else:
    ok = {}
    for k in DC.STARTS:
        x = d[f"axmm_s{k}"]; bad = set(x.get("failed") or []); ok[f"axmm_s{k}"] = [c for c in W.maximin_design(lib, k, 45, "antipodal") if c not in bad][:len(x["axes"])] == x["axes"]
    for k in DC.RANDOM_SEEDS:
        x = d[f"random_{k}"]; bad = set(x.get("failed") or []); ok[f"random_{k}"] = [c for c in W.random_design(len(AX.directions), k, 45) if c not in bad][:len(x["axes"])] == x["axes"]
    out["R3"] = dict(all=all(ok.values()), detail=ok)
out["passed"] = bool(all(v for k in ("R1", "R2") for v in out[k].values() if isinstance(v, bool)) and (not isinstance(out["R3"], dict) or out["R3"]["all"]))
out["frozen_hashes"] = hashes; out["wrapper"] = W.VERSION
(ROOT / "controller/regression_result.json").write_text(json.dumps(out, indent=1)); print(json.dumps({k: v for k, v in out.items() if k != "frozen_hashes"}, indent=1))
