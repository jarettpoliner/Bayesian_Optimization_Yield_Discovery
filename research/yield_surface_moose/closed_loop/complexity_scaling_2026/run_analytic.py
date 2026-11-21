#!/usr/bin/env python3
"""One ANALYTIC (synthetic-oracle) campaign of the complexity-scaling study.  No MOOSE.  Every observation is an exact analytic evaluation.

    python run_analytic.py <level> BO <start>        frozen V2Controller (V2Settings schedule A20, antipodal, budget N) on the 4,096-axis library
    python run_analytic.py <level> AXMM <start>      frozen bomes_v2.axial_maximin, same first axis as BO start <start>
    python run_analytic.py <level> RANDOM <start>    the same 8 axial-maximin initial axes as BO start <start>, then default_rng(9000+start) unique axes
    python run_analytic.py ANY BOFIXED <start>       DIAGNOSTIC arm: frozen V2Controller with schedule FIXED (response-independent; same for all levels)

Frozen pieces used unmodified: V1 AnalyticalOracle, CandidateLibrary, acquisition.select; bomes_v2.V2Controller / GPv2 / axial_maximin / entropy.
RecordingController only wraps V2Controller.fit (calls super().fit unchanged and returns the same GP) to store the posterior mean on the
candidates and the hyper-parameters after every fit.  It never alters a decision (checked in the self-test: identical to plain V2Controller).
Outputs: campaigns/<level>/<method>_s<start>/{design.json, steps.csv, scores_f32.npy, mean_f32.npy, fits.csv, evaluations.jsonl, run_info.json}"""
import os, sys
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "1")
import json, time, argparse
from pathlib import Path
import numpy as np, pandas as pd
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import cs_common as C                                                              # noqa: E402
import bomes_v2 as B                                                               # noqa: E402  frozen
from yield_surface_discovery.oracle import AnalyticalOracle                        # noqa: E402  frozen V1


class LoggingOracle(AnalyticalOracle):
    def __init__(self, fn, name, log):
        super().__init__(fn, name); self.log = open(log, "w"); self.n = 0

    def evaluate_many(self, directions):
        out = super().evaluate_many(directions)
        for o in out:
            self.n += 1; self.log.write(json.dumps(dict(call=self.n, source="analytic_synthetic_oracle", direction=[float(x) for x in o.requested_direction],
                                                        yield_radius=o.yield_radius, success=o.success)) + "\n")
        self.log.flush(); return out


class RecordingController(B.V2Controller):
    rec_mean, rec_fit = [], []

    def fit(self, obs_dirs, radii, schedule=None):
        t = time.perf_counter(); gp = super().fit(obs_dirs, radii, schedule); dt = time.perf_counter() - t
        m, _ = gp.predict_log(self.C); hp = gp.hyperparameters(); self.rec_mean.append(m.astype(np.float32))
        self.rec_fit.append(dict(n_real=hp["n_real"], mode=hp["mode"], amplitude=hp["amplitude"], length_scale=hp["length_scale"], white_noise=hp["white_noise"],
                                 at_bounds=";".join(hp["at_bounds"]), lml=hp["log_marginal_likelihood"], fit_s=dt)); return gp


def main(level, method, start, budget, outroot):
    lib = C.library(); fn = C.radius_fn("A1_J2" if level == "ANY" else level)
    out = outroot / level / f"{method}_s{start}"; out.mkdir(parents=True, exist_ok=False)
    oracle = LoggingOracle(fn, f"analytic:{level}", out / "evaluations.jsonl"); t0 = time.time()
    st = B.V2Settings(schedule="FIXED" if method == "BOFIXED" else "A20", symmetry="antipodal", budget=budget)
    ctrl = RecordingController(lib, oracle, st); RecordingController.rec_mean, RecordingController.rec_fit = [], []
    assert np.array_equal(ctrl.C, lib.directions) and len(ctrl.C) == 4096, "axis library must be the 4,096 Fibonacci axes (no antipodes)"
    info = dict(level=level, method=method, start=start, budget=budget, library_hash=C.library_hash(lib), settings=st.__dict__,
                frozen_sources=C.file_hashes([C.CL / "v2/bomes_v2.py", *sorted((C.REL / "yield_surface_discovery").glob("*.py"))]))
    if method in ("BO", "BOFIXED"):
        res = ctrl.run(start); seq = res["sequence"]; radii = [o.yield_radius for o in res["observations"]]
        pd.DataFrame([{k: v for k, v in s.items() if k != "tie_set"} for s in res["steps"]]).to_csv(out / "steps.csv", index=False)
        np.save(out / "scores_f32.npy", np.array([np.where(np.isfinite(s), s, np.nan) for s in res["scores"]], np.float32))
        np.save(out / "mean_f32.npy", np.array(RecordingController.rec_mean, np.float32)); pd.DataFrame(RecordingController.rec_fit).to_csv(out / "fits.csv", index=False)
    else:
        f = ctrl.first(start); init = ctrl.initial_design(start)
        if method == "AXMM":
            seq = [int(c) for c in B.axial_maximin(ctrl.C, f, budget)]
            assert seq[:8] == init, "AXMM must share the BO initial design"
        else:
            rest = [c for c in np.random.default_rng(9000 + start).permutation(len(ctrl.C)) if c not in set(init)]
            seq = list(init) + [int(c) for c in rest[: budget - 8]]
        assert len(set(seq)) == budget
        radii = [o.yield_radius for o in oracle.evaluate_many(ctrl.C[seq])]
    assert len(seq) == budget and len(set(seq)) == budget
    json.dump(dict(axes=[int(c) for c in seq], radii=[float(r) for r in radii]), open(out / "design.json", "w"))
    info["wall_s"] = time.time() - t0; info["n_oracle_calls"] = oracle.n; json.dump(info, open(out / "run_info.json", "w"), indent=1, default=str)
    print(f"done {level} {method} s{start}: {budget} obs, {info['wall_s']:.0f} s", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("level"); ap.add_argument("method", choices=["BO", "AXMM", "RANDOM", "BOFIXED"]); ap.add_argument("start", type=int)
    ap.add_argument("--budget", type=int, default=1000); ap.add_argument("--out", default=str(HERE / "campaigns")); a = ap.parse_args()
    main(a.level, a.method, a.start, a.budget, Path(a.out))
