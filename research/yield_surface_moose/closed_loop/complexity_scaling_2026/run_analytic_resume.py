#!/usr/bin/env python3
"""RESUMABLE continuation of the analytic BO campaigns (complexity-scaling study), written after the 2026-11-19 17:09 reboot.  No MOOSE.

Why this is exact: with schedule A20 (every = 1) the frozen V2Controller re-fits ML-II from scratch at every step (fixed restarts,
random_state 0), so the selection at step n depends only on the n recorded (axis, radius) pairs; the analytic oracle is deterministic.
ResumableController.run_from(prefix) therefore executes the SAME adaptive loop as bomes_v2.V2Controller.run (copied line for line; the
frozen class is subclassed, not modified) starting from the recorded prefix instead of from the initial design.  Verified by
--selftest: (a) run_from(prefix of a fresh frozen run) reproduces the fresh frozen run bitwise; (b) every recorded prefix equals a fresh
frozen run over its first 40 tests.

Checkpointing: every new observation is appended (flushed) to evaluations_cont.jsonl and its step row to steps_cont.jsonl, so a further
interruption is resumed by re-running this script (prefix = original log + continuation log).
    python run_analytic_resume.py <level> <start> [--budget 1000]          ->  campaigns_r2/<level>/BO_s<start>/
    python run_analytic_resume.py --selftest"""
import os, sys
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "1")
import argparse, hashlib, json, time
from pathlib import Path
import numpy as np
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import cs_common as C                                                              # noqa: E402
import bomes_v2 as B                                                               # noqa: E402  frozen
from yield_surface_discovery.oracle import AnalyticalOracle                        # noqa: E402  frozen V1


class ResumableController(B.V2Controller):
    def run_from(self, start, prefix_axes, prefix_radii, on_step=None):
        """identical to V2Controller.run (no failures occur with the analytic oracle) but starting from a recorded prefix."""
        s = self.s; seq, radii = list(prefix_axes), list(prefix_radii); steps, scores, unavailable = [], [], set()
        if not seq:                                                    # initial design exactly as V2Controller.run (analytic oracle never fails)
            seq = list(self.initial_design(start)); radii = [o.yield_radius for o in self.oracle.evaluate_many(self.C[seq])]
            if on_step:
                for c, r in zip(seq, radii):
                    on_step(c, r, None)
        while len(seq) < s.budget:
            gp = self.fit(self.C[seq], radii)
            mean, var = gp.predict_log(self.C); sc = var.copy(); sc[seq] = -np.inf; sc[list(unavailable)] = -np.inf
            j, tied = B.select(sc, "deterministic", s.tie_rel, s.tie_abs); hp = gp.hyperparameters()
            step = dict(test_number=len(seq) + 1, candidate=int(j), latent_var=float(var[j]), entropy=float(B.entropy(var[j])), gp_mean_log_radius=float(mean[j]),
                        tie_set_size=int(len(tied)), **{k: v for k, v in hp.items() if k != "at_bounds"}, at_bounds=";".join(hp["at_bounds"]))
            o = self.oracle.evaluate(self.C[j]); seq.append(int(j)); radii.append(o.yield_radius); steps.append(step); scores.append(sc)
            if on_step:
                on_step(int(j), o.yield_radius, step)
        return seq, radii, steps, scores


def axis_index(lib, d):
    i = int(np.argmin(np.abs(lib.directions - np.asarray(d)).sum(1))); assert np.array_equal(lib.directions[i], np.asarray(d, float)); return i


def read_log(p):
    return [json.loads(l) for l in open(p)] if p.exists() else []


def selftest():
    lib = C.library(); out = {}
    for level, start in (("A5_rand10_p32", 1), ("A2_hill", 3)):
        fn = C.radius_fn(level); st = B.V2Settings(schedule="A20", symmetry="antipodal", budget=45)
        fresh = B.V2Controller(lib, AnalyticalOracle(fn, "t"), st).run(start); fseq = fresh["sequence"]
        rc = ResumableController(lib, AnalyticalOracle(fn, "t"), st); rr = [o.yield_radius for o in fresh["observations"]]
        seq, rad, steps, sc = rc.run_from(start, fseq[:25], rr[:25])
        out[f"{level}_s{start}_resume_equals_fresh"] = seq == fseq and np.array_equal(np.array(sc), np.array(fresh["scores"][17:])) and rad == rr
        rc0 = ResumableController(lib, AnalyticalOracle(fn, "t"), st); seq0, *_ = rc0.run_from(start, [], [])
        out[f"{level}_s{start}_from_scratch_equals_fresh"] = seq0 == fseq
    st = B.V2Settings(schedule="A20", symmetry="antipodal", budget=40)
    for d in sorted((HERE / "campaigns").glob("A*/BO_s*")):
        level, start = d.parent.name, int(d.name.split("_s")[1]); rec = read_log(d / "evaluations.jsonl")
        fresh = B.V2Controller(lib, AnalyticalOracle(C.radius_fn(level), "t"), st).run(start)
        out[f"log_{level}_s{start}_first40"] = [axis_index(lib, r["direction"]) for r in rec[:40]] == fresh["sequence"] and \
            all(abs(r["yield_radius"] - o.yield_radius) == 0 for r, o in zip(rec[:40], fresh["observations"]))
    out["PASS"] = all(out.values()); (HERE / "logs").mkdir(exist_ok=True); json.dump(out, open(HERE / "logs/resume_selftest.json", "w"), indent=1); print(json.dumps(out, indent=1))


def main(level, start, budget):
    lib = C.library(); fn = C.radius_fn(level); out = HERE / "campaigns_r2" / level / f"BO_s{start}"; out.mkdir(parents=True, exist_ok=True)
    orig = HERE / "campaigns" / level / f"BO_s{start}" / "evaluations.jsonl"; recs = read_log(orig) + read_log(out / "evaluations_cont.jsonl")
    if orig.exists():
        (out / "prefix_source.json").write_text(json.dumps(dict(path=str(orig), n=len(read_log(orig)), sha256=hashlib.sha256(orig.read_bytes()).hexdigest()), indent=1))
    seq = [axis_index(lib, r["direction"]) for r in recs]; radii = [float(r["yield_radius"]) for r in recs]
    assert all(r["yield_radius"] == float(fn(np.asarray(r["direction"])[None])[0]) for r in recs[-5:]), "recorded radii must equal the analytic oracle"
    ev, sl = open(out / "evaluations_cont.jsonl", "a"), open(out / "steps_cont.jsonl", "a"); n0 = len(seq); t0 = time.time()

    def on_step(c, r, step):
        ev.write(json.dumps(dict(call=None, source="analytic_synthetic_oracle", direction=[float(x) for x in lib.directions[c]], yield_radius=r, success=True)) + "\n"); ev.flush()
        if step:
            sl.write(json.dumps(dict(step, wall_unix=time.time())) + "\n"); sl.flush()
    st = B.V2Settings(schedule="A20", symmetry="antipodal", budget=budget)
    seq, radii, steps, _ = ResumableController(lib, AnalyticalOracle(fn, f"analytic:{level}"), st).run_from(start, seq, radii, on_step)
    json.dump(dict(axes=seq, radii=radii), open(out / "design.json", "w"))
    json.dump(dict(level=level, method="BO", start=start, budget=budget, resumed_from_n=n0, wall_s_this_session=time.time() - t0, library_hash=C.library_hash(lib),
                   frozen_sources=C.file_hashes([C.CL / "v2/bomes_v2.py", *sorted((C.REL / "yield_surface_discovery").glob("*.py"))]),
                   resume_code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()), open(out / "run_info.json", "w"), indent=1, default=str)
    print(f"done {level} s{start}: {len(seq)} obs (resumed from {n0})", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("level", nargs="?"); ap.add_argument("start", nargs="?", type=int); ap.add_argument("--budget", type=int, default=1000)
    ap.add_argument("--selftest", action="store_true"); a = ap.parse_args()
    selftest() if a.selftest else main(a.level, a.start, a.budget)
