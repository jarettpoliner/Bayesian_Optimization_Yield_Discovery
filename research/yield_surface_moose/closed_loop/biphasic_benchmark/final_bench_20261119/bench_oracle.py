"""BO-side oracle adapter for the final B1/B2 benchmark (BLIND: never imports the simulator-side oracle or the material table).
Each trial is one subprocess of oracle/biphasic_worker_r2.py in the moose environment (corners pin, h40, dr_rel 0.05 = the validated
configuration); only the public observation comes back.  Exact (variant, direction) reuse from the listed trial logs of THIS benchmark.
Every request -- new, reused or failed -- is appended to the campaign's trials.jsonl."""
import json, subprocess, threading, time, sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
FB = Path(__file__).resolve().parent; BB = FB.parent; CL = BB.parent
for p in (CL, CL / "v2"):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))
from run_bo_moose import to_yield_observation, _key                                  # noqa: E402
from yield_surface_discovery import YieldOracle                                     # noqa: E402
MARK = "@@OBSERVATION@@"
MM = [str(Path.home() / ".micromamba/bin/micromamba"), "run", "-r", str(Path.home() / "micromamba"), "-n", "moose", "python"]
WORKER = BB / "oracle/biphasic_worker_r2.py"
CONFIG = dict(pin_mode="corners", nel=40, dr_rel=0.05)


class BiphasicOracle(YieldOracle):
    def __init__(self, variant, run_prefix, log_path, reuse_paths=(), max_parallel=8):
        self.variant, self.prefix, self.log, self.maxp = variant, run_prefix, Path(log_path), max_parallel
        self.lock, self.records, self.n = threading.Lock(), [], 0; self.table = {}
        for p in reuse_paths:
            for l in open(p):
                r = json.loads(l)
                if r["variant"] == variant and r["config"] == CONFIG and "reused_from" not in r:
                    self.table.setdefault(_key(r["direction"]), (str(p), r))

    def _write(self, rec):
        with self.lock:
            self.records.append(rec)
            with open(self.log, "a") as fh:
                fh.write(json.dumps(rec, default=float) + "\n")

    def _run(self, d, call, batch):
        k = _key(d)
        if k in self.table:
            src, r = self.table[k]; rec = dict(r, call_no=call, batch=batch, reused_from=src, original_tag=r["tag"]); self._write(rec)
            return to_yield_observation(d, r["observation"])
        tag = f"{self.prefix}_t{call:02d}"; t0 = time.time()
        p = subprocess.run(MM + [str(WORKER), "--variant", self.variant, "--direction", *[repr(float(x)) for x in d], "--tag", tag, "--pin-mode", CONFIG["pin_mode"],
                                 "--nel", str(CONFIG["nel"]), "--dr-rel", repr(CONFIG["dr_rel"])], capture_output=True, text=True, cwd=WORKER.parent)
        line = next((l for l in p.stdout.splitlines()[::-1] if l.startswith(MARK)), None)
        obs = json.loads(line[len(MARK):]) if line else dict(trial_id=f"{self.variant}/{tag}", yield_radius=None, accepted=False, converged=False,
                                                              error=f"worker exit {p.returncode}: {p.stderr[-2000:]}")
        rec = dict(call_no=call, batch=batch, tag=tag, variant=self.variant, config=CONFIG, direction=[float(x) for x in d], started_unix=t0,
                   worker_wall_s=time.time() - t0, worker_returncode=p.returncode, observation=obs)
        self._write(rec); print(f"  [{time.strftime('%H:%M:%S')}] {tag} r={obs.get('yield_radius')} accepted={obs.get('accepted')}", flush=True)
        return to_yield_observation(d, obs)

    def _next(self):
        with self.lock:
            self.n += 1; return self.n

    def evaluate(self, direction):
        return self._run(np.asarray(direction, float), self._next(), "sequential")

    def evaluate_many(self, directions):
        D = np.atleast_2d(directions); nos = [self._next() for _ in D]
        with ThreadPoolExecutor(min(self.maxp, len(D))) as ex:
            return list(ex.map(lambda a: self._run(a[0], a[1], "batch"), zip(D, nos)))
