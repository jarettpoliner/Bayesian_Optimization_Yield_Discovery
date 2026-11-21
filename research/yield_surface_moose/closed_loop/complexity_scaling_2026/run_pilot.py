#!/usr/bin/env python3
"""Physical pilot (PROTOCOL.md): 24 area-uniform axes per configuration P2-P6 (+ 2 h=1/80 axes for P6).  NEW MOOSE experiments, logged to
pilot/<config>/trials.jsonl with source = "moose_new".  Resumable: a (config, nel, direction) already in the log with a record is not re-run.
Exact reuse only for identical (config, nel, direction)."""
import json, subprocess, sys, time, threading
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import cs_common as C                                                              # noqa: E402
MM = [str(Path.home() / ".micromamba/bin/micromamba"), "run", "-r", str(Path.home() / "micromamba"), "-n", "moose", "python"]
CONFIGS = ["P2_single", "P3_regular", "P4_clustered_s4401", "P5_multiscale_s5501", "P6_aligned_hiAR_s6601"]
F = C.fibonacci_sphere(48); AX = F[F[:, 2] > 0]; assert len(AX) == 24
lock = threading.Lock()


def job(cfg, i, d, nel):
    log = HERE / "pilot" / cfg / "trials.jsonl"; tag = f"pilot_h{nel}_a{i:02d}"
    with lock:
        done = {json.loads(l)["tag"] for l in open(log)} if log.exists() else set()
    if tag in done:
        return
    t0 = time.time()
    p = subprocess.run(MM + [str(HERE / "pilot_worker.py"), "--config", cfg, "--direction", *[repr(float(x)) for x in d], "--tag", tag, "--nel", str(nel)],
                       capture_output=True, text=True, cwd=HERE)
    line = next((l for l in p.stdout.splitlines()[::-1] if l.startswith("@@OBSERVATION@@")), None)
    obs = json.loads(line[15:]) if line else dict(yield_radius=None, accepted=False, error=f"exit {p.returncode}: {p.stderr[-1500:]}")
    rec = dict(config=cfg, tag=tag, nel=nel, axis_index=i, direction=[float(x) for x in d], source="moose_new", wall_s=time.time() - t0, observation=obs)
    with lock:
        log.parent.mkdir(parents=True, exist_ok=True)
        with open(log, "a") as fh:
            fh.write(json.dumps(rec, default=float) + "\n")
    print(time.strftime("%H:%M:%S"), cfg, tag, obs.get("yield_radius"), obs.get("accepted"), f"{rec['wall_s']:.0f}s", flush=True)


if __name__ == "__main__":
    par = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    jobs = [(c, i, d, 40) for i, d in enumerate(AX) for c in CONFIGS] + [("P6_aligned_hiAR_s6601", i, AX[i], 80) for i in (0, 12)]
    with ThreadPoolExecutor(par) as ex:
        list(ex.map(lambda a: job(*a), jobs))
    print("PILOT_DONE", flush=True)
