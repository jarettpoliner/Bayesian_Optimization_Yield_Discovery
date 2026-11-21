#!/usr/bin/env python3
"""Checkpoint monitor: when a resumed BO campaign has >= 512 recorded observations (original log + continuation log), send SIGTERM to its
process right after the observation is written (the next write is >= several seconds away, so no partial line).  Resume later with
    python run_analytic_resume.py <level> <start> --budget 1000
Writes logs_r2/stop_at_512_status.json."""
import json, os, signal, subprocess, time
from pathlib import Path
HERE = Path(__file__).resolve().parent; status = {}
count = lambda p: sum(1 for _ in open(p)) if p.exists() else 0
while True:
    ps = subprocess.run(["ps", "-eo", "pid,args"], capture_output=True, text=True).stdout.splitlines(); alive = 0
    for line in ps:
        if "run_analytic_resume.py" in line and "--budget" in line and "python" in line:
            pid, *args = line.split(); i = args.index("run_analytic_resume.py") if "run_analytic_resume.py" in args else [a.endswith("run_analytic_resume.py") for a in args].index(True)
            level, start = args[i + 1], args[i + 2]; alive += 1
            n = count(HERE / "campaigns" / level / f"BO_s{start}" / "evaluations.jsonl") + count(HERE / "campaigns_r2" / level / f"BO_s{start}" / "evaluations_cont.jsonl")
            if n >= 512 and f"{level}_s{start}" not in status:
                os.kill(int(pid), signal.SIGTERM); status[f"{level}_s{start}"] = dict(n_at_stop=n, time=time.strftime("%H:%M:%S"))
                (HERE / "logs_r2/stop_at_512_status.json").write_text(json.dumps(status, indent=1))
    if alive == 0:
        break
    time.sleep(2)
