"""Batch execution of MOOSE truss observations with (i) a guard against other MOOSE workloads, (ii) a CPU-time ledger and (iii) a HARD stage cap.

* Guard: before every batch, any combined-opt / *-opt MOOSE process not started by this runner stops the batch (feasibility-stage rule:
  no other significant MOOSE workload during the pilot).
* Ledger: logs/cpu_ledger.jsonl, one line per observation (child CPU user+sys from getrusage).  Sum of the ledger = stage CPU use.
* Cap: a batch refuses to start, and stops submitting, once ledger total + (in-flight estimate) would exceed the cap (default 2 CPU-h).
* Nothing is ever deleted; every observation gets its own directory.
"""
import json, os, subprocess, time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "logs" / "cpu_ledger.jsonl"
STAGE_CAP_CPU_S = 2 * 3600.0


def ledger_total():
    tot = 0.0
    if LEDGER.exists():
        for l in LEDGER.read_text().splitlines():
            x = json.loads(l); tot += x.get("cpu_s", x.get("cpu_s_estimate_upper", 0.0))
    return tot


def ledger_add(item, cpu_s, **kw):
    with open(LEDGER, "a") as f:
        f.write(json.dumps(dict(time=time.strftime("%Y-%m-%dT%H:%M:%S"), item=item, cpu_s=cpu_s, **kw)) + "\n")


def foreign_moose():
    """MOOSE processes (any *-opt app) not descending from this process."""
    out = subprocess.run(["ps", "-eo", "pid,ppid,args"], capture_output=True, text=True).stdout.splitlines()[1:]
    mine = {os.getpid()}; rows = [l.split(None, 2) for l in out if l.strip()]
    changed = True
    while changed:
        changed = False
        for pid, ppid, _ in rows:
            if int(ppid) in mine and int(pid) not in mine:
                mine.add(int(pid)); changed = True
    return [a for pid, _, a in rows if int(pid) not in mine and ("-opt " in a + " ") and ("combined-opt" in a or "solid_mechanics-opt" in a or "moose" in a.split()[0])]


def _one(args):
    import truss_oracle as T
    d, wd, cal, opt = args
    try:
        rec = T.observe(np.asarray(d, float), wd, cal, **opt)
    except Exception as e:                                       # recorded as a failed observation, never silently dropped
        rec = dict(direction=list(np.asarray(d, float) / np.linalg.norm(d)), accepted=False, reason=f"exception: {e!r}", cpu_s=0.0, wall_s=0.0)
    rec["workdir"] = str(wd); (Path(wd) / "record.json").write_text(json.dumps(rec, indent=1, default=float)); return rec


def run_batch(jobs, cal, workers=4, cap_cpu_s=STAGE_CAP_CPU_S, label="batch", opt=None):
    """jobs: list of (direction, workdir).  returns list of records (same order); None for jobs not run because of the cap."""
    f = foreign_moose()
    if f:
        raise SystemExit(f"GUARD: other MOOSE workload running, refusing to start {label}: {f[:3]}")
    if ledger_total() >= cap_cpu_s:
        raise SystemExit(f"CAP: ledger {ledger_total():.0f} s >= cap {cap_cpu_s:.0f} s")
    for _, wd in jobs:
        if Path(wd).exists():
            raise SystemExit(f"refusing to reuse existing directory {wd}")
    out = [None] * len(jobs); t0 = time.time(); est = 10.0
    with ProcessPoolExecutor(workers) as ex:
        futs = {}; it = iter(enumerate(jobs))
        def submit():
            for k, (d, wd) in it:
                if ledger_total() + est * (len(futs) + 1) > cap_cpu_s:
                    return False
                futs[ex.submit(_one, (d, wd, cal, opt or {}))] = k
                if len(futs) >= workers:
                    return True
            return True
        submit()
        while futs:
            for fu in as_completed(list(futs)):
                k = futs.pop(fu); r = fu.result(); out[k] = r; ledger_add(f"{label}:{Path(r['workdir']).name}", r.get("cpu_s", 0.0), wall_s=r.get("wall_s"))
                est = max(est, 1.5 * r.get("cpu_s", 0.0))
                break
            if foreign_moose():
                print("GUARD: foreign MOOSE workload appeared; not submitting more"); it = iter(())
            submit()
    print(f"{label}: {sum(r is not None for r in out)}/{len(jobs)} run, wall {time.time() - t0:.0f} s, ledger total {ledger_total() / 3600:.3f} CPU-h")
    return out
