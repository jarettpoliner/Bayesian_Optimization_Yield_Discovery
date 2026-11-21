#!/usr/bin/env python3
"""Stages 2-4 pipeline with gates (Stage 5 evaluation is run separately afterwards).  Logs: logs/*.log; status: logs/pipeline_status.json"""
import hashlib, json, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
FB = Path(__file__).resolve().parent; PY = sys.executable; LOG = FB / "logs"; LOG.mkdir(exist_ok=True); st = {}
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
sys.path.insert(0, str(FB)); import designs_common as D                              # noqa: E402


def save():
    (LOG / "pipeline_status.json").write_text(json.dumps(st, indent=1, default=str))


def run(args, name):
    with open(LOG / f"{name}.log", "w") as fh:
        return subprocess.run([PY] + args, cwd=FB, stdout=fh, stderr=subprocess.STDOUT).returncode


def stop(msg):
    st["STOPPED"] = msg; save(); print("STOP:", msg, flush=True); sys.exit(2)


assert sha(FB / "PROTOCOL_FINAL.md") == (FB / "PROTOCOL_FINAL.sha256").read_text().split()[0]
st["protocol_sha256"] = sha(FB / "PROTOCOL_FINAL.md"); st["start"] = time.strftime("%H:%M:%S")
# Stage 2: references (both variants concurrently, 11 parallel each)
with ThreadPoolExecutor(2) as ex:
    rc = list(ex.map(lambda v: run(["run_references.py", v, "11"], f"refs_{v}"), D.VARIANTS))
for v in D.VARIANTS:
    recs = [json.loads(l) for l in open(FB / "refs" / v / "trials.jsonl")]; acc = sum(bool(r["observation"].get("accepted")) for r in recs)
    st[f"stage2_{v}"] = dict(returncode=rc[D.VARIANTS.index(v)], records=len(recs), accepted=acc, sha256=sha(FB / "refs" / v / "trials.jsonl"))
save()
if not all(st[f"stage2_{v}"]["accepted"] == 72 and st[f"stage2_{v}"]["records"] == 72 for v in D.VARIANTS):
    stop("stage 2 reference gate failed")
st["stage2_done"] = time.strftime("%H:%M:%S"); save()
for stage, v in ((3, "B1"), (4, "B2")):
    if run(["run_nonadaptive.py", v, "22"], f"nonadaptive_{v}") != 0:
        stop(f"stage {stage}: non-adaptive run failed")
    d = json.load(open(FB / "campaigns" / v / "nonadaptive/designs.json"))
    st[f"stage{stage}_{v}_nonadaptive"] = dict(n_new=d["n_new"], n_failed=d["n_failed"], all_33=all(len(x["axes"]) == 33 for x in d["designs"].values()), time=time.strftime("%H:%M:%S"))
    save()
    if not st[f"stage{stage}_{v}_nonadaptive"]["all_33"]:
        stop(f"stage {stage}: non-adaptive designs incomplete")
    with ThreadPoolExecutor(len(D.STARTS)) as ex:
        rcs = list(ex.map(lambda k: run(["run_bo.py", v, str(k)], f"bo_{v}_s{k}"), D.STARTS))
    ok = []
    for k, c in zip(D.STARTS, rcs):
        t = FB / "campaigns" / v / f"bo_s{k}" / "tests.csv"; ok.append(c == 0 and t.exists() and sum(1 for _ in open(t)) == 34)
    rep = [run(["run_bo.py", v, str(k), "--replay"], f"replay_{v}_s{k}") for k in D.STARTS]
    same = []
    for k in D.STARTS:
        a, b = FB / "campaigns" / v / f"bo_s{k}", FB / "campaigns" / v / f"bo_s{k}_replay"
        same.append((a / "tests.csv").read_bytes() == (b / "tests.csv").read_bytes() and (a / "steps.csv").read_bytes() == (b / "steps.csv").read_bytes()
                    and (a / "acquisition_scores.npy").read_bytes() == (b / "acquisition_scores.npy").read_bytes())
    st[f"stage{stage}_{v}_bo"] = dict(complete=ok, replay_rc=rep, replay_identical=same, time=time.strftime("%H:%M:%S")); save()
    if not (all(ok) and all(same)):
        stop(f"stage {stage}: BO campaigns incomplete or replay mismatch")
st["ALL_DONE"] = time.strftime("%H:%M:%S"); save(); print("ALL_DONE", flush=True)
