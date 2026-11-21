#!/usr/bin/env python3
"""Stage 3/4: axial-maximin (9 starts) and random-unique-axis (10 seeds) designs for one variant, 33 successful physical tests each.
Response-INDEPENDENT: the design lists are fixed before any measurement (designs_common.nonadaptive_designs); a failed trial is recorded and
replaced by the design's next spare.  One pooled BiphasicOracle runs each distinct axis once (exact reuse within the pool).  BLIND.
Output: campaigns/<variant>/nonadaptive/{trials.jsonl, designs.json}"""
import json, sys, time
from pathlib import Path
FB = Path(__file__).resolve().parent; sys.path.insert(0, str(FB))
import designs_common as D                                                          # noqa: E402
from bench_oracle import BiphasicOracle                                             # noqa: E402
v = sys.argv[1]; par = int(sys.argv[2]) if len(sys.argv) > 2 else 22
lib, AX = D.library(); des = D.nonadaptive_designs(lib, AX); out = FB / "campaigns" / v / "nonadaptive"; out.mkdir(parents=True, exist_ok=False)
orc = BiphasicOracle(v, f"FB_{v}_pool", out / "trials.jsonl", (), par); done = {}; want = {k: s[:D.BUDGET] for k, s in des.items()}; t0 = time.time()
while True:
    todo = sorted({c for s in want.values() for c in s if c not in done})
    if todo:
        for c, o in zip(todo, orc.evaluate_many(AX.directions[todo])):
            done[c] = o
    more = False
    for k, seq in des.items():
        if sum(done[c].success for c in want[k]) < D.BUDGET:
            if len(want[k]) >= len(seq):
                raise SystemExit(f"{k}: spares exhausted")
            want[k] = want[k] + [seq[len(want[k])]]; more = True
    if not more:
        break
res = dict(variant=v, budget=D.BUDGET, wall_s=time.time() - t0, n_new=sum("reused_from" not in r for r in orc.records), n_failed=sum(not o.success for o in done.values()), designs={})
for k in des:
    good = [c for c in want[k] if done[c].success][:D.BUDGET]
    res["designs"][k] = dict(axes=good, directions=[AX.directions[c].tolist() for c in good], radii=[done[c].yield_radius for c in good],
                             trial_ids=[done[c].metadata.get("trial_id") for c in good], failed=[c for c in want[k] if not done[c].success])
(out / "designs.json").write_text(json.dumps(res, indent=1, default=float)); print(f"done {v}: {res['n_new']} new, {res['n_failed']} failed, {res['wall_s']:.0f} s", flush=True)
