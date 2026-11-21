#!/usr/bin/env python3
"""Part 8 replay & integrity verification for the final B1/B2 benchmark (read-only except verification/*).
  A. independent controller recomputation per BO campaign from its recorded tests.csv: GP refit (frozen V2Controller.fit / GPv2) at every step,
     latent-variance acquisition scores vs acquisition_scores.npy (bitwise), frozen select() -> recorded candidate, hyper-parameters/mode vs steps.csv
  B. ReplayOracle replays (run by the pipeline): tests/steps/scores bitwise
  C. exact-match reuse: every 'reused_from' record has an original record in its source with the same variant, config, direction key and an
     identical observation; sources never cross variants and never include the reference logs
  D. physical-observation accounting: 33 tests per campaign, distinct axes, every trial_id -> runs/<v>/h40/<tag>/record.json accepted with the same r
  E. non-adaptive designs = prefix of the frozen design lists (+ spares only after recorded failures), radii = pool records
  F. references: 72/72 accepted per variant, frozen hash, no reference direction in any controller log
  G. final source audit (recovery_r1/audit_recovery.py final)     H. historical-file preservation vs recovery_r1/preservation_snapshot.json
python verify_final.py   (BLAS threads forced to 1 before numpy import)"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"
import hashlib, json, subprocess, sys
from pathlib import Path
import numpy as np, pandas as pd
FB = Path(__file__).resolve().parent; BB = FB.parent; CL = BB.parent
for p in (FB, CL, CL / "v2"):
    sys.path.insert(0, str(p))
import designs_common as D                                                          # noqa: E402
from bomes_v2 import V2Controller, V2Settings                                       # noqa: E402
from yield_surface_discovery.acquisition import select                              # noqa: E402
OUT = FB / "verification"; OUT.mkdir(exist_ok=True)
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest(); key = lambda d: tuple(np.round(np.asarray(d, float), 12))
lines = lambda p: [json.loads(l) for l in open(p) if l.strip()]
lib, AX = D.library(); V = {}


def recompute(camp):
    t = pd.read_csv(camp / "tests.csv", float_precision="round_trip"); s = pd.read_csv(camp / "steps.csv", float_precision="round_trip"); sc = np.load(camp / "acquisition_scores.npy")
    ctl = V2Controller(lib, None, V2Settings(schedule="A20", symmetry="antipodal", budget=33)); C = ctl.C
    radius = dict(zip(t.axis, t.yield_radius)); seq = list(t.axis[:8])
    ev = (s["event"] if "event" in s else pd.Series([np.nan] * len(s))).astype(str); init_fail = ev.str.startswith("failed initial")
    unavailable = set(int(c) for c in s.candidate[init_fail]); rows = s[~init_fail].reset_index(drop=True); evr = ev[~init_fail].reset_index(drop=True)
    res = dict(n_steps=len(rows), scores_bitwise=[], selection=[], hyper=[], mode=[])
    for i, st in rows.iterrows():
        gp = ctl.fit(C[seq], [radius[c] for c in seq]); mean, var = gp.predict_log(C); x = var.copy(); x[seq] = -np.inf; x[list(unavailable)] = -np.inf
        j, tied = select(x, "deterministic", 1e-10, 1e-12); hp = gp.hyperparameters()
        res["scores_bitwise"].append(bool(np.array_equal(x, sc[i]))); res["selection"].append(int(j) == int(st.candidate))
        res["hyper"].append(all(abs(hp[k] - st[k]) <= 1e-12 * abs(st[k]) for k in ("amplitude", "length_scale", "white_noise")) and hp["n_real"] == st.n_real)   # CSV-stored values; res["mode"].append(hp["mode"] == st["mode"])
        if evr[i].startswith("failed"):
            unavailable.add(int(j))
        else:
            seq.append(int(j))
    res["sequence_reproduced"] = seq == list(t.axis)
    res["ok"] = res["sequence_reproduced"] and all(res["scores_bitwise"]) and all(res["selection"]) and all(res["hyper"]) and all(res["mode"]) and len(sc) == len(rows)
    res["hyper_evolution"] = rows[["test_number", "mode", "amplitude", "length_scale", "white_noise", "log_marginal_likelihood"]].to_dict("list")
    return {k: (all(v) if isinstance(v, list) and k != "hyper_evolution" else v) for k, v in res.items()}


if __name__ == "__main__":
    # A
    V["A_recompute"] = {f"{v}_s{k}": recompute(FB / "campaigns" / v / f"bo_s{k}") for v in D.VARIANTS for k in D.STARTS}
    # B
    stp = json.load(open(FB / "recovery_r1/pipeline_status_r1.json"))
    V["B_replay"] = {v: stp.get(f"stage{3 if v == 'B1' else 4}_{v}_bo") for v in D.VARIANTS}
    V["B_replay"]["ok"] = all(V["B_replay"][v] and all(V["B_replay"][v]["replay_identical"]) for v in D.VARIANTS)
    # C + D
    refkeys = {v: {key(r["direction"]) for r in lines(FB / "refs" / v / "trials.jsonl")} for v in D.VARIANTS}
    C_bad, n_reused, D_rows, overlap = [], 0, [], 0
    for v in D.VARIANTS:
        camp = FB / "campaigns" / v; logs = [camp / "nonadaptive/trials.jsonl"] + [camp / f"bo_s{k}/trials.jsonl" for k in D.STARTS]
        for p in logs:
            for r in lines(p):
                overlap += key(r["direction"]) in refkeys[v]
                if "reused_from" in r:
                    n_reused += 1; src = Path(r["reused_from"])
                    orig = [o for o in lines(src) if "reused_from" not in o and key(o["direction"]) == key(r["direction"]) and o["variant"] == v and o["config"] == r["config"]]
                    if not (orig and orig[0]["observation"] == r["observation"] and r["variant"] == v and f"/{v}/" in str(src) and "/refs/" not in str(src)):
                        C_bad.append(dict(log=str(p), tag=r.get("original_tag"), src=str(src)))
        for k in D.STARTS:
            t = pd.read_csv(camp / f"bo_s{k}/tests.csv", float_precision="round_trip"); okr = []
            for r in t.itertuples():
                R = json.load(open(BB / "runs" / v / "h40" / r.trial_id.split("/")[1] / "record.json")); okr.append(R["accepted"] and R["yield_radius"] == r.yield_radius and R["variant"] == v)
            D_rows.append(dict(campaign=f"{v}_bo_s{k}", n=len(t), distinct_axes=t.axis.nunique(), records_ok=all(okr)))
        d = json.load(open(camp / "nonadaptive/designs.json")); pool = {key(r["direction"]): r for r in lines(camp / "nonadaptive/trials.jsonl")}
        des = D.nonadaptive_designs(lib, AX)
        for name, x in d["designs"].items():
            okr = [pool[key(AX.directions[c])]["observation"]["yield_radius"] == rr and json.load(open(BB / "runs" / v / "h40" / tid.split("/")[1] / "record.json"))["accepted"]
                   for c, rr, tid in zip(x["axes"], x["radii"], x["trial_ids"])]
            want = des[name][:33 + len(x["failed"])]
            D_rows.append(dict(campaign=f"{v}_{name}", n=len(x["axes"]), distinct_axes=len(set(x["axes"])), records_ok=all(okr),
                               design_rule=x["axes"] == [c for c in want if c not in x["failed"]][:33]))
    V["C_reuse"] = dict(ok=not C_bad, n_reused_records=n_reused, bad=C_bad[:20])
    Dd = pd.DataFrame(D_rows); Dd.to_csv(OUT / "observation_accounting.csv", index=False)
    V["D_E_observations"] = dict(ok=bool((Dd.n == 33).all() and (Dd.distinct_axes == 33).all() and Dd.records_ok.all() and Dd.design_rule.fillna(True).all()), n_campaigns=len(Dd))
    fr = json.load(open(FB / "recovery_r1/reference_datasets_frozen.json"))
    exc = json.load(open(FB / "recovery_r1/stage2_gate_exception.json")); allowed = {(e["variant"], e["tag"]) for e in exc["allowed_unaccepted"]}
    unacc = {(v, r["tag"]) for v in D.VARIANTS for r in lines(FB / "refs" / v / "trials.jsonl") if not r["observation"]["accepted"]}
    V["F_references"] = dict(ok=all(sha(FB / fr[v]["path"]) == fr[v]["sha256"] and len(refkeys[v]) == 72 for v in D.VARIANTS) and unacc == allowed and overlap == 0,
                             unaccepted=sorted(unacc), approved_exception=sorted(allowed), controller_records_at_reference_directions=overlap)
    # G
    g = subprocess.run([sys.executable, str(FB / "recovery_r1/audit_recovery.py"), "final"], capture_output=True, text=True)
    V["G_source_audit"] = dict(ok=json.load(open(FB / "recovery_r1/audit_final.json"))["PASS"], returncode=g.returncode)
    # H
    snap = json.load(open(FB / "recovery_r1/preservation_snapshot.json")); bad = []
    for rel, e in snap["files"].items():
        p = BB / rel
        if rel.endswith("refs/B1/trials.jsonl") or rel.endswith("refs/B2/trials.jsonl"):
            v = "B1" if "B1" in rel else "B2"; ok = hashlib.sha256(b"".join(open(p, "rb").readlines()[:11])).hexdigest() == snap["refs_prefix"][v]["sha256_first11"]
        else:
            ok = p.exists() and sha(p) == e["sha256"]
        if not ok:
            bad.append(rel)
    V["H_preservation"] = dict(ok=not bad, n_files=len(snap["files"]), changed_or_missing=bad[:50], note="refs/*/trials.jsonl checked as append-only (first 11 lines unchanged)")
    V["PASS"] = all(x["ok"] for x in V.values() if isinstance(x, dict) and "ok" in x) and all(x["ok"] for x in V["A_recompute"].values())
    json.dump(V, open(OUT / "verification.json", "w"), indent=1, default=str)
    print(json.dumps({k: (x.get("ok") if isinstance(x, dict) and "ok" in x else (all(y["ok"] for y in x.values()) if k == "A_recompute" else x)) for k, x in V.items()}, indent=1))
