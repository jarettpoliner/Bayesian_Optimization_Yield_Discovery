#!/usr/bin/env python3
"""Part 9 deliverables: consolidated experiment table/JSON, per-material and per-method histories, reference datasets, runtime/cost.
Read-only on all inputs; writes deliverables/*."""
import json, sys, time
from pathlib import Path
import numpy as np, pandas as pd
FB = Path(__file__).resolve().parent; BB = FB.parent; sys.path.insert(0, str(FB))
import designs_common as D                                                          # noqa: E402
OUT = FB / "deliverables"; OUT.mkdir(exist_ok=True); (OUT / "histories").mkdir(exist_ok=True)
lines = lambda p: [json.loads(l) for l in open(p) if l.strip()]; key = lambda d: tuple(np.round(np.asarray(d, float), 12))
rows, uses = [], []
for v in D.VARIANTS:
    logs = {"reference": FB / "refs" / v / "trials.jsonl", "nonadaptive_pool": FB / "campaigns" / v / "nonadaptive/trials.jsonl"} | \
           {f"bo_s{k}": FB / "campaigns" / v / f"bo_s{k}/trials.jsonl" for k in D.STARTS}
    for src, p in logs.items():
        for r in lines(p):
            o = r["observation"]
            if "reused_from" in r:
                uses.append(dict(variant=v, log=src, original_tag=r["original_tag"], reused_from=str(Path(r["reused_from"]).relative_to(FB)), yield_radius=o.get("yield_radius"))); continue
            R = json.load(open(BB / "runs" / v / "h40" / r["tag"] / "record.json")) if (BB / "runs" / v / "h40" / r["tag"] / "record.json").exists() else {}
            ay = R.get("at_yield", {})
            rows.append(dict(variant=v, stage="reference" if src == "reference" else ("nonadaptive" if src == "nonadaptive_pool" else "bo"), executed_by=src, tag=r["tag"], status="completed" if o.get("accepted") else "failed",
                             z1=r["direction"][0], z2=r["direction"][1], z3=r["direction"][2], yield_radius=o.get("yield_radius"), log_yield_radius=o.get("log_yield_radius"),
                             accepted=o.get("accepted"), converged=o.get("converged"), worker_wall_s=r["worker_wall_s"], started_unix=r["started_unix"],
                             eqps_rve_avg=ay.get("eqps_rve_avg"), eqps_m_avg=ay.get("eqps_m_avg"), eqps_i_avg=ay.get("eqps_i_avg"), pvf_m=ay.get("pvf_m"), pvf_i=ay.get("pvf_i"),
                             inclusion_share=(R["volumes"]["f_i"] * ay["eqps_i_avg"] / ay["eqps_rve_avg"]) if ay else None,
                             r1_matrix=R.get("first_local_yield_radius_probe", {}).get("matrix"), r1_inclusion=R.get("first_local_yield_radius_probe", {}).get("inclusion"),
                             first_plastic_matrix=R.get("first_plastic_ramp_step_radius", {}).get("matrix"), first_plastic_inclusion=R.get("first_plastic_ramp_step_radius", {}).get("inclusion"),
                             nl_iterations=R.get("nonlinear_iterations"), n_steps=R.get("n_steps"), flags=";".join(R.get("flags", [])), error=o.get("error")))
led = pd.read_csv(FB / "recovery_r1/ledger_references.csv")
for r in led[led.status == "interrupted"].itertuples():
    rows.append(dict(variant=r.variant, stage="reference", executed_by="reference (pre-reboot attempt)", tag=r.tag, status="interrupted_by_reboot", worker_wall_s=r.partial_wall_s,
                     error=f"killed by the 13:11 reboot during {r.stages_present}; superseded by the recovery_r1 run of the same design row"))
E = pd.DataFrame(rows); E.to_csv(OUT / "experiments_all.csv", index=False, float_format="%.17g"); U = pd.DataFrame(uses); U.to_csv(OUT / "reuse_records.csv", index=False)
for v in D.VARIANTS:
    E[E.variant == v].to_csv(OUT / "histories" / f"experiments_{v}.csv", index=False, float_format="%.17g")
    for st in ("references",):
        pd.read_csv(FB / "evaluation" / f"references_{v}.csv").to_csv(OUT / f"reference_dataset_{v}.csv", index=False)
    d = json.load(open(FB / "campaigns" / v / "nonadaptive/designs.json"))["designs"]
    for meth, names in (("BO", [f"bo_s{k}" for k in D.STARTS]), ("AXMM", [f"axmm_s{k}" for k in D.STARTS]), ("RANDOM", [f"random_{k}" for k in D.RANDOM_SEEDS])):
        H = []
        for nm in names:
            if meth == "BO":
                t = pd.read_csv(FB / "campaigns" / v / nm / "tests.csv", float_precision="round_trip"); t.insert(0, "campaign", nm); H.append(t)
            else:
                x = d[nm]; H.append(pd.DataFrame(dict(campaign=nm, test_number=np.arange(1, 34), axis=x["axes"], z1=[q[0] for q in x["directions"]], z2=[q[1] for q in x["directions"]],
                                                      z3=[q[2] for q in x["directions"]], yield_radius=x["radii"], trial_id=x["trial_ids"])))
        pd.concat(H).to_csv(OUT / "histories" / f"campaigns_{v}_{meth}.csv", index=False, float_format="%.17g")
ok = E[E.status != "interrupted_by_reboot"]
cost = dict(executions_total=int(len(E)), completed=int((E.status == "completed").sum()), failed=int((E.status == "failed").sum()), interrupted=int((E.status == "interrupted_by_reboot").sum()),
            unique_completed_or_failed_executions=int(len(ok)), distinct_variant_directions=int(ok.groupby("variant").apply(lambda g: len({key(r) for r in g[["z1", "z2", "z3"]].to_numpy()})).sum()),
            reuse_records=int(len(U)), cpu_h_completed_and_failed=float(ok.worker_wall_s.sum() / 3600), cpu_h_interrupted=float(E[E.status == "interrupted_by_reboot"].worker_wall_s.sum() / 3600),
            per_stage={f"{v}/{s}": dict(n=int(len(g)), cpu_h=float(g.worker_wall_s.sum() / 3600), mean_s=float(g.worker_wall_s.mean()), median_s=float(g.worker_wall_s.median()))
                       for (v, s), g in ok.groupby(["variant", "stage"])},
            first_start=time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(ok.started_unix.min())),
            last_finish=time.strftime("%Y-%m-%d %H:%M:%S", time.localtime((ok.started_unix + ok.worker_wall_s).max())),
            wall_h_first_start_to_last_finish=float(((ok.started_unix + ok.worker_wall_s).max() - ok.started_unix.min()) / 3600))
cost["cpu_h_total_incl_interrupted"] = cost["cpu_h_completed_and_failed"] + cost["cpu_h_interrupted"]
pre = ok[ok.started_unix < 1795111800]; post = ok[ok.started_unix >= 1795111800]
cost["concurrency"] = {"22-way (pre-reboot wave 1)": dict(n=int(len(pre)), mean_s=float(pre.worker_wall_s.mean())), "12-way (recovery)": dict(n=int(len(post)), mean_s=float(post.worker_wall_s.mean()))}
json.dump(cost, open(OUT / "runtime_cost.json", "w"), indent=1)
json.dump(dict(experiments=E.replace({np.nan: None}).to_dict("records"), reuse=uses, cost=cost), open(OUT / "experiments_all.json", "w"), indent=0, default=float)
print(json.dumps(cost, indent=1))
