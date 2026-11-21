#!/usr/bin/env python3
"""Integrity audit of holdout H1 (read-only except integrity_audit.json).   python audit_h1.py  (ysd interpreter)"""
import hashlib, json, subprocess, sys, time, warnings
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; V3D = HERE.parent; CL = V3D.parent; PS = CL.parent; sys.dont_write_bytecode = True; sys.path.insert(0, str(HERE))
import run_h1 as RH                                                                  # noqa: E402  (imports only)
from yield_surface_discovery.oracle import AnalyticalOracle                          # noqa: E402
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest(); chk = lambda cwd, f: subprocess.run(["sha256sum", "--quiet", "-c", f], cwd=cwd, capture_output=True).returncode
out = dict(protocol_freeze_rc=chk(HERE, "PROTOCOL.sha256"))
t_f = (HERE / "PROTOCOL.sha256").stat().st_mtime; t_r = min(p.stat().st_mtime for p in (HERE / "runs").rglob("*") if p.is_file())
out["frozen_before_first_run"] = dict(freeze=time.strftime("%F %T", time.localtime(t_f)), first_run_file=time.strftime("%F %T", time.localtime(t_r)), ok=t_f < t_r)
out["v3_is_frozen_v3_0"] = sha(V3D / "v3_lrv.py") == RH.V3_FROZEN
infos = sorted((HERE / "runs").glob("*/r*/run_info.json")); bad, cur = [], {}
for ri in infos:
    I = json.load(open(ri))
    for p, h in I["frozen_sources"].items():
        cur.setdefault(p, sha(p)); bad += [(str(ri), p)] if cur[p] != h else []
out["sources"] = dict(campaigns=len(infos), files=len(cur), mismatches=bad)
LIB = RH.C.library().directions; dz = []
for i in range(24):
    F = RH.H.Yld91(RH.POP[i]); first8 = None
    for arm in RH.ARMS:
        d = json.load(open(HERE / f"runs/{arm}/r{i:02d}/design.json")); X = LIB[d["axes"]]; fn = RH.J2 if arm == "V3_ISO" else F.radius
        exp = np.r_[fn(X[:8]), [fn(x[None])[0] for x in X[8:]]] if arm in ("V3", "V2", "V3_ISO") else fn(X)
        first8 = first8 or d["axes"][:8]
        dz.append(dict(arm=arm, i=i, n64=len(d["axes"]) == 64, unique=len(set(d["axes"])) == 64, radii_exact=bool(np.array_equal(np.asarray(d["radii"]), exp)), same_init=d["axes"][:8] == first8))
out["designs"] = dict(n=len(dz), all_ok=all(all(v for k, v in x.items() if k not in ("arm", "i")) for x in dz))
rep = {}
for arm in ("V3", "V2"):
    F = RH.H.Yld91(RH.POP[0]); st = RH.B.V2Settings(schedule="A20", symmetry="antipodal", budget=64)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore"); r = (RH.V3.V3Controller if arm == "V3" else RH.B.V2Controller)(RH.C.library(), AnalyticalOracle(F.radius, "replay"), st).run(RH.POP[0]["start"])
    sc = np.load(HERE / f"runs/{arm}/r00/scores_f32.npy"); mine = np.array([np.where(np.isfinite(x), x, np.nan) for x in r["scores"]], np.float32)
    rep[arm] = bool(r["sequence"] == json.load(open(HERE / f"runs/{arm}/r00/design.json"))["axes"] and np.array_equal(mine, sc, equal_nan=True))
out["replay_bitwise_r00"] = rep
out["references_sha_ok"] = hashlib.sha256(np.ascontiguousarray(np.load(HERE / "references_h1.npy")).tobytes()).hexdigest() == json.load(open(HERE / "population.json"))["references_sha256"]
out["earlier_work"] = dict(V3_pilot_manifest=chk(V3D, "MANIFEST.sha256"), V3_pilot_protocol=chk(V3D, "PROTOCOL.sha256"), methodology_diagnosis=chk(CL / "methodology_diagnosis_2026", "MANIFEST.sha256"),
                           DBU=chk(CL / "adaptivity_dbu_2026", "MANIFEST.sha256"), LSF_v2=chk(CL / "adaptivity_stress_test_2026", "MANIFEST.sha256"),
                           LSF_v2_addendum=chk(CL / "adaptivity_stress_test_2026", "MANIFEST_ADDENDUM_1.sha256"), release_B1_B2=chk(PS / "release_B1_B2_BO_benchmark_20261119", "MANIFEST.sha256"),
                           release_closeout=chk(PS / "release_closeout_20261119", "MANIFEST.sha256"), complexity_protocol=chk(CL / "complexity_scaling_2026", "PROTOCOL.sha256"))
t0 = min(p.stat().st_mtime for p in HERE.rglob("*") if p.is_file()) - 1; CONC = "closed_loop/truss_lattice_benchmark_2026/"
ch = [str(p.relative_to(PS)) for p in PS.rglob("*") if p.is_file() and p.stat().st_mtime > t0 and HERE not in p.parents and "__pycache__" not in p.parts]
out["outside_modified_since_start"] = dict(start=time.strftime("%F %T", time.localtime(t0)), n_excluding_concurrent_truss_study=len([f for f in ch if not f.startswith(CONC)]),
                                           files=[f for f in ch if not f.startswith(CONC)][:20], concurrent_truss_study_files=len([f for f in ch if f.startswith(CONC)]))
out["PASS"] = bool(out["protocol_freeze_rc"] == 0 and out["frozen_before_first_run"]["ok"] and out["v3_is_frozen_v3_0"] and not bad and out["designs"]["all_ok"] and all(rep.values())
                   and out["references_sha_ok"] and all(v == 0 for v in out["earlier_work"].values()) and out["outside_modified_since_start"]["n_excluding_concurrent_truss_study"] == 0)
json.dump(out, open(HERE / "integrity_audit.json", "w"), indent=1, default=str); print(json.dumps(out, indent=1, default=str))
