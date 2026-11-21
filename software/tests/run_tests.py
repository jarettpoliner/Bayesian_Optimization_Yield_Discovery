#!/usr/bin/env python3
"""Portable-release test suite (Part D checks).  Uses only archived observations and analytic surfaces; never launches MOOSE unless
--with-moose is given AND $YSBO_MOOSE_EXE is set (then only `--version` and `--check-input`, no simulation).
    python -B software/tests/run_tests.py [--outdir DIR] [--quick] [--with-moose]
Writes <outdir>/tests_report.json; exit code 0 iff every executed test passed."""
import os, sys
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "1")
sys.dont_write_bytecode = True
import argparse, json, platform, shutil, tempfile, time, traceback, warnings
from pathlib import Path
import numpy as np, pandas as pd
SW = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(SW))
import yield_surface_bo as Y
from yield_surface_bo import _frozen as F, core, campaign, moose_truss as MT

ap = argparse.ArgumentParser(); ap.add_argument("--outdir", default=str(F.ROOT / "test_outputs")); ap.add_argument("--quick", action="store_true")
ap.add_argument("--with-moose", action="store_true"); ap.add_argument("--only", default="", help="comma-separated test-name prefixes, e.g. T06,T08")
A = ap.parse_args()
OUT = Path(A.outdir); OUT.mkdir(parents=True, exist_ok=True); RESULTS = []; EX = F.ROOT / "examples"


def test(name):
    def deco(fn):
        def run():
            t0 = time.time()
            try:
                det = fn() or {}; st = det.pop("_status", "pass")
            except Exception as e:
                st, det = "fail", dict(error=repr(e), trace=traceback.format_exc()[-1500:])
            RESULTS.append(dict(test=name, status=st, seconds=round(time.time() - t0, 1), **det)); print(f"[{st.upper():4s}] {name} ({time.time() - t0:.0f} s)", flush=True)
        run.__name__ = name; return run
    return deco


def cfg_from(example, **over):
    c = campaign.load_config(EX / example); c.update(over); return c


def must(cond, msg):
    if not cond:
        raise AssertionError(msg)


@test("T00_imports_dependencies_frozen_hashes")
def t00():
    import scipy, sklearn, matplotlib, yaml
    lock = {l.split("==")[0]: l.split("==")[1].strip() for l in (F.ROOT / "requirements-lock.txt").read_text().splitlines() if "==" in l and not l.startswith("#")}
    have = dict(numpy=np.__version__, scipy=scipy.__version__, pandas=pd.__version__, **{"scikit-learn": sklearn.__version__}, matplotlib=matplotlib.__version__, pyyaml=yaml.__version__)
    mism = {k: (v, lock.get(k)) for k, v in have.items() if lock.get(k) and lock[k] != v}
    fz = F.verify(strict=False); M = F.load()
    must(fz["ok"], f"frozen sources: {fz}")
    return dict(python=platform.python_version(), packages=have, lock_mismatches=mism, frozen=fz,
                note="version mismatches do not fail the test; bitwise reproduction is only guaranteed with the locked versions")


@test("T01_library_uniqueness_and_symmetry_checks")
def t01():
    reps = {}
    for b in ("tl1_signed_8192", "analytic_axes_4096", "b1b2_axes_207"):
        reps[b] = core.library_report(core.load_library({"builtin": b}))
    must(reps["tl1_signed_8192"]["n"] == 8192 and reps["tl1_signed_8192"]["duplicate_pairs"] == 0 and reps["tl1_signed_8192"]["antipodally_closed"], "TLIB")
    must(reps["analytic_axes_4096"]["duplicate_pairs"] == 0 and reps["analytic_axes_4096"]["antipode_fraction"] == 0.0, "analytic axes")
    must(reps["b1b2_axes_207"]["n"] == 414 and reps["b1b2_axes_207"]["antipodally_closed"], "B1/B2 library")
    lib = core.load_library({"builtin": "tl1_signed_8192"}); base = dict(n_initial=8, budget=33, library={"builtin": "tl1_signed_8192"})
    c1 = core.check_compat(dict(base, method="v2_mes", specimen=dict(kind="moose_truss", symmetry="antipodal")), lib, reps["tl1_signed_8192"])
    c2 = core.check_compat(dict(base, method="v3_signed", specimen=dict(kind="analytic", symmetry="antipodal", symmetry_validated=True, symmetry_evidence="x")), lib, reps["tl1_signed_8192"])
    c3 = core.check_compat(dict(base, method="v3_weighted", specimen=dict(kind="moose_truss", symmetry="none")), lib, reps["tl1_signed_8192"])
    c4 = core.check_compat(dict(base, method="v2_mes", specimen=dict(kind="analytic", symmetry="none")), None, reps["analytic_axes_4096"])
    c4 = core.check_compat(dict(base, method="v2_mes", specimen=dict(kind="analytic", symmetry="none")), core.load_library({"builtin": "analytic_axes_4096"}), reps["analytic_axes_4096"])
    from yield_surface_discovery.candidates import CandidateLibrary
    D = lib.directions[:20].copy(); D[5] = D[3]; dup = CandidateLibrary(D, np.arange(20), [str(i) for i in range(20)], ["x"] * 20)
    c5 = core.check_compat(dict(base, budget=10, method="v2_mes", specimen=dict(kind="analytic", symmetry="none")), dup)
    must(not c1.ok and any("symmetry_validated" in e for e in c1.errors), "unvalidated antipodal symmetry must be refused")
    must(not c2.ok, "v3_signed + antipodal must be refused"); must(c3.ok and len(c3.warnings) >= 2, "v3_weighted signed must warn (symmetric risk; physical superiority)")
    must(not c4.ok, "half-sphere library in a signed campaign must be refused"); must(not c5.ok, "duplicate directions must be refused")
    return dict(libraries=reps, refused_unvalidated_antipodal=c1.errors[0][:120], v3_weighted_signed_warnings=c3.warnings)


def analytic_smoke(method):
    c = cfg_from(f"analytic_{method}.yaml", outdir=str(OUT / f"analytic_{method}"))
    if (OUT / f"analytic_{method}").exists():
        shutil.rmtree(OUT / f"analytic_{method}")                                          # test output only (never research data)
    r = campaign.run_campaign(c, quiet=True); e = r["info"]["error_curve"]; seq = r["sequence"]
    must(r["info"]["status"] == "complete" and len(seq) == 33 and len(set(seq)) == 33, "33 unique observations")
    must(seq[:8] == r["controller"].initial_design(0), "shared initial design"); must(e[-1]["mean_rel_error"] < e[0]["mean_rel_error"], "error must decrease")
    must((OUT / f"analytic_{method}" / "reconstruction.png").stat().st_size > 10000, "figure written")
    return dict(error_curve={x["n"]: round(100 * x["mean_rel_error"], 4) for x in e}, maturity=r["info"]["maturity"], first_selected=seq[8:12])


t02 = test("T02_v2_mes_analytic_smoke")(lambda: analytic_smoke("v2_mes"))
t03 = test("T03_v3_weighted_analytic_smoke")(lambda: analytic_smoke("v3_weighted"))
t04 = test("T04_angular_maximin_analytic_smoke")(lambda: analytic_smoke("angular_maximin"))


def scores_cmp(mine, saved, as_f32):
    a = np.array([np.where(np.isfinite(x), x, np.nan) for x in mine], np.float32 if as_f32 else float); b = np.asarray(saved); b = np.where(np.isfinite(b), b, np.nan)   # saved -inf = tested
    if a.shape != b.shape:
        return dict(shape_mismatch=[list(a.shape), list(b.shape)])
    fin = np.isfinite(a) & np.isfinite(b); same_nan = bool(np.array_equal(np.isnan(a), np.isnan(b)))
    rel = float(np.max(np.abs(a[fin] - b[fin]) / np.maximum(np.abs(b[fin]), 1e-300))) if fin.any() else 0.0
    return dict(bitwise=bool(np.array_equal(a, b, equal_nan=True)), max_rel_diff=rel, nan_pattern_identical=same_nan)


@test("T05_frozen_H1_V3_regression_replay")
def t05():
    F.load()                                                                          # frozen import paths (independent of test order)
    import run_h1 as RH
    from yield_surface_discovery.oracle import AnalyticalOracle
    M = F.load(); rec = RH.POP[0]; Fh = RH.H.Yld91(rec)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = M["V3"].V3Controller(RH.C.library(), AnalyticalOracle(Fh.radius, "replay"), M["B"].V2Settings(schedule="A20", symmetry="antipodal", budget=64)).run(rec["start"])
    sv = json.load(open(F.H1D / "runs/V3/r00/design.json")); sc = scores_cmp(r["scores"], np.load(F.H1D / "runs/V3/r00/scores_f32.npy"), True)
    seq_ok = r["sequence"] == sv["axes"]; must(seq_ok, "H1 V3 sequence differs from the frozen record"); must(sc.get("max_rel_diff", 1) < 1e-5, f"scores {sc}")
    return dict(sequence_identical=seq_ok, radii_identical=[o.yield_radius for o in r["observations"]] == sv["radii"], scores=sc, n=64)


def replay_release(v):
    c = cfg_from(f"replay_{v.lower()}_v2.yaml", outdir=str(OUT / f"replay_{v.lower()}_v2"))
    if (OUT / f"replay_{v.lower()}_v2").exists():
        shutil.rmtree(OUT / f"replay_{v.lower()}_v2")
    r = campaign.run_campaign(c, quiet=True); camp = F.PS / f"release_B1_B2_BO_benchmark_20261119/data/campaigns/{v}/bo_s0"
    tests = pd.read_csv(camp / "tests.csv", float_precision="round_trip"); seq_ok = r["sequence"] == tests.axis.tolist()
    rad_ok = bool(np.array_equal(r["radii"], tests.yield_radius.to_numpy())); sc = scores_cmp(r["scores"], np.load(camp / "acquisition_scores.npy"), False)
    must(seq_ok and rad_ok, f"{v} replay sequence/radii differ"); must(sc.get("max_rel_diff", 1) < 1e-6, f"scores {sc}")
    return dict(sequence_identical=seq_ok, radii_identical=rad_ok, scores=sc, n=len(r["sequence"]), new_moose_runs=0)


t06 = test("T06_recorded_B1_MOOSE_replay_V2")(lambda: replay_release("B1"))
t07 = test("T07_recorded_B2_MOOSE_replay_V2")(lambda: replay_release("B2"))


def replay_tl1(method, saved_dir, frozen_err=None):
    c = cfg_from(f"replay_tl1_{method}.yaml", outdir=str(OUT / f"replay_tl1_{method}"))
    if (OUT / f"replay_tl1_{method}").exists():
        shutil.rmtree(OUT / f"replay_tl1_{method}")
    r = campaign.run_campaign(c, quiet=True); sv = json.load(open(saved_dir / "design.json")); seq_ok = r["sequence"] == sv["axes"]
    sc = scores_cmp(r["scores"], np.load(saved_dir / "scores_f32.npy"), True); e33 = next(x for x in r["info"]["error_curve"] if x["n"] == 33)["mean_rel_error"]
    must(seq_ok, f"TL1 {method} replay sequence differs"); must(sc.get("max_rel_diff", 1) < 1e-5, f"scores {sc}")
    out = dict(sequence_identical=seq_ok, radii_identical=bool(np.array_equal(r["radii"], sv["radii"])), scores=sc, error_n33=e33, new_moose_runs=0)
    if frozen_err is not None:
        out["frozen_evaluation_n33"] = frozen_err; out["error_matches_frozen_evaluation"] = abs(e33 / frozen_err - 1) < 1e-9; must(out["error_matches_frozen_evaluation"], "n33 error")
    return out


def frozen_metric(evdir, method):
    m = pd.read_csv(evdir / "evaluation/metrics.csv", float_precision="round_trip")
    return float(m[(m.method == method) & (m.start == 0) & (m.n == 33) & (m.surrogate == "A20")].area_mean.iloc[0])


t08 = test("T08_recorded_TL1_MOOSE_replay_V2")(lambda: replay_tl1("v2_mes", F.C3D / "runs/V2/s0", frozen_metric(F.C3D, "V2")))
t09 = test("T09_recorded_TL1_MOOSE_replay_V3_frozen")(lambda: replay_tl1("v3_weighted", F.V3X / "runs/V3F/s0", frozen_metric(F.V3X, "V3F")))
t10 = test("T10_recorded_TL1_MOOSE_replay_V3_signed")(lambda: replay_tl1("v3_signed", F.V3X / "runs/V3S/s0", frozen_metric(F.V3X, "V3S")))


@test("T11_checkpoint_save_and_resume")
def t11():
    ref_dir, cp_dir = OUT / "ckpt_reference", OUT / "ckpt_interrupted"
    for d in (ref_dir, cp_dir):
        if d.exists():
            shutil.rmtree(d)
    base = cfg_from("analytic_v3_weighted.yaml"); base["budget"] = 24; base["evaluation"] = dict(checkpoints=[8, 24], n_references=2000)
    a = campaign.run_campaign(dict(base, outdir=str(ref_dir)), quiet=True)
    c = dict(base, outdir=str(cp_dir), limits=dict(stop_file="STOP")); orig = core.AnalyticOracle.evaluate; calls = {"n": 0}

    def evaluate(self, d):                                                          # simulate an operator stopping the campaign after 13 experiments
        calls["n"] += 1
        if calls["n"] == 13:
            (cp_dir / "STOP").write_text("stop")
        return orig(self, d)
    core.AnalyticOracle.evaluate = evaluate
    try:
        b1 = campaign.run_campaign(c, quiet=True)
    finally:
        core.AnalyticOracle.evaluate = orig
    must(b1["info"]["status"] == "stopped", "campaign must stop at the STOP file"); n_saved = sum(1 for _ in open(cp_dir / "observations.jsonl"))
    (cp_dir / "STOP").unlink(); b2 = campaign.run_campaign(c, quiet=True)
    same = b2["sequence"] == a["sequence"] and np.array_equal(b2["radii"], a["radii"])
    must(same, "resumed campaign differs from the uninterrupted one"); must(b2["info"]["n_resumed_observations"] >= n_saved, "stored observations must be reused")
    return dict(stopped_after=n_saved, resumed_reused=b2["info"]["n_resumed_observations"], resumed_new=b2["info"]["n_new_observations"], identical_to_uninterrupted=same)


@test("T12_truss_specimen_workspace_and_input_dry_run")
def t12():
    ws = OUT / "truss_ws_tl1"; ws2 = OUT / "truss_ws_new"
    for d in (ws, ws2, OUT / "truss_ws_bad"):
        if d.exists():
            shutil.rmtree(d)
    orc = MT.TrussMooseOracle(ws); inp = ws / "dryrun" / "run.i"; inp.parent.mkdir(parents=True); w = orc.write_input([0.3, -0.8, 0.5], inp)
    # frozen reference: the research copy of the oracle with its own lattice writes the input for the same direction
    TOf = MT.load_oracle_module(F.TRUSS, None); TOf.MESH_DIR = F.TRUSS / "oracle" / "mesh_TL1"; ref = OUT / "truss_frozen_ref.i"
    cal = json.loads((F.TRUSS / "oracle/calibration_TL1.json").read_text()); sch, _ = TOf.schedule(np.array([0.3, -0.8, 0.5]), cal, dict(TOf.DEFAULTS))
    TOf.write_input([0.3, -0.8, 0.5], ref, C=np.asarray(cal["C"]), Sbar=cal["Sbar"], **sch)
    norm = lambda t: "\n".join(l for l in t.splitlines() if "mesh.e" not in l and "constraint_matrix" not in l)
    same = norm(inp.read_text()) == norm(ref.read_text()); must(same, "workspace input differs from the frozen oracle's input")
    v_new = MT.prepare_workspace(ws2, EX / "specimens/example_new_truss_specimen.json")
    try:
        MT.prepare_workspace(OUT / "truss_ws_bad", EX / "specimens/example_INVALID_unequal_bars.json"); bad = "accepted (WRONG)"
    except ValueError as e:
        bad = str(e)[:160]
    must("equal bar lengths" in bad, "invalid specimen must be rejected")
    return dict(input_identical_to_frozen_oracle_except_mesh_paths=same, predicted_first_yield=w["predicted_first_yield"], new_specimen_check=v_new, invalid_specimen=bad,
                note="dry run: no MOOSE executed")


@test("T13_live_moose_adapter_dependency_detection")
def t13():
    old = os.environ.pop("YSBO_MOOSE_EXE", None)
    try:
        none = MT.find_moose(None) if not shutil.which("combined-opt") else dict(found=None, note="combined-opt on PATH")
        bogus = MT.find_moose("/nonexistent/combined-opt"); nonexec = MT.find_moose(str(EX / "specimens/example_new_truss_specimen.json"))
    finally:
        if old:
            os.environ["YSBO_MOOSE_EXE"] = old
    must(bogus["found"] is False and nonexec["found"] is False, "invalid executables must be rejected")
    out = dict(no_executable=none, nonexistent=bogus, not_executable=nonexec)
    if A.with_moose and os.environ.get("YSBO_MOOSE_EXE"):
        exe = os.environ["YSBO_MOOSE_EXE"]; out["configured"] = MT.find_moose(exe); out["version"] = MT.moose_version(exe, MT._env())
        orc = MT.TrussMooseOracle(OUT / "truss_ws_check", None, exe); inp = OUT / "truss_ws_check" / "check" / "run.i"; inp.parent.mkdir(parents=True, exist_ok=True)
        orc.write_input([1.0, 0.0, 0.0], inp); import subprocess
        p = subprocess.run([exe, "-i", "run.i", "--check-input", "--color", "off"], cwd=inp.parent, capture_output=True, text=True, timeout=600, env=MT._env())
        out["check_input"] = dict(returncode=p.returncode, tail=(p.stdout + p.stderr)[-500:]); must(p.returncode == 0, "MOOSE --check-input failed")
    else:
        out["_status"] = "pass"; out["live_check"] = "not executed (no --with-moose / YSBO_MOOSE_EXE): detection logic only"
    return out


@test("T14_external_command_adapter_plumbing")
def t14():
    d = OUT / "external_command_v2"
    if d.exists():
        shutil.rmtree(d)
    r = campaign.run_campaign(cfg_from("external_command_v2.yaml", outdir=str(d)), quiet=True)
    must(r["info"]["status"] == "complete" and len(r["sequence"]) == 16, "external-command campaign")
    return dict(n=len(r["sequence"]), note="analytic stand-in command (not MOOSE): verifies the per-direction command/result.json protocol")


@test("T15_independence_from_original_workstation")
def t15():
    mods = F.first_party_module_files(); root = str(F.ROOT)
    fp = {k: v for k, v in mods.items() if any(s in v for s in ("yield_surface_bo", "bomes_v2", "v3_lrv", "v3s_lrv", "yield_surface_discovery", "run_bo_moose", "cs_common",
                                                                 "surfaces_h1", "run_h1", "truss_oracle", "lattice_def"))}
    outside = {k: v for k, v in fp.items() if not v.startswith(root)}
    home = Path.home(); orig = home / "projects" / "yield_surface_moose"
    must(not outside, f"first-party modules loaded from outside the release: {outside}")
    return dict(release_root=root, first_party_modules=len(fp), outside_release=outside, HOME=str(home), original_project_path_visible_under_HOME=orig.exists(),
                note="run via verify_installation.sh with HOME redirected to an empty directory for a strict check")


tests = [t00, t01, t02, t03, t04, t05, t06, t07, t08, t09, t10, t11, t12, t13, t14, t15]
if A.quick:
    tests = [t00, t01, t03, t06, t08, t11, t15]
if A.only:
    tests = [t for t in tests if any(t.__name__.startswith(p.strip()) for p in A.only.split(","))]
for t in tests:
    t()
rep = dict(created=time.strftime("%Y-%m-%dT%H:%M:%S"), release_root=str(F.ROOT), python=platform.python_version(), platform=platform.platform(),
           summary={s: sum(r["status"] == s for r in RESULTS) for s in ("pass", "fail", "skip")}, results=RESULTS)
(OUT / "tests_report.json").write_text(json.dumps(rep, indent=1, default=str)); print(json.dumps(rep["summary"]))
sys.exit(0 if rep["summary"]["fail"] == 0 else 1)
