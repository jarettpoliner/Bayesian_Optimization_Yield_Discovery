#!/usr/bin/env python3
"""Check that the installation can reproduce the frozen method.  Exit code 0 = all required checks pass.
    python scripts/validate_install.py"""
import os, sys
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")                     # single-threaded BLAS: bit-reproducible results
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT / "src"))
import json, hashlib, importlib, tempfile, traceback
import numpy as np

RESULTS = []


def check(name, fn, required=True):
    try:
        msg = fn() or ""; RESULTS.append((name, True, required, msg)); print(f"[PASS] {name} {msg}")
    except Exception as e:
        RESULTS.append((name, False, required, repr(e))); print(f"[{'FAIL' if required else 'WARN'}] {name}: {e!r}")


def versions():
    lock = {l.split("==")[0].lower(): l.split("==")[1].strip() for l in (ROOT / "requirements-lock.txt").read_text().splitlines() if "==" in l and not l.startswith("#")}
    names = {"numpy": "numpy", "scipy": "scipy", "pandas": "pandas", "scikit-learn": "sklearn", "matplotlib": "matplotlib", "pyyaml": "yaml"}
    diff = [f"{k} {importlib.import_module(m).__version__} (lock {lock[k]})" for k, m in names.items() if k in lock and importlib.import_module(m).__version__ != lock[k]]
    if diff:
        raise RuntimeError("versions differ from requirements-lock.txt: " + ", ".join(diff) + " -- deterministic mode still works; bit-exact legacy reproduction may not")
    return "(" + ", ".join(f"{k}=={lock[k]}" for k in names if k in lock) + ")"


def data_files():
    sums = [l.split() for l in (ROOT / "data" / "SHA256SUMS").read_text().splitlines() if l.strip()]
    bad = [f for h, f in sums if hashlib.sha256((ROOT / "data" / f).read_bytes()).hexdigest() != h]
    if bad:
        raise RuntimeError(f"checksum mismatch: {bad}")
    return f"({len(sums)} files)"


def gp_reference():
    from yield_surface_discovery import DirectGP, CandidateLibrary, SHEAR_MIRROR, NO_SYMMETRY
    ref = json.loads((ROOT / "data" / "reproduction" / "gp_reference.json").read_text())
    P2 = CandidateLibrary.from_csv(ROOT / "data" / "candidates" / "expanded_physical_loading_library.csv").directions
    Q = CandidateLibrary.from_csv(ROOT / "data" / "candidates" / "arbitrary_direction_reference.csv").directions[:200]; worst = 0
    for c in ref["cases"]:
        gp = DirectGP(SHEAR_MIRROR if c["symmetry"] else NO_SYMMETRY).fit(P2[c["candidate_ids"]], np.array(c["radii"]))
        worst = max(worst, float(np.max(np.abs(gp.predict_radius(Q) / np.array(c["pred_radius"]) - 1))))
    if worst > 1e-6:
        raise RuntimeError(f"GP predictions differ from the registered implementation by {worst:.2e}")
    return f"(max relative deviation {worst:.1e})"


def registered_replay():
    sys.path.insert(0, str(ROOT / "tests")); import test_reproduction as T
    n = T.replay_registered("C_entropy_sym", 0, max_steps=6); return f"(registered physics-informed choices reproduced within tie tolerance for {n} steps)"


def minimal_controller():
    from yield_surface_discovery import CandidateLibrary, AdaptiveController, ControllerSettings, SHEAR_MIRROR, CachedOracle
    from yield_surface_discovery.oracle import analytical_oracle
    lib = CandidateLibrary.from_csv(ROOT / "data" / "candidates" / "expanded_physical_loading_library.csv", deduplicate=True)
    o = CachedOracle(analytical_oracle("von_mises", sigma_y=250.0)); R = AdaptiveController(lib, o, SHEAR_MIRROR, ControllerSettings(budget=10)).run(0)
    assert len(R.sequence) == 10 == o.calls; return f"(10 real tests, {o.calls} oracle calls)"


def moose_adapter_fake():
    sys.path.insert(0, str(ROOT / "tests")); import test_moose_adapter as M
    with tempfile.TemporaryDirectory() as d:
        err = M.run_fake_check(Path(d)); return f"(fake executable; radius error {err:.1e}; real MOOSE not required)"


if __name__ == "__main__":
    import platform; print(f"python {platform.python_version()} on {platform.platform()}")
    check("package import", lambda: importlib.import_module("yield_surface_discovery").__version__)
    check("dependency versions match requirements-lock.txt", versions, required=False)
    check("bundled data checksums", data_files)
    check("direct GP equals the registered implementation", gp_reference)
    check("registered run replay (first acquisition steps)", registered_replay)
    check("controller end-to-end (analytical oracle)", minimal_controller)
    check("MOOSE adapter against a fake executable", moose_adapter_fake)
    failed = [r for r in RESULTS if not r[1] and r[2]]; print("\nVALIDATION", "PASSED" if not failed else f"FAILED ({len(failed)} required checks)"); sys.exit(1 if failed else 0)
