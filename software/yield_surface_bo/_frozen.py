"""Locate the frozen research code inside the portable release and import it READ-ONLY.

Release layout (mirrors the original ~/projects tree so every relative path inside the frozen code keeps working unchanged):
    <ROOT>/research/yield_surface_moose/plane_stress/...      complete project snapshot (frozen controllers, oracles, data, campaigns)
    <ROOT>/research/yield_surface_discovery_release/...       vendored first-party V1 package (CandidateLibrary, GP helpers, maximin)
    <ROOT>/software/yield_surface_bo/                         this public interface (wrappers only; no algorithm code)
ROOT = $YSBO_ROOT if set, else two levels above this file.

Frozen modules are verified against FROZEN_SOURCES.json (sha256 recorded from the original workstation at packaging time) before import.
Several frozen modules insert Path.home()/projects/... into sys.path; on another machine those paths do not exist and are harmless, because
the vendored copies below are placed first on sys.path.
"""
import hashlib, json, os, sys
from pathlib import Path

ROOT = Path(os.environ.get("YSBO_ROOT", Path(__file__).resolve().parents[2])).resolve()
RESEARCH = ROOT / "research"
PS = RESEARCH / "yield_surface_moose" / "plane_stress"
CL = PS / "closed_loop"
V1_SRC = RESEARCH / "yield_surface_discovery_release" / "software" / "yield_surface_discovery" / "src"
TRUSS = CL / "truss_lattice_benchmark_2026"
V3D = CL / "response_adaptive_v3_2026"
V3X = TRUSS / "final_comparison_2026" / "v3transfer_20261120"
C3D = TRUSS / "final_comparison_2026" / "campaign3_20261120"
H1D = V3D / "holdout_h1_2026"
CSD = CL / "complexity_scaling_2026"
FROZEN_JSON = Path(__file__).resolve().parent / "FROZEN_SOURCES.json"
sys.dont_write_bytecode = True                                    # never write caches into the frozen snapshot

_PATHS = [V1_SRC, CL / "v2", V3D, V3X, CSD, H1D, CL]


def setup_paths():
    for p in reversed(_PATHS):
        sp = str(p)
        if sp in sys.path:
            sys.path.remove(sp)
        sys.path.insert(0, sp)


def sha256(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def verify(strict=True):
    """compare every frozen source listed in FROZEN_SOURCES.json with the extracted copy.  returns dict(ok, mismatches, missing)."""
    rec = json.loads(FROZEN_JSON.read_text())["files"]; bad, missing = [], []
    for rel, h in rec.items():
        p = ROOT / rel
        if not p.exists():
            missing.append(rel)
        elif sha256(p) != h:
            bad.append(rel)
    out = dict(ok=not bad and not missing, n=len(rec), mismatches=bad, missing=missing)
    if strict and not out["ok"]:
        raise RuntimeError(f"frozen sources differ from the release record: {out}")
    return out


_MODS = {}


def load():
    """import the frozen modules (after hash verification).  returns a namespace dict."""
    if _MODS:
        return _MODS
    if not PS.exists() or not V1_SRC.exists():
        raise RuntimeError(f"release snapshot not found under {RESEARCH} (set YSBO_ROOT to the extracted Yield_Surface_BO directory)")
    verify(strict=True); setup_paths()
    import bomes_v2, v3_lrv, v3s_lrv                                                  # noqa: E401  frozen controllers
    from yield_surface_discovery import candidates, initialization, oracle, acquisition   # frozen V1 package
    _MODS.update(B=bomes_v2, V3=v3_lrv, V3S=v3s_lrv, candidates=candidates, initialization=initialization, oracle=oracle, acquisition=acquisition)
    return _MODS


def first_party_module_files():
    """files of every loaded module that lives inside ROOT (used by the independence test)."""
    out = {}
    for name, m in list(sys.modules.items()):
        f = getattr(m, "__file__", None)
        if f:
            out[name] = str(Path(f).resolve())
    return out
