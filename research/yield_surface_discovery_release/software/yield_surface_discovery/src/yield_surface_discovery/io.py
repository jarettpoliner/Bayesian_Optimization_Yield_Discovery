"""Configuration (YAML) and file I/O.  Relative paths in a config file are resolved against the config file's directory."""
import json, os, platform, sys
from pathlib import Path
import numpy as np, pandas as pd
from .candidates import CandidateLibrary
from .controller import ControllerSettings, AdaptiveController
from .symmetry import symmetry_from_name
from .oracle import analytical_oracle, CachedOracle

PATH_KEYS = {"path", "directions_path", "directory", "work_dir", "input_template", "command_script", "expected_path"}


def package_root():
    return Path(__file__).resolve().parents[2]


def default_data_dir():
    return Path(os.environ.get("YSD_DATA_DIR", package_root() / "data"))


def _resolve(obj, base):
    if isinstance(obj, dict):
        return {k: (str((base / v).resolve()) if (k in PATH_KEYS and isinstance(v, str) and not Path(v).is_absolute()) else _resolve(v, base)) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_resolve(v, base) for v in obj]
    return obj


def load_config(path):
    import yaml
    path = Path(path).resolve(); cfg = yaml.safe_load(path.read_text()); cfg = _resolve(cfg, path.parent); cfg["_config_path"] = str(path); validate(cfg); return cfg


FROZEN = dict(kernel="matern32", target="log_radius", initialization="angular_maximin", acquisition="predictive_variance")


def validate(cfg):
    g, i, a = cfg.get("gp", {}), cfg.get("initialization", {}), cfg.get("acquisition", {})
    if g.get("kernel", "matern32") != "matern32" or g.get("target", "log_radius") != "log_radius":
        raise ValueError("the frozen method uses gp.kernel = matern32 on gp.target = log_radius")
    if i.get("method", "angular_maximin") != "angular_maximin" or a.get("method", "predictive_variance") != "predictive_variance":
        raise ValueError("the frozen method uses initialization.method = angular_maximin and acquisition.method = predictive_variance")


def build_library(cfg):
    c = cfg["candidate_library"]; return CandidateLibrary.from_csv(c["path"], deduplicate=c.get("deduplicate", True), dedup_tol=float(c.get("dedup_tolerance", 1e-9)))


def build_symmetry(cfg):
    return symmetry_from_name("shear_mirror" if cfg.get("symmetry", {}).get("shear_mirror", False) else "none")


def build_settings(cfg):
    i, a, b = cfg.get("initialization", {}), cfg.get("acquisition", {}), cfg.get("budget", {})
    return ControllerSettings(n_initial=int(i.get("n_initial", 8)), budget=int(b.get("total_real_tests", 28)), start_rule=i.get("start_rule", "seeded"),
                              seed_base=int(i.get("seed_base", 7000)), first_index=int(i.get("first_index", 0)), maximin_ties=i.get("ties", "deterministic"),
                              acquisition_ties=a.get("ties", "deterministic"), tie_rel=float(a.get("tie_relative_tolerance", 1e-10)),
                              tie_abs=float(a.get("tie_absolute_tolerance", 1e-12)), fit_direction=cfg.get("gp", {}).get("fit_direction", "requested"))


def build_gp_kwargs(cfg):
    g = cfg.get("gp", {}); kw = dict(optimize_hyperparameters=bool(g.get("optimize_hyperparameters", True)))
    for k in ("length_scale", "white_noise", "ml2_min_points", "n_restarts", "random_state", "jitter"):
        if k in g:
            kw[k] = g[k]
    return kw


def build_oracle(cfg):
    o = cfg["oracle"]; t = o["type"]
    if t == "analytical":
        return CachedOracle(analytical_oracle(o["model"], **o.get("parameters", {})))
    if t == "external_command":
        from .solvers.external import ExternalCommandOracle
        base = Path(cfg["_config_path"]).parent
        cmd = [sys.executable if c == "{python}" else (str((base / c).resolve()) if ("{" not in c and (base / c).exists()) else c) for c in o["command"]]
        return CachedOracle(ExternalCommandOracle(cmd, o["work_dir"], o.get("timeout_s", 3600)))
    if t == "moose":
        from .solvers.moose import MooseYieldOracle
        kw = {k: v for k, v in o.items() if k not in ("type",)}
        return CachedOracle(MooseYieldOracle(**kw))
    raise ValueError(f"unknown oracle type {t!r}")


def build_controller(cfg, oracle=None):
    return AdaptiveController(build_library(cfg), oracle or build_oracle(cfg), build_symmetry(cfg), build_settings(cfg), build_gp_kwargs(cfg))


def load_eval_directions(path):
    return pd.read_csv(path, float_precision="round_trip")[["z1", "z2", "z3"]].to_numpy(float)


def environment_info():
    import numpy, scipy, sklearn, matplotlib
    return dict(python=platform.python_version(), platform=platform.platform(), numpy=numpy.__version__, scipy=scipy.__version__,
                scikit_learn=sklearn.__version__, pandas=pd.__version__, matplotlib=matplotlib.__version__)


def write_json(obj, path):
    Path(path).write_text(json.dumps(obj, indent=1, default=lambda x: x.tolist() if isinstance(x, np.ndarray) else float(x) if isinstance(x, np.floating) else str(x)))
