import json
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]; DATA = ROOT / "data"


def library(name="expanded_physical_loading_library", dedup=False):
    from yield_surface_discovery import CandidateLibrary
    return CandidateLibrary.from_csv(DATA / "candidates" / f"{name}.csv", deduplicate=dedup)


def registered():
    return json.loads((DATA / "reproduction" / "registered_arms_A_C.json").read_text())["runs"]


def eval_set():
    from yield_surface_discovery.io import load_eval_directions
    from yield_surface_discovery.analytical import cb2004_radius
    from yield_surface_discovery.metrics import region_masks
    D = load_eval_directions(DATA / "benchmark" / "heldout_eval_rays.csv"); r = cb2004_radius(D); return D, r, region_masks(r, D)


def pinned_environment():
    import numpy, scipy, sklearn
    lock = {l.split("==")[0].lower(): l.split("==")[1].strip() for l in (ROOT / "requirements-lock.txt").read_text().splitlines() if "==" in l and not l.startswith("#")}
    return numpy.__version__ == lock["numpy"] and scipy.__version__ == lock["scipy"] and sklearn.__version__ == lock["scikit-learn"]
