#!/usr/bin/env python3
"""POST-ANALYSIS DESCRIPTIVE (written after analyze_h1.py; not part of the frozen protocol; no claim): why does V3 beat maximin when its gain
is not attributable to the material's responses?  Geometry of the n = 33 designs and the GP hyper-parameters fitted to them.
    python post_analysis_geometry.py  ->  analysis/post_geometry.json"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd
HERE = Path(__file__).resolve().parent; sys.dont_write_bytecode = True; sys.path.insert(0, str(HERE.parents[1] / "complexity_scaling_2026"))
import cs_common as C                                                                # noqa: E402
LIB = C.library().directions; SPH = C.fibonacci_sphere(20000); SPH = SPH[SPH[:, 2] > 0]; M = pd.read_csv(HERE / "evaluation/metrics_all.csv"); out = {}
ax = lambda A, B: np.degrees(np.arccos(np.clip(np.abs(A @ B.T), 0, 1)))
EQB = np.array([[1, 1, 0]]) / np.sqrt(2); rows = []
for arm in ("V3", "V3_ISO", "V2", "AXMM", "RANDOM"):
    for i in range(24):
        X = LIB[json.load(open(HERE / f"runs/{arm}/r{i:02d}/design.json"))["axes"][:33]]; d = ax(X, X); np.fill_diagonal(d, 180)
        rows.append(dict(arm=arm, i=i, nn_median=float(np.median(d.min(1))), nn_min=float(d.min()), covering=float(ax(SPH, X).min(1).max()),
                         share_within_30deg_of_equal_biaxial=float((ax(X[8:], EQB)[:, 0] < 30).mean()), share_shear_dominant=float((np.abs(X[8:, 2]) > 0.7).mean())))
G = pd.DataFrame(rows); out["design_geometry_median"] = G.groupby("arm")[["nn_median", "nn_min", "covering", "share_within_30deg_of_equal_biaxial", "share_shear_dominant"]].median().round(3).to_dict("index")
m = M[M.n == 33]; out["gp_n33_median"] = m.groupby("arm")[["length_scale", "amplitude"]].median().round(4).to_dict("index")
out["gp_n33_share_ls_at_upper_bound"] = m.groupby("arm").length_scale.apply(lambda x: float((x > 4.999).mean())).to_dict()
out["note"] = "descriptive, post-analysis; area share within 30 deg of equal biaxial on the axis hemisphere = %.3f" % float((ax(SPH, EQB)[:, 0] < 30).mean())
json.dump(out, open(HERE / "analysis/post_geometry.json", "w"), indent=1); print(json.dumps(out, indent=1))
