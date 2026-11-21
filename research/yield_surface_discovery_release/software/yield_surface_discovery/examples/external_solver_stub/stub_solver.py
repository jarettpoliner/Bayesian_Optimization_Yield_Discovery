#!/usr/bin/env python3
"""Stand-in for an external mechanics solver (replace with your own code).  Usage:
    python stub_solver.py --z1 Z1 --z2 Z2 --z3 Z3 --out result.json
It 'runs' one proportional loading test along the unit Mandel direction (z1, z2, z3) and writes the first-yield stress as JSON:
    {"stress_at_yield": [s_xx, s_yy, s_xy], "status": "ok", "metadata": {...}}
The hidden material here is an orthotropic quadratic (Hill-type) criterion with yield stress 250 MPa; the adaptive controller never sees it."""
import argparse, json, math
ap = argparse.ArgumentParser(); [ap.add_argument(f"--{k}", type=float, required=True) for k in ("z1", "z2", "z3")]; ap.add_argument("--out", required=True)
a = ap.parse_args(); sxx, syy, sxy = a.z1, a.z2, a.z3 / math.sqrt(2.0)          # unit direction in Cauchy components
F, G, H, N, sy = 0.45, 0.55, 0.50, 1.40, 250.0                                    # hidden Hill-48 coefficients (oracle side only)
q = G * sxx ** 2 + F * syy ** 2 + H * (sxx - syy) ** 2 + 2 * N * sxy ** 2
lam = sy / math.sqrt(q)                                                            # load factor at first yield
json.dump(dict(stress_at_yield=[lam * sxx, lam * syy, lam * sxy], status="ok", metadata=dict(solver="stub_solver.py", load_factor=lam)), open(a.out, "w"))
