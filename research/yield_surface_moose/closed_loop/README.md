# Closed-loop BO ↔ MOOSE (plane-stress composite benchmark)

The existing BO package, **unmodified**, drives the psbench MOOSE oracle.

* **BO package:** `~/projects/yield_surface_discovery_release/software/yield_surface_discovery` (v1.0.0).
  `~/Desktop/archives/BO_benchmark` does not exist. This is the only Matérn-GP / angular-maximin code on the machine. Its BO source
  files match the release `MANIFEST.csv` sha256 hashes, and the runner checks this before every run.
* **Python env:** `~/micromamba/envs/ysd` (Python 3.12.9 with the package's pinned `requirements-lock.txt`). MOOSE trials run in the `moose` env.

```
python run_bo_moose.py --config configs/S2_blind.yaml --out results/S2_blind_run1           # real MOOSE (~7 min)
python run_bo_moose.py --config configs/S2_blind.yaml --out results/S2_blind_replay \
       --replay results/S2_blind_run1/trials.jsonl                                           # BO re-run on recorded observations
python analyze_run.py results/S2_blind_run1 --replay results/S2_blind_replay                # checks, reference comparison, figures
```
(use `~/micromamba/envs/ysd/bin/python`)

## How the pieces connect

| BO package expects | psbench provides | Mapping |
|---|---|---|
| `YieldOracle.evaluate(d)`, d = unit Mandel 3-vector [σxx, σyy, √2 τxy] | `observe(direction, specimen_id)`, same Mandel convention, normalised internally | direction passed unchanged |
| `yield_radius` r > 0 (the GP takes log r itself) | `yield_radius` [MPa] (and `log_yield_radius`) | `yield_radius`. `log_yield_radius` is recorded but not used, because the GP computes log r itself |
| `achieved_direction`, `stress_at_yield` (Cauchy: s_xy tensor shear) | `stress_mandel` | achieved = stress_mandel / ‖·‖; Cauchy = stress_mandel · [1, 1, 1/√2] |
| `success` | `accepted` (converged + direction and plane-stress checks) | success = accepted and r finite > 0. Otherwise the package marks the candidate unavailable and never fits it |

* The BO process never imports `psbench` (asserted in the runner). Each trial runs in a separate `moose_worker.py` process.
  That process receives only (specimen_id, direction, tag) and returns only `observe()`'s public keys.
* Units: MPa in, MPa out. The GP is scale-free in r because it centres log r on the median.

## Files

`run_bo_moose.py` (runner + `PsbenchMooseOracle`, `ReplayOracle`), `moose_worker.py`, `analyze_run.py`, `configs/S2_blind.yaml`.
`configs/homog_smoke.yaml` is a plumbing check on the homogeneous material.
Run outputs are in `results/<run>/`, and MOOSE work dirs in `../runs/S2/h40/S2_blind_run1_tNN/`.
