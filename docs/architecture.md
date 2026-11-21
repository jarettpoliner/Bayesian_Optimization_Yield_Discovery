# Architecture

The code is the authority for everything below. Abbreviations:
- `CL` = `research/yield_surface_moose/plane_stress/closed_loop`;
- `V1` = `research/yield_surface_discovery_release/software/yield_surface_discovery/src/yield_surface_discovery`.

## Closed loop

```
 candidate library ──► controller.run(start)
                          │ initial design: 8 directions, angular/axial maximin (response-independent)
                          ▼
                  ┌─► GPv2.fit(directions, radii)      Matérn-3/2 on log r, median mean, A20 schedule,
                  │        │                           antipodal pseudo-observations where validated
                  │        ▼ predict_log(candidates) → μ, latent σ²
                  │   score: V2 σ² | V3 σ²·ρ(LOO risk) | V3-Signed (experimental); tested/failed → −∞
                  │        ▼ select(): deterministic argmax (rel. tie tolerance 1e-10, lowest index)
                  │   oracle.evaluate(d̂) → YieldObservation(yield_radius, success)
                  └── append (failure → candidate unavailable) … until budget
 evaluator (separate): independent references → |r̂/r − 1|, paired statistics
```

## Components

| Component | Code | Interface |
|---|---|---|
| Matérn-3/2 GP (frozen) | `CL/v2/bomes_v2.py::GPv2` | `fit(D, r)`, `predict_log(D)` → (mean log r, latent variance), `hyperparameters()` |
| V2 controller, maximin | `CL/v2/bomes_v2.py::V2Controller`, `V2Settings`, `axial_maximin`, `AxisLibrary` | `run(start)` → sequence, observations, steps, scores |
| V3 controller | `CL/response_adaptive_v3_2026/v3_lrv.py::V3Controller`, `loo_standardised`, `risk` | same |
| V3-Signed (experimental) | `CL/truss_lattice_benchmark_2026/final_comparison_2026/v3transfer_20261120/v3s_lrv.py` | same |
| Candidate library, symmetry, tie rule, signed maximin, oracle interface | `V1/candidates.py`, `symmetry.py`, `acquisition.py::select`, `initialization.py::angular_maximin`, `oracle.py` | |
| Analytical surfaces | `CL/complexity_scaling_2026/cs_common.py` (A1–A6, 4,096-axis library), `CL/response_adaptive_v3_2026/holdout_h1_2026/surfaces_h1.py` (Yld91) | `radius(D)` |
| Truss MOOSE oracle | `CL/truss_lattice_benchmark_2026/oracle/truss_oracle.py`, `controller/truss_v2_compat.py` | writes MOOSE input, runs, extracts r at ‖E_p‖ = 1e-4 |
| Composite MOOSE oracle | `research/yield_surface_moose/plane_stress/psbench/`, `CL/biphasic_benchmark/oracle/`, `CL/biphasic_benchmark/final_bench_20261119/bench_oracle.py`, `CL/run_bo_moose.py`, `CL/moose_worker.py` | plane-stress periodic RVE, r at RVE-average eqps = 1e-4 |
| Public API/CLI | `software/yield_surface_bo/`: `core.py` (methods, libraries, compatibility checks, analytic/replay/checkpoint/external-command oracles), `campaign.py` (`run_campaign`, `fit_gp`, `predict`, `evaluate_errors`), `moose_truss.py` (live TL1), `_frozen.py` (paths + SHA-256 of 43 frozen files), `cli.py` | `./ysbo`, `from yield_surface_bo import run_campaign` |
| Evaluation scripts | `CL/biphasic_benchmark/final_bench_20261119/evaluate_final.py`, `verify_final.py`; `…/holdout_h1_2026/evaluate_h1.py`, `analyze_h1.py`; `CL/complexity_scaling_2026/evaluate_analytic_r3.py`; TL1 `evaluate_c3.py`, `evaluate_v3x.py` | |

## Dependencies

- **Python ≥ 3.10:** numpy, scipy, pandas, scikit-learn, matplotlib, pyyaml. The pins are in `requirements-lock.txt`.
- **MOOSE:** only for live physical runs.

## Data flow of a recorded physical observation

d̂ → MOOSE input decks (`probe.i`, `ramp.i`, `refine_*.i`) → MOOSE → CSV histories → yield extraction → `record.json` → controller.

A complete example trial is included: `CL/biphasic_benchmark/runs/B1/h40/FB_B1_pool_t01/` (composite) and `CL/truss_lattice_benchmark_2026/validation/runs/V2_20261120T170634/obs_00/` (truss).
