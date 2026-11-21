# Reproducibility

## 1. Exact frozen algorithm

| Step | Specification |
|---|---|
| Representation | Plane stress, Mandel coordinates z = [σxx, σyy, √2 σxy]. A test is a unit direction d; the oracle returns the first-yield stress, whose Mandel norm is the yield radius r(d). |
| Initial design | 8 tests by greedy angular maximin over the candidate library, using geometry only (no oracle values). Start 0 begins at candidate 0; start s ≥ 1 begins at a random candidate (`default_rng(7000 + s)`). |
| Surface model | Direct Matérn-3/2 GP on y = log r: `ConstantKernel(amp) × Matern(ν = 1.5) + WhiteKernel`, constant mean = median(y), jitter 1e-8. Observations at the same stress point are merged. |
| Hyper-parameters | Below 20 training points (pseudo-observations count): amp = max(var(y), 1e-4), length scale 0.5 and white noise 1e-3, all fixed. From 20 points: ML-II over amp ∈ [1e-5, 10], length scale ∈ [0.05, 5], white noise ∈ [1e-6, 0.1]; 2 restarts, `random_state` 0. |
| Acquisition | Maximum predictive uncertainty: argmax over untested candidates of the GP predictive variance of log r, from the same GP. |
| Ties (deployment) | Tie set = {score ≥ max − (\|max\|·1e-10 + 1e-12)}; the lowest candidate index wins. |
| Ties (registered runs) | `np.argmax`, i.e. round-off decided ties. |
| Loop | One oracle call per selected test; refit; repeat to the budget of **real** tests. |
| Physics-informed option | Shear mirror M = diag(1, 1, −1): after every real test, a pseudo-observation at M d with noise 1e-4, added iff M d is more than 1e-6 rad from every real test. Pseudo-observations are never counted as tests. |
| Generic blind baseline | Identical, without the shear mirror. |
| Reconstruction | r̂(d) = exp(m + μ_GP(d)), with predictive sd in log space. The registered scoring used a pinned icosphere mesh; the package uses the continuous prediction (≤ 0.04 pp difference on the frozen runs). |

## 2. Exact environment

**Interpreter:** Python 3.12.9 (CPython), Linux x86_64.

**Pinned packages** (`software/yield_surface_discovery/requirements-lock.txt`): numpy 1.26.4, scipy 1.15.2, pandas 2.2.3, scikit-learn 1.6.1, matplotlib 3.10.1, pyyaml 6.0.2, pytest 8.3.5, plus their transitive pins.

**Other runtime requirements:** no GPU, no MOOSE, no mesh library.

**Threads:** scripts and tests set `OMP/OPENBLAS/MKL_NUM_THREADS=1` by default.

**Exact reproduction:** bit-exact reproduction of the registered round-off-tie runs requires these exact versions. The portability test (§8) confirms it in a fresh virtual environment built only from the lock file. With other versions, deterministic mode, the registered-choice replay and all metrics still agree; only legacy re-runs may diverge (pytest marks them xfail).

## 3. Reproduce the analytical benchmarks

From `software/yield_surface_discovery/`, after installing the environment:

```
python scripts/validate_install.py
python -m pytest -q                                                     # 32 tests (~10 s)
python scripts/run_benchmark.py --config configs/reproduce_registered_physics_informed.yaml   # registered Arm C, bit-exact
python scripts/run_benchmark.py --config configs/reproduce_registered_blind.yaml              # registered Arm A, bit-exact
python scripts/run_benchmark.py --config configs/benchmark_physics_informed.yaml              # deployment configuration
python scripts/run_benchmark.py --config configs/benchmark_blind.yaml
python scripts/run_minimal_example.py --check                                                 # 12-test example vs stored result
```

**What to expect from the reproduction configs:**
- They reproduce the registered sequences and radii exactly.
- Their reported metrics are continuous-prediction values. The registered *continuous* global p95 is reproduced to ≤ 2e-16; the registered tables report the pinned-mesh values, which are within 0.04 pp.
- For bitwise comparisons of the written CSV files, read them with `pandas.read_csv(..., float_precision="round_trip")`. The default parser can be off by one ulp; the JSON observation files are exact.

## 4. Reproduce the figures and tables

```
python scripts/reproduce_figures.py --output ./figures --tables ./tables
```

This reads `data/frozen_results/` (authoritative saved values extracted by `release/build/extract_frozen_results.py`, with source hashes in `data/frozen_results/PROVENANCE.json`). It writes the 12 figures (PDF + PNG + plotted values) and the 6 tables. In the portability test, all regenerated plotted-value CSVs and tables were byte-identical to `release/figures` and `release/tables`.

## 5. Configure a new candidate library

1. Create a CSV with columns `id, label, family, z1, z2, z3`: one row per physically executable loading test, giving the unit direction of its stress state in Mandel coordinates. Use `s_xx_dir, s_yy_dir, s_xy_dir` if you prefer Cauchy components.
2. Point `candidate_library.path` to it. Set `deduplicate: true` to merge exact duplicate directions (tolerance 1e-9).
3. Set `symmetry.shear_mirror` to `true` only if the material is known to be orthotropic in the test axes.

## 6. Connect an external solver

**Interface:** subclass `YieldOracle` and implement `evaluate(direction) -> YieldObservation`. A `YieldObservation` holds:
- requested direction;
- achieved direction;
- stress at yield;
- yield radius;
- success / status;
- metadata.

**Command-line codes:** use `oracle: {type: external_command, command: [...], work_dir: ...}`.
- Placeholders: `{z1} {z2} {z3} {s_xx} {s_yy} {s_xy} {run_dir} {run_id} {python}`.
- The command writes `result.json` (see `examples/external_solver_stub/`).

**MOOSE:** use `oracle: {type: moose, executable, input_template, work_dir, material_parameters, ...}` (see `examples/moose_adapter/README.md`).
1. Run `python scripts/check_solver_interface.py --config <your config> [--reference von_mises --sigma-y 250]` first. It runs the five deterministic loading families.
2. Then run `python scripts/run_closed_loop.py --config <your config>`.

**What the controller receives:** only the observation fields. Material parameters stay on the oracle side.

**Failures:** a failed test (`success = False`) marks that candidate unavailable, and the next candidate is used.

## 7. Registered primary vs exploratory results

| Status | Results |
|---|---|
| **Registered, primary** | Same-surrogate Arms A–D (`protocol/same_surrogate_ablation_protocol.json`); the blind GP study's primary arm and decomposition (`protocol/blind_gp_mes_protocol.json`); the symmetry-aware cleanup (`protocol/symmetry_aware_cleanup_protocol.json`; not adopted); the expanded-library pool experiment and the initialisation experiment (template-informed era) |
| **Frozen existing** | Same-points Fuhg comparison (`frozen_protocol.json`) |
| **Exploratory / post hoc** | Post-hoc decomposition arms; Tsai-Wu direct GP + mirror; reachability split (README §7); arbitrary stress-direction experiment; edge follow-up; controller tolerance check (analytical, not MOOSE) |

## 8. Portability test

**Script:** `release/build/portability_test.sh`; log in `release/build/portability_test_log.txt`.

**Procedure:**
1. Copy only `software/yield_surface_discovery/` to a fresh temporary directory.
2. Scan it for workstation paths and symlinks (none found).
3. Build a new virtual environment from `requirements-lock.txt` and `pip install .`.
4. Run:
   - pytest (32 passed, including the bit-exact registered re-runs);
   - `validate_install.py` (all checks passed);
   - the minimal example from outside the package (PASS against the expected result);
   - both benchmarks (one start);
   - the external-solver loop;
   - figure and table reproduction.

## 9. Known limitations

- **Targets:** primary results are for one target (Mg–4%Li), 5 starts (4 distinct designs for the frozen method), one candidate library per comparison.
- **Pools 1 and 3:** the frozen practical method is registered on the expanded library only. Pool 1/3 comparisons use the same reconstructor and symmetry on the frozen symmetry-assisted test sequences.
- **Response-awareness:** selections respond to measured values only through GP hyper-parameter learning (≥ 20 training points).
- **Tsai-Wu:** needle-shaped surfaces are poorly represented by the smooth stationary GP.
- **MOOSE:** the adapter and input template were never run with MOOSE. The first MOOSE run must verify object names and the yield-detection output.
- **Mesh:** the package mesh export is a continuous radial mesh, not the pinned mesh of the registered scoring.
