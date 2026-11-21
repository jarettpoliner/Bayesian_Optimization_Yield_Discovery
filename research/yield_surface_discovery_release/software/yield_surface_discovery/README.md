# yield_surface_discovery

Given a library of possible **physical loading tests**, this software chooses a short sequence of experiments to reconstruct an **unknown plane-stress yield surface**. Each test is run by an *oracle*:

* a laboratory test;
* MOOSE or another mechanics solver;
* an analytical model, for benchmarking.

The adaptive algorithm never sees the material model: it sees only the requested loading directions and the measured first-yield points.

## Frozen practical workflow

```
8 angular-maximin initial tests (geometry only)
→ direct Matérn-3/2 Gaussian process on log yield radius (constant mean; no analytical yield criterion)
→ maximum predictive-uncertainty sampling over the untested candidates, using the same GP
→ run that test, update the GP, repeat until the budget of real tests is used
```

**Optional physics-informed enhancement:** the known shear-mirror symmetry (σxy → −σxy). It is exact for orthotropic materials tested in their symmetry axes, and it adds a free pseudo-observation at the mirror image of every real test.

| Mode | Analytical yield criterion | Symmetry | Config |
|---|---|---|---|
| **Generic blind** | none | none | `symmetry: {shear_mirror: false}` |
| **Physics-informed** (the practical method) | none | known shear mirror | `symmetry: {shear_mirror: true}` |

The physics-informed mode is **not fully blind**: it assumes the shear-mirror symmetry. "Number of tests" always means **real experiments (oracle calls)**; symmetry pseudo-observations are never counted.

## Inputs and outputs

**Inputs:**
* **Candidate loading directions:** a CSV of unit plane-stress directions in Mandel coordinates z = [σxx, σyy, √2 σxy] (columns `id, label, family, z1, z2, z3`).
* **An oracle** that can execute one proportional loading test along a direction and report first yield.

**Outputs:**
* the selected experiment sequence (`tests.csv`);
* every yield observation: requested and achieved direction, stress at yield, radius, status, metadata (`observations_start*.json`);
* the GP reconstruction and its uncertainty (`section_sxy0_start*.csv`), the reconstructed surface mesh (`surface_start*.obj`) and GP hyper-parameters (`steps.csv`);
* evaluation metrics when reference data are available (`per_run.csv`, `summary.csv`).

**Bundled candidate libraries** (`data/candidates/`):

| File | Contents |
|---|---|
| `conventional_loading_library.csv` | 140 uniaxial + aligned-biaxial tests |
| `expanded_physical_loading_library.csv` | 416 labels (414 unique directions): the conventional tests + in-plane pure shear + combined normal–shear tests. The recommended library. |
| `arbitrary_direction_reference.csv` | 1,504 directions; a relaxed reference, **not** a physical test plan |

## Install

Python 3.12 is recommended; ≥ 3.10 works. The pinned versions reproduce the registered results bit-for-bit.

```
python3.12 -m venv .venv
. .venv/bin/activate              # Windows: .venv\Scripts\activate
pip install -r requirements-lock.txt
pip install -e .                  # optional; every script also works from the source tree
```

Alternatively: `conda env create -f environment.yml`, or `docker build -t ysd .` (Python only; no MOOSE inside).

## Quick start

```
python scripts/validate_install.py                                   # checks environment, data, GP equivalence, controller, adapters
python -m pytest -q                                                  # unit + reproduction tests (~10 s)
python scripts/run_minimal_example.py --check                        # 12-test analytical example + figure, compared with the stored result
python scripts/run_benchmark.py --config configs/benchmark_physics_informed.yaml
python scripts/run_benchmark.py --config configs/benchmark_blind.yaml
python scripts/reproduce_figures.py --output ./figures --tables ./tables   # release figures/tables from the bundled frozen results
python scripts/run_closed_loop.py --config configs/example_solver.yaml     # closed loop with an external command-line solver (stub)
python scripts/check_solver_interface.py --config configs/example_solver.yaml   # 5 deterministic loading families through any oracle (do this first on a new solver)
```

Outputs go to the `output.directory` of the config (default `outputs/`) or to `--output`.

## Connecting a solver

The controller calls only `YieldOracle.evaluate(direction) -> YieldObservation`. There are three ways to connect a solver:

1. **Any command-line code** (`oracle.type: external_command`): the command receives the direction through placeholders (`{z1} {z2} {z3} {s_xx} {s_yy} {s_xy} {run_dir}`) and writes `result.json` (see `examples/external_solver_stub/`).
2. **MOOSE** (`oracle.type: moose`; optional, never a dependency): a stress-controlled single-element test. A ramp plus bisection finds the first increase of equivalent plastic strain (see `examples/moose_adapter/` and `configs/example_moose.yaml`). The adapter is tested only against a fake executable; **it has not been run with MOOSE.**
3. **Python:** subclass `yield_surface_discovery.YieldOracle` (e.g. for laboratory hardware) and pass it to `AdaptiveController`.

Failed tests (`success = False`) mark the candidate unavailable, and the controller picks the next one.

## Configuration

Everything is set in YAML (`configs/`). Relative paths are resolved against the config file.

| Section | Contents |
|---|---|
| `candidate_library` | path, deduplicate |
| `symmetry` | shear_mirror |
| `gp` | frozen: matern32 on log_radius; hyper-parameter optimisation from 20 training points |
| `initialization` | angular_maximin, n_initial, start rule, ties |
| `acquisition` | predictive_variance, tie tolerance |
| `budget` | total_real_tests, report_at |
| `starts` | — |
| `oracle` | — |
| `evaluation` | optional reference for scoring |
| `output` | — |

The configs contain no target-specific template information. The `oracle` section describes the hidden "material", which only the oracle sees.

**Reproducing the registered runs exactly:** `configs/reproduce_registered_{physics_informed,blind}.yaml` use the original 416-label library and the registered legacy round-off tie handling. They reproduce the registered sequences bit-for-bit **only** with the pinned environment. The default (deterministic) configs are robust to small numerical differences.

## Numerical determinism

Scripts default to one BLAS thread (`OMP/OPENBLAS/MKL_NUM_THREADS=1`). In deterministic mode, ties among acquisition scores (relative tolerance 1e-10) go to the lowest candidate index.

## Layout

| Path | Contents |
|---|---|
| `src/yield_surface_discovery/` | `geometry`, `candidates`, `symmetry`, `gp`, `initialization`, `acquisition`, `oracle`, `analytical`, `controller`, `reconstruction`, `metrics`, `io`, `workflows`, `figures`, `tables`, `solvers/{external,moose}` |
| `scripts/` | `validate_install`, `run_benchmark`, `run_closed_loop`, `run_minimal_example`, `reproduce_figures` |
| `configs/` | blind / physics-informed benchmarks; registered reproduction; external solver; MOOSE |
| `data/` | candidate libraries, held-out evaluation directions, frozen results (for figures/tables), reproduction references, minimal example; checksums in `data/SHA256SUMS` |
| `tests/` | pytest suite, including registered-result reproduction and a fake MOOSE executable |
| `examples/` | analytical oracle, external-solver stub, MOOSE adapter template |

See `../../REPRODUCIBILITY.md` and `../../RELEASE_REPORT.md` in the release bundle for the scientific provenance.
