# Release report — yield_surface_discovery 1.0.0

This release closes out the analytical benchmark and packages the frozen method. No scientific experiment was re-run or changed; every number below comes from saved registered results (sources in `AUDIT.md` and `software/yield_surface_discovery/data/frozen_results/PROVENANCE.json`).

## 1. What exact algorithm was frozen?

1. 8 angular-maximin initial tests from the candidate library (geometry only).
2. Direct Matérn-3/2 GP on log yield radius (constant median mean; frozen hyper-parameter schedule: fixed below 20 training points, ML-II from 20).
3. Maximum predictive-uncertainty sampling over the untested candidates, using the same GP.
4. One real test per selection; refit; repeat.
5. **Physics-informed (practical) version:** add known shear-mirror pseudo-observations, and use the expanded physical loading library.

This is registered same-surrogate Arm C, confirmed as the frozen method by the registered symmetry-aware cleanup.

**Deployment details** (bookkeeping, recommended by the cleanup's keep outcome): deduplicated candidate list and deterministic ties.

## 2. Which analytical templates were removed?

All of them. The surface model is a GP on log r with a constant mean: no Hill-1948, Drucker-Prager, Gurson, Hosford, Cazacu-Barlat or any other criterion, and no residual relative to one.

At identical tests and symmetry, adding the calibrated Hill-1948 template changes median global p95 by +10.9 / +0.9 / +0.0 pp at 12 / 16 / 28 tests, i.e. never better. The template helps only when the target lies in or near its family (post-hoc cross-surface evidence).

## 3. What role does shear-mirror symmetry play?

It is the largest single prior. It is exact for orthotropic materials tested in their axes, and it doubles the information per real test through free pseudo-observations. Removing it raises median global p95 by:

| Tests | 12 | 16 | 20 | 28 |
|---|---|---|---|---|
| Δ median global p95 (pp) | +19.1 | +7.7 | +1.4 | +2.4 |

All 5 starts are worse without it from 12 tests on. Registered decision: Case 4, symmetry essential.

Because it is material-specific knowledge, the method that uses it is **physics-informed, not blind**.

## 4. Why predictive uncertainty instead of IVR?

**Registered result:** IVR (integrated variance reduction from the same GP) and predictive uncertainty tie by the pre-specified rule at both symmetry levels. IVR minus predictive uncertainty, with symmetry: +0.0 pp at 16 tests and +0.6 pp at 28 (median difference). Predictive uncertainty is simpler, so it is kept (Case 2).

**Neither criterion is response-aware below 20 training points:** a GP's posterior variance does not depend on the observed values at fixed hyper-parameters (registered diagnostic, confirmed).

## 5. What is the generic blind baseline?

The same algorithm without the shear mirror, and without any analytical template (registered Arm A).

On the expanded physical library, median global p95 is 29.8 / 26.7 / 11.5 / 4.9 / 4.3 / 4.0% at 8 / 12 / 16 / 20 / 24 / 28 tests.

## 6. What is the practical physics-informed method?

The frozen algorithm plus the shear mirror on the expanded physical library (Arm C). Median global p95 (5 starts, 4 distinct designs):

| Tests | 8 | 12 | 16 | 20 | 24 | 28 |
|---|---|---|---|---|---|---|
| Median global p95 (%) | 30.1 | 7.6 | 3.7 | 3.5 | 3.2 | 1.6 |

At 28 tests: weak-region maximum overprediction 0.11%, sharp-region 2.13%, high-shear p95 0.62%.

## 7. What does Pool 2 contain?

The **expanded physical loading library**: 416 labels (414 unique directions; 2 exact duplicate labels).

| Family | Count |
|---|---|
| Uniaxial tension/compression (5° steps) | 72 |
| Aligned-cruciform biaxial | 68 |
| In-plane pure shear (5° orientations, both signs) | 36 |
| Combined normal + shear proportional loading, axis-aligned (butterfly/Arcan type, 15° grid) | 40 |
| Combined normal + shear proportional loading, rotated | 200 |

The conventional library (Pool 1) is the first two families (140 tests). Its maximum shear content is |z3|/|z| = 1/√2.

## 8. Headline benchmark results

**Loading space matters most at higher budgets.** Same direct GP + mirror reconstruction on the frozen symmetric test sequences, median global p95 at 28 tests:

| Library | Global p95 at 28 tests |
|---|---|
| Conventional | 7.1% |
| Expanded physical | 1.2% |
| Arbitrary-direction reference | 1.2% (not an upper bound) |

The conventional floor sits in the high-shear directions that uniaxial/biaxial tests cannot reach (Figure 2).

**Ranking of ingredients** (Figure 4; Δ median global p95 at 16 / 28 tests when one ingredient changes):

| Change | 16 tests | 28 tests |
|---|---|---|
| Conventional instead of expanded library | +3.9 pp | +5.8 pp |
| No symmetry | +7.7 pp | +2.4 pp |
| IVR instead of predictive uncertainty | +0.0 pp | +0.6 pp |
| Hill template added | +0.9 pp | +0.0 pp |

**Physics-informed vs generic blind at 28 tests:** 1.6% vs 4.0% global p95.

## 9. What comparison to Fuhg et al. (2023) is defensible?

Only the **information-matched reconstruction comparison on identical published points** (D8 / D16 / D28, no symmetry for anyone, single runs):

| Tests | Result |
|---|---|
| 8 | Comparable (direct GP 21.5% vs Fuhg ICNN 21.0% / GP 22.4%) |
| 16 | Fuhg ICNN better (6.1 vs 9.4%) |
| 28 | Direct GP lowest global p95 (6.90 vs 8.6–9.0%) and sharp-region overprediction; weak-region overprediction higher than the Fuhg GP/SVR (2.5 vs 1.0 / 0.8%) |

**Not defensible:**
- general superiority over Fuhg;
- treating the adaptive Pool-2 results as a model comparison against Fuhg (they combine model, prior, test design and action space).

## 10. Major limitations

- **Scope:** one target (Mg–4%Li); 5 starts (4 distinct designs for the frozen method).
- **Pools 1 and 3:** the frozen practical method is registered on the expanded library only. The Pool 1/3 comparisons use the same reconstructor and symmetry on frozen symmetry-assisted sequences.
- **Response-awareness:** selection reacts to measurements only through hyper-parameter learning.
- **Symmetry must be known:** the shear mirror requires orthotropy in the test axes.
- **Surface shape:** the smooth stationary GP struggles with sharp or needle-like surfaces (Tsai-Wu, Tresca corners).
- **Plane stress only:** stress directions are 3-D Mandel vectors.
- **Mesh export:** the package's surface mesh is continuous, not pinned.
- **MOOSE:** the adapter is untested against real MOOSE.

## 11. Is Tsai-Wu still a limitation?

**Yes.** It is outside the primary benchmark, but it is unresolved. On the needle-shaped Tsai-Wu surfaces (r_min/r_max ≈ 0.01–0.03), the blind GP's global p95 at 28 tests is 51–69% on the expanded library.

At identical points the direct GP + mirror equals the template-informed bank (5.6 vs 5.6% at 100 tests), so the template is not the cause. The limitation is the representation: the narrow tip is underpredicted and the flanks overshoot (Figure S6).

## 12. Is the package reproducible from a fresh environment?

**Yes**, verified by `build/portability_test.sh`. Only the package was copied to a new temporary directory: no workstation paths or symlinks. A new virtual environment was built from `requirements-lock.txt` and the package installed. Then:

| Check | Result |
|---|---|
| Unit and reproduction tests | 32/32 passed, including bit-exact re-runs of all 10 registered Arm A/C sequences |
| `validate_install.py` | all checks passed; GP identical to the registered implementation |
| Minimal example | matched the stored result |
| Both benchmarks and the external-solver loop | ran |
| Regenerated figure data and tables | byte-identical to the release copies |

With other library versions, deterministic mode and metrics agree; only the legacy round-off-tie re-runs may diverge.

## 13. Can it run without MOOSE?

**Yes.** MOOSE is neither imported nor a dependency. The analytical oracle, the external-command oracle and custom `YieldOracle` subclasses need only Python. The Dockerfile is Python-only.

## 14. Is there a clean external-solver interface?

**Yes.** The controller calls only `YieldOracle.evaluate(direction) → YieldObservation`: requested and achieved direction, stress at yield, radius, success / status, metadata.

**Adapters:**
- `ExternalCommandOracle` for any command-line code (tested with a stub solver);
- `MooseYieldOracle`: stress-controlled single element, eqps > 1e-8 onset, ramp + bisection to 1e-4, direction and plane-stress checks, stdout / stderr / CSV kept per run. Tested against a fake MOOSE executable only.

Failed tests mark candidates unavailable without breaking the loop. `scripts/check_solver_interface.py` runs the five canonical loading families through any configured oracle.

## 15. What exact steps are required on another computer with MOOSE installed?

1. Copy and extract `yield_surface_discovery_release.tar.gz`. In `software/yield_surface_discovery/`, run `python3.12 -m venv .venv && . .venv/bin/activate && pip install -r requirements-lock.txt`.
2. Run `python scripts/validate_install.py && python -m pytest -q`. All must pass.
3. Have a MOOSE application with the Solid Mechanics module (source build or the official `idaholab/moose` container). Note its executable (or container wrapper → `launcher`).
4. Check `examples/moose_adapter/single_element_yield.i` against that revision: `<app> --registry`, and run `<app> --check-input` on one hand-filled copy. Fix object names if needed. This is the first step that has never been executed.
5. Edit `configs/example_moose.yaml`:
   - `executable` / `launcher`;
   - `material_parameters` (J2 first: e.g. E = 70 GPa, ν = 0.3, σy = 250 MPa, H = 700 MPa);
   - `work_dir`.
6. Run `python scripts/check_solver_interface.py --config configs/example_moose.yaml --reference von_mises --sigma-y 250 --tol 1e-3`. All five families must pass. Inspect `work_dir/test_*/` (input, CSV, stdout, stderr).
7. Register the MOOSE validation protocol before scoring. Proposed acceptance: with J2, the controller should reproduce the analytical von Mises behaviour within detection tolerance. Then run `python scripts/run_closed_loop.py --config configs/example_moose.yaml` with the generic blind and physics-informed settings.

## 16. Is the project ready for solver-in-the-loop validation?

**The software side is ready:**
- the frozen controller is solver-independent, verified and portable;
- the MOOSE adapter, input template, interface check and configuration exist;
- MOOSE is not required for anything else.

**Not yet demonstrated:**
- No real MOOSE run has occurred on any machine.
- The input template and the eqps output path are unverified.

So the next step is step 4 above on a machine with MOOSE. That may require template corrections, but no change to the adaptive algorithm.

## Release status: **READY FOR PORTABLE ANALYTICAL RELEASE**

The analytical release is complete, reproducible from a fresh environment, and verified bit-exact against the registered results.

The MOOSE path is packaged but unexecuted, so "READY FOR MOOSE VALIDATION" is not claimed. It will be the status once the template passes `check_solver_interface.py` on a machine with MOOSE.
