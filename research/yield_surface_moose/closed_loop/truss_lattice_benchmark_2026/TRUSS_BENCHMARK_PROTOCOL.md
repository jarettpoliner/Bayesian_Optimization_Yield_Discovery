# TL1 matched comparison: proposed protocol (DRAFT, not frozen, not executed)

Status: proposal for a future run. **It is not authorized.** Nothing below has been run. Before any campaign, the user must authorize the budget, and this file must be frozen and hashed.

## Fixed inputs (already frozen and hashed)

| Item | Artefact |
|---|---|
| RVE | TL1 (`lattice/TL1_frozen.json`, sha256 37cec8d4…) |
| Oracle | `oracle/truss_oracle.py` as validated in V2, with `oracle/calibration_TL1.json`. Criterion: Levin macroscopic plastic strain of 1e-4. Solver settings as in `DEFAULTS` |
| Candidates | TLIB-4096 v1: **8,192 signed directions**, `symmetry='none'` (antipodal not admissible) |
| References | TREF v1: 482 MOOSE references, already measured (`pilot/refs_TREF_v1_20261120T171229/`), evaluator-only |
| Controller | Frozen `bomes_v2.V2Controller` through TV2C 1.0, `V2Settings(symmetry='none', schedule='A20')`. Matérn-3/2 kernel, maximum latent predictive variance, ML-II from n = 20 at every step |

## Methods and matching (per start k)

- **A. BO+MES:** the frozen V2-A20 controller.
- **B. Angular maximin:** frozen V1 deterministic rule on the 8,192 signed directions, with the same first direction as BO start k.
- **C. Random:** unique directions drawn by `default_rng(9000 + k)`, after the shared initialization.

All three methods share the identical 8-point initialization of start k, the angular-maximin initial design that BO itself uses. Those 8 points are reused exactly.

**Starts.** k = 0..8 for budgets up to n = 64, about 3.4 CPU-h. Starts {0, 1, 2} only for the n = 1,000 extension, about 21 CPU-h.

**Bookkeeping.**
- A failed MOOSE trial is not counted. BO marks the direction unavailable; the baselines take their next direction.
- Exact (TL1, oracle version, direction) reuse is allowed within the benchmark only, never from the references.
- Every request, whether new, reused or failed, is logged.

## Checkpoints and endpoints

**Checkpoints:** n = 8, 16, 33, 64, 128, 256, 512, 1000.

**Primary endpoint:** mean relative yield-radius error over the 338 area-uniform references at **n = 33**. The set is area-uniform, so the mean is area-weighted. Surrogate: the frozen GPv2 (A20, no symmetry), fitted identically to each design.

**Secondary endpoints:**
- maximum and P95 error;
- error on the two sweeps;
- error near cusps and edges: references within 10° of a mechanism switch identified from the reference data;
- 90% and 95% coverage;
- error against cumulative CPU cost;
- hyperparameter trajectories.

The **long-budget endpoint** is n = 1,000.

**Sensitivity:** the FIXED surrogate (ℓ = 0.5).

## Claim rule (proposed)

BO is claimed better than maximin at n = 33 only if all of the following hold:
1. it wins in ≥ 8 of 9 paired starts;
2. the paired median improvement exceeds 0.2 percentage points, which is roughly 10× the oracle's numerical uncertainty and the tension/compression asymmetry scale;
3. conditions 1 and 2 also hold under the FIXED surrogate.

Wilcoxon signed-rank results are reported. At n = 1,000 (3 starts) only a case-study description is allowed.

**Response dependence:** report the first test at which BO's sequence diverges from a run on a response-free surrogate (constant r). This separates kernel-aware geometric design from response adaptivity.
