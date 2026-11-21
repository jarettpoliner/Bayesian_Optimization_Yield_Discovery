# Final B1/B2 BO benchmark: frozen protocol (2026-11-19, written and hashed before any benchmark MOOSE run)

## Materials, oracle and observation

| Item | Specification |
|---|---|
| Materials | B1 and B2 from `oracle/materials_frozen.json` (sha256 499cde8e…), unchanged. Controlled synthetic benchmarks |
| Oracle | `oracle/biphasic_oracle_r1.py` and writer r1, through `oracle/biphasic_worker_r2.py` |
| Run settings | Corners pin, S2 mesh h40, dr_rel 0.05, n_refine 5, PJFNK/LU, nl 1e-9: the S1/S1b-validated configuration |
| Observation (primary) | Whole-cell volume-average eqps = 1e-4, giving the yield radius r. The BO models log r |
| Secondary diagnostics (simulator side, `record.json`) | Matrix- and inclusion-average eqps; phase shares f_p ε̄_p / ε̄_RVE; first local yield per phase; plastic volume fractions; convergence and cutbacks |

## Controller and baselines

| Item | Specification |
|---|---|
| BO | Frozen `v2/bomes_v2.V2Controller` with `V2Settings(schedule="A20", symmetry="antipodal", budget=33)`: 8 axial-maximin initial tests, then 25 maximum-predictive-variance selections. ML-II from 20 real observations; deterministic ties |
| Symmetry | Antipodal: 207 axes and pseudo-observations (validated in S1b-min). Shear mirror: off |
| Candidate library | 414 V1 directions mapped to 207 axes, with lowest-index representatives (`designs_common.py`) |
| BO and axial-maximin starts | k ∈ {0, 1, 2, 3, 5, 6, 7, 8, 9}: the 9 distinct initial designs. Axial maximin at start k uses the same first axis and is response-independent |
| Random | Seeds 1–10: `default_rng(9000 + k).choice(207, 45)`, the first 33 successful axes |
| Budget | 33 successful physical tests per campaign. Pseudo-observations are not counted. Failures are recorded and replaced (BO: axis marked unavailable; baselines: next spare) |
| Reuse | Exact (variant, direction) reuse only, within this benchmark. Never across variants, and never from the references |

## References (evaluator side only; `refs/`; never read by a controller)

- **Design:** frozen `refs/reference_design.csv` (sha256 in `.sha256`).
  - 60 area-stratified axes, each ≥ 3° from every library axis.
  - 12 continuity-sweep axes, each ≥ 2° from every library axis.
- **Stage-2 gate:** 72/72 accepted per variant. Otherwise stop.

## Evaluation (pre-specified)

**Surrogates:**
- **Primary:** the frozen GPv2 (antipodal, A20), applied identically to every design and budget.
- **Sensitivity:** the FIXED surrogate (ℓ = 0.5).

**Budgets:** n = 8, 12, 16, 20, 24, 28, 33.

**Metrics on the 60 references:** mean, RMS, P95 and max relative error; sphere-area-weighted error; error by region; 90%/95% coverage; RMS standardised residual. The 12 sweep references are reported separately.

**Primary endpoint:** area-weighted mean relative error at n = 33 under the primary surrogate.

**BO advantage over axial maximin is claimed for a material only if all four hold:**
1. BO is better in ≥ 8 of 9 paired starts at n = 33 (primary endpoint).
2. The paired median improvement exceeds δ_num = the S1b-min load-step |Δr/r| for that material (B1 0.010%, B2 0.034%), in percentage points.
3. Conditions 1 and 2 also hold under the FIXED sensitivity surrogate.
4. The unweighted mean error on the 60 references is also lower in the median.

**Statistics:** Wilcoxon signed-rank on the 9 pairs (reported, no claim threshold). BO vs random: Mann–Whitney, 9 vs 10.

**Response dependence:** first divergent test between B1 and B2 BO sequences at the same start, and against the homog V2-A20 sequences and axial maximin.

**Replay:** ReplayOracle must reproduce `tests.csv`, `steps.csv` and the acquisition scores bitwise.

## Stage gates

| Stage | Gate |
|---|---|
| 2 | Every reference accepted |
| 3 / 4 | Non-adaptive: all 19 designs reach 33 successful tests. BO: all 9 campaigns reach 33 tests. Replays identical |

Any failure stops the pipeline at that stage. No configuration may be changed after launch.
