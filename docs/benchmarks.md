# Benchmarks

**Physical** = MOOSE simulations. **Analytical** = exact closed-form oracle evaluations; these are *not* experiments.

All numbers are read from the frozen evaluation records and copied into `results/`. The summary table is `results/benchmark_summary.csv` (+ `.json`).

## Physical (MOOSE)

| Benchmark | Specimen | Methods | Budget, starts | References | Files |
|---|---|---|---|---|---|
| **B1 / B2** | 2-D plane-stress periodic RVE, 11 aligned ellipses (f = 0.165), J2 phases with linear hardening. B1: stiff, slightly stronger inclusion. B2: compliant, much weaker inclusion. h = 1/40. | V2 BO+MES, angular maximin, random | 33 tests (8 initial + 25 BO), 9 starts (s = 0, 1, 2, 3, 5, 6, 7, 8, 9); random: 10 seeds | 72 (B1) / 71 (B2) independent MOOSE directions | `results/physical/B1_B2_composites/` |
| **TL1** | Periodic pin-jointed elastoplastic truss: frame + kagome infill, 6 bar families | V2, angular maximin, random (campaign 3); frozen V3, V3-Signed (transfer study) | 64 tests (primary endpoint n = 33), 9 starts | 482 independent MOOSE directions | `results/physical/TL1_truss/` |

- **B1/B2 metric:** area-weighted mean |r̂/r − 1| at n = 33. The pre-specified claim rule needed ≥ 8/9 paired wins, a median improvement > δ, robustness under fixed hyperparameters, and mean support.
- **B1/B2 result:** V2 beats random (Mann–Whitney p = 0.0017 / 0.0004). An advantage over maximin was **not claimed** (B1 7/9, B2 6/9).
- **TL1 result:** V2 vs random 1.28× [1.07, 1.59], 9/9. Neither V2 vs maximin (1.11×, 7/9) nor frozen V3 vs maximin (1.11×, 7/9, Holm p 0.41) is a supported difference.

## Analytical

| Benchmark | Surfaces | Methods | Budget | Reference | Files |
|---|---|---|---|---|---|
| **H1 holdout** (pre-registered, independent) | 24 Barlat Yld91 surfaces: m ∈ {2, 6, 8} × anisotropy s ∈ {0.15, 0.30} × 4 realisations; random in-plane orientation | V3, V2, angular maximin, random, V3_ISO ablation | 64 (primary n = 33) | 20,000 rotated equal-area directions, exact | `results/analytical/H1_Yld91_holdout/` |
| **A1–A6 complexity study** | von Mises (A1), Hill (A2), icosahedral p-norms (A3, A4), random-axis p-norms (A5, A6; most localized) | V2 BO, angular maximin, random, BOFIXED (diagnostic) | up to 1,000 sequential observations, 4 starts per level | 20,000 Fibonacci directions, exact | `results/analytical/A1_A6_complexity/` |

**H1 results:**
- V3/maximin: 1.34× [1.13, 1.65] (19/24, Holm p = 0.0012) at n = 33, and 2.01× [1.69, 2.39] (24/24) at n = 64.
- V3/V2: 1.48× (23/24).
- Median errors at n = 33: V3 0.24 %, maximin 0.27 %, V2 0.29 %, random 0.47 %.
- V3_ISO ablation: 1.01× (12/24), so V3's gain is not attributable to material-specific responses.

**A1–A6:** the BO error decreases monotonically on every level. For example, A6 goes from 4.31 % (n = 8) to 0.028 % (n = 1,000). An advantage over maximin was not shown at large budgets (n ≥ 256; see `conclusion.json`).

## Selection trajectories

| File | Content |
|---|---|
| `results/physical/B1_B2_composites/bo_selection_sequences.csv` | V2, 18 campaigns × 33 tests (incl. recorded failed tests) |
| `results/physical/TL1_truss/selection_sequences.csv` | 5 methods × 9 starts × 64 |
| `results/analytical/H1_Yld91_holdout/selection_sequences.csv` | 5 arms × 24 surfaces × 64 |
| `results/analytical/A1_A6_complexity/BO_1000_step_selection_sequences.csv` | 6 levels × 4 starts × 1,000 |

## GP posterior replays

`results/gp_posterior_replay/` holds:
- the full per-step posterior mean and latent variance on a 2,562-direction rendering grid and at the references, for B1, B2, TL1 V2 and TL1 V3 (start 0), as `.npz`;
- the per-step error curves for these and for H1/A6 (`error_vs_n.csv`);
- the verification records: sequences reproduced, score vectors bitwise equal to the recorded files.

## Seeds and configurations

| Rule | Value |
|---|---|
| Start s = 0 | first candidate index 0 |
| Start s ≥ 1 | `default_rng(7000 + s)` |
| Random baselines | `default_rng(9000 + k)` |
| H1 starts | 2000 + i |
| GP | ML-II with `random_state = 0` and 2 restarts |

Configurations: `research/.../closed_loop/configs/`, the protocols (`PROTOCOL*.md`), and `examples/*.yaml`.
