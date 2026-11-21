# Complexity-scaling study 2026: frozen protocol

Written and hashed on 2026-11-19, **before any comparative campaign was run**. It is separate from the final B1/B2 benchmark, which is complete, verified and untouched.

## Research question

Does the relative sample efficiency of BO+MES (frozen V2-A20) over axial (angular) maximin improve as the measured effective yield-surface complexity and the experimental budget increase?

## Frozen methods (imported, never modified)

| Arm | Definition |
|---|---|
| **BO** | `closed_loop/v2/bomes_v2.V2Controller` with `V2Settings(schedule="A20", symmetry="antipodal", budget=1000)`.<br>Matérn-3/2 GP on log r; antipodal pseudo-observations; fixed hyper-parameters below 20 real observations; ML-II at every observation from 20 on (2 restarts).<br>Acquisition: maximum latent predictive variance (= max entropy).<br>Ties: deterministic lowest index. |
| **AXMM** | `bomes_v2.axial_maximin` on the same library, with the same first axis as BO start k. Its first 8 axes equal the BO initial design. |
| **RANDOM** | The same 8 initial axes as BO start k, then `default_rng(9000+k).permutation` of the remaining axes. Unique axes only; no antipodes are possible (see the library). |
| **BOFIXED** (diagnostic only, not a primary arm) | V2Controller with schedule FIXED: ℓ = 0.5 and fixed amplitude/noise forever. Its selections are response-independent by construction, so it isolates the effect of ML-II hyper-parameter adaptation. |

## Candidate library (frozen)

- **Size:** 4,096 unoriented axes = the z3 > 0 half of a 8,192-point Fibonacci sphere in Mandel space.
- **sha256** of the rounded coordinates: `8e1654f80a0240d0857de2abbbac455501c5fd1fabbbdb633f8b1e368ef0a545`.
- **Spacing:** minimum axial separation 1.16°.
- **No antipodes:** no two library members are antipodal.
- **Shared:** all three arms use exactly this pool.

## Budgets, seeds and evaluation

- **Budgets and checkpoints:** 8, 16, 33, 64, 128, 256, 512, 1000 real observations. Thresholds are additionally evaluated on the grid n ∈ {8, 12, 16, 20, 24, 28, 33, 40, 48, 64, 80, 96, 128, 160, 192, 256, 320, 384, 512, 640, 768, 896, 1000}.
- **Starts / seeds:** k = 0, 1, 2, 3 for every level and every arm. This gives matched pairs (BO_k vs AXMM_k vs RANDOM_k share the first 8 axes).
- **Surrogate used to score every design:** the same frozen GPv2(ANTIPODAL, A20), refit from scratch on the first n observations of each design. A FIXED-schedule surrogate is reported as a sensitivity analysis.

## Analytic stress-test family (exact ground truth; synthetic oracle, never presented as MOOSE)

The radius is r(d) = R0 / g(d), with the gauge

g(z) = (q(z)^p + λ Σ_k (w_k |a_k·z|)^p)^(1/p),  q(z) = √(zᵀMz).

g is a norm, so each surface is convex and exactly antipodally symmetric. R0 sets the median r to 200 MPa.

| Level | M | p | λ | a_k |
|---|---|---|---|---|
| A1_J2 | J2 | 2 | 0 | – |
| A2_hill | Hill | 2 | 0 | – |
| A3_ico_p6 | J2 | 6 | 0.6 | icosahedral 6 |
| A4_ico_p16 | J2 | 16 | 1 | icosahedral 6 |
| A5_rand10_p32 | Hill | 32 | 1 | 10 random (seed 11) |
| A6_rand16_p64 | Hill | 64 | 1 | 16 random (seed 23) |

## Complexity descriptors and classes (predefined; applied identically to analytic and MOOSE surfaces)

The descriptors are computed on a 20,000-point equal-area reference sphere:

- max/min r;
- std of log r;
- dimensionless principal curvatures κ·r (κ by finite differences of the star surface; κ·r = 1 for a sphere);
- p99 and maximum of κ₁·r;
- coefficient of variation of κ₁·r;
- share of total κ₁·r held by the top 10% of directions (concentration);
- angular correlation length of log r (axial separation where the correlation drops below 1/e);
- antipodal asymmetry;
- convexity (min κ₂·r ≥ −1e-3, the finite-difference tolerance).

**Class = by p99 of κ₁·r:**

| Class | p99 of κ₁·r |
|---|---|
| C1_smooth | < 3 |
| C2_moderate | 3–6 |
| C3_localized | 6–12 |
| C4_highly_localized | ≥ 12 |

The measured classes (`analytic_descriptors.json`) are:

| Class | Levels |
|---|---|
| C1 | A1, A3 |
| C2 | A2, A4 |
| C3 | A5, A6 |
| C4 | none |

The level *names* are not a complexity order; the measured class is.

## Physical (MOOSE) configurations (predefined; same frozen two-phase J2 + linear hardening oracle r1; B1 phase properties)

| Config | Geometry |
|---|---|
| P1 | Homogeneous (analytic J2 is its exact reference; also HOMOG2 MOOSE runs exist) |
| P2 | Single circular inclusion, f = 0.15 |
| P3 | 3 × 3 regular circles, f = 0.15 |
| P4 | Clustered circles (S3-type RSA, geometry seed 4401), f ≈ 0.15 |
| P5 | Multiscale ellipses (sizes 0.03–0.11, aspect 1–4, random orientations, seed 5501), f ≈ 0.18 |
| P6 | Aligned high-aspect ellipses (a = 0.18, b = 0.03, θ = 0°, seed 6601), f ≈ 0.15 |

- **Pilot reference per configuration:** 24 area-uniform axes (z3 > 0 half of a 48-point Fibonacci sphere), solved by the frozen oracle r1 (corners pin, h = 1/40, dr_rel 0.05, n_refine 5). The complexity class is assigned from the GP fitted to these 24 observations (a coarse estimate, reported as such).
- **Mesh check:** P6 (thinnest features) at h = 1/80 on 2 axes.
- **Second geometry seeds** (P4: 4402, P5: 5502, P6: 6602) are predefined but deferred: they are part of the full campaign, not the pilot.
- **Launch rule:** the full physical BO / AXMM / RANDOM campaign to 1,000 observations is **NOT** launched in this session unless the measured pilot cost fits the compute already used by the approved B1/B2 benchmark (≈ 23 CPU-h). Otherwise only an exact resource estimate is reported.

## Metrics (independent held-out references, never selected by any arm)

- **Analytic:** the 20,000 Fibonacci directions with exact r. Their minimum angle to any library axis is > 0.
  - mean |relative error| (equal-area, hence area-weighted);
  - RMS, P95 and max |error|;
  - error on the high-curvature region (top 10% of κ₁·r);
  - 90% coverage and RMS z;
  - mean predictive sd.
- **n to reach** 1%, 0.5%, 0.2% and 0.1% area-weighted error: the first grid n at which the median over starts falls below the threshold.

## Statistics and success criterion (pre-specified)

- **Paired comparison** of BO vs AXMM per level and budget over the 4 starts: the sign of the difference (count of wins) and the median paired improvement. With n = 4, a per-level Wilcoxon test cannot reach p < 0.125, so **no per-level significance is claimed**.
- **Pooled test:** over all 6 × 4 = 24 matched pairs at each budget, a Wilcoxon signed-rank test.
- **Trend test:** Spearman correlation between the level's κ₁·r p99 and the paired log error ratio log(err_AXMM / err_BO).
- **Conclusion: "BO's relative advantage increases with complexity"** is stated only if both hold at n ≥ 256:
  1. the trend correlation is positive with p < 0.05 (24 pairs);
  2. the pooled paired improvement on C3 levels is positive in ≥ 7 of 8 pairs.

  Otherwise the result is reported as **not supported**.

## Response dependence (recorded for every BO campaign)

- **Recorded per step:** hyper-parameters, predictive variance (the score vectors), the posterior mean on the candidates, selected axes and radii.
- **Divergence analysis:** the first divergent selection between levels at the same start, and against BOFIXED and AXMM.
- **Expected below 20 observations:** with fixed hyper-parameters the variance does not depend on y, so the selections must be identical across levels.
