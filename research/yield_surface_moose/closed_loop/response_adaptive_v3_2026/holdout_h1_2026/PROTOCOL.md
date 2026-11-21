# Holdout H1: independent confirmatory test of V3's global reconstruction advantage (frozen protocol)

**Written and hashed:** 2026-11-20, **before any controller was run on an H1 surface**.

## 1. Purpose and status

**Origin.** The V3 pilot (`../REPORT.md`) failed its frozen feature-localisation gate (0/4). That result stands unmodified. A post-gate, exploratory look at those 4 realisations suggested lower global error than angular maximin (4/4). This study tests that suggestion on **independent** surfaces.

**Fixed elements.** The frozen V3 (`../v3_lrv.py`, sha256 `6eaba9e5…`, checked at import), the frozen V2-A20, and all their rules are used unchanged.

**Hypothesis origin.** The hypothesis was formed on the DBU/LSF-v2 realisations. None of them, nor their orientations, starts or references, is used here.

## 2. Population (`surfaces_h1.py` → `population.json`, `surface_verification.json`)

**Family:** Barlat Yld91 plane stress (a Hosford function of a linearly transformed stress).
- It is convex for m ≥ 1, and even, hence antipodally symmetric.
- m = 2 is Hill-type quadratic; isotropic m = 2 is von Mises (checked to 4e-16).
- **It differs from every family used to design V3:** the cs_common p-norm-of-norms family A1–A6 and LSF-v2.

**Levels** (fixed from physical convention, not from algorithm performance):

| Factor | Values |
|---|---|
| Hosford exponent m (directional complexity) | 2, 6 (BCC convention), 8 (FCC convention) |
| Anisotropy spread s: log a, b, c, h ~ N(0, s²) | 0.15 (mild), 0.30 (strong) |
| Material orientation θ | U[0, π) |

- **Realisations:** 24, i.e. 6 levels × 4 realisations, with seeds 52000 + i.
- **Starts:** 2000 + i. These are fresh and have never been used before.
- **Replacement rule:** geometry only.
- **Verification:** all 24 were accepted on the first draw.
  - Midpoint convexity holds on 200,000 pairs. Antipodal symmetry is exact and every radius is positive.
  - Minimum principal curvature κ₂r is 0.021–0.403.
  - Anisotropy r max/min is 1.23–3.69, and κ₁r p99 is 2.1–6.9.
- **R0:** median radius 200.

**Candidate library:** the frozen 4,096-axis Fibonacci half-sphere (hash `8e1654f8…`), as in all analytic studies.

**References:** `references_h1.npy`, 20,000 equal-area Fibonacci points rotated by `Rotation.random(random_state=52999)`.
- This point set has never been used before.
- The minimum axial angle to the library is 0.017°, so no reference coincides with a candidate.
- sha256 is recorded in `population.json`.

**High-curvature region HC:** the top 10 % of references by κ₁r, per surface.

## 3. Arms (`run_h1.py`)

All arms use budget 64, the same library, the same start, and the identical initial 8 (`V2Controller.initial_design(start)`, which equals the first 8 of axial maximin; asserted).

| Arm | Description |
|---|---|
| **V3** | frozen v3.0-LRV, actual radii |
| **V2** | frozen V2-A20 BO+MES, actual radii |
| **AXMM** | frozen angular (axial) maximin |
| **RANDOM** | the initial 8, then `default_rng(9000 + start)` unique axes |
| **V3_ISO** (diagnostic, not a competing method) | frozen V3 observing **isotropic von Mises** radii, i.e. no information about the material. Its design is re-observed with the actual radii for evaluation. It answers whether the measured responses account for V3's behaviour |

## 4. Evaluation (`evaluate_h1.py`)

- **Fit:** the frozen GPv2 (antipodal, A20) is refit from scratch on each design's first n ∈ {8, 12, 16, 20, 24, 28, 33, 40, 48, 56, 64} actual observations.
- **Prediction and metrics:** prediction at the 20,000 H1 references; |r̂/r − 1|.
- **Recorded per n:**
  - global mean (area-weighted, since the references are equal-area) and max;
  - HC mean and max; out-of-HC mean;
  - 90 %/95 % coverage, global and HC;
  - RMS z; hyper-parameters.
- **Self-check:** the metric code reproduces the frozen LSF-v2 evaluator to 1e-16 (`evaluation/selfcheck.json`).

## 5. Statistics (`analyze_h1.py`; code identical to the frozen LSF-v2 analysis)

**Paired statistic:** L = log(err_other / err_V3) per realisation, 24 pairs. Tests: exact two-sided Wilcoxon signed-rank; the Hodges–Lehmann ratio exp(HL) with its exact 95 % CI; wins.

**Confirmatory family** (Holm, α = 0.05):

| Hypothesis | Comparison | Metric |
|---|---|---|
| **H1 (primary)** | V3 vs AXMM | global mean, n = 33 |
| H2 | V3 vs V2 | global mean, n = 33 |
| H3 | V3 vs AXMM | global mean, n = 64 |

**Verdict per hypothesis:**
- **CONFIRMED** (practical superiority) iff the Holm-adjusted p < 0.05 **and** exp(HL) ≥ 1.10 **and** the lower CI > 1 **and** wins ≥ 16/24.
- **INFERIOR** iff the upper CI < 1.
- Otherwise **NOT ESTABLISHED**.

**Trade-off guards** (V3 vs AXMM, n = 33): the HC mean error and the global max error each need a lower CI of exp(HL) ≥ 0.90. Failing either is an "unacceptable trade-off".

**Descriptive only** (nominal p, no claim):
- all metrics × {33, 64} × {AXMM, V2, RANDOM, V3_ISO};
- response attribution (V3 vs the V3_ISO design);
- moderators (Spearman of L_H1 with anisotropy, κ₁r p99 and m; per-level medians);
- error vs n;
- coverage;
- selection differences;
- computational cost.

**Sample size:** 24 pairs. The smallest attainable exact two-sided p is 1.2e-7, so Holm over 3 hypotheses is not limiting.

## 6. Resources (`preflight.json`)

| Item | Value |
|---|---|
| Campaigns | 120 |
| Estimated CPU | ≈ 12 CPU-min |
| Workers | 4, at `nice 19`, single-threaded BLAS |
| Machine | 24 logical CPUs |
| Stop rule | stop after the protocol if the estimate exceeds 2 CPU-h |

The concurrent truss workstream is never read or written.

## 7. Integrity

- **Freeze:** `PROTOCOL.sha256` freezes:
  - this file;
  - `surfaces_h1.py`, `population.json`, `surface_verification.json`, `references_h1.npy`;
  - `run_h1.py`, `evaluate_h1.py`, `analyze_h1.py`;
  - `preflight.json`, `evaluation/selfcheck.json`.
- **No changes after results:** no realisation is excluded, and nothing in the method, metric or criterion changes after any result.
- **Negative results are reported as such.**
- **Audit (`audit_h1.py`):**
  - frozen sources and the V3 hash;
  - designs exact;
  - a bitwise replay;
  - earlier manifests intact;
  - nothing outside this directory changed, other than the concurrent truss study, which is recorded and untouched.
