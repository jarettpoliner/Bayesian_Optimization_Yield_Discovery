# Mathematical method

_Source: the frozen project methods document (`research/yield_surface_moose/plane_stress/project_closeout_v1_20261120/METHODS.md`), reproduced unchanged below._


All notation refers to plane stress in Mandel coordinates.

## 1. Yield surface representation

- **Loading direction:** a unit Mandel vector **d** = (σ_xx, σ_yy, √2 σ_xy)/‖·‖.
- **Yield radius:** r(**d**) > 0 is the stress magnitude at initial yield along **d**.
- **Antipodal symmetry:** r(**d**) = r(−**d**). It was verified for the MOOSE oracle (|Δr/r| ≤ 2.3e-14) and holds exactly for every analytic family.
- **Loading axes:** candidates are therefore unoriented axes {**d**, −**d**}.

**MOOSE oracle (B1/B2).**

- **Model:** 2-D plane-stress periodic RVE (L = 1), 11 aligned ellipses (inclusion volume fraction 0.165), h = 1/40, corner-pinned.
- **Constitutive law:** J2 phases with linear isotropic hardening.
- **Solver:** PJFNK/LU.
- **Yield definition:** r is the load radius at which the volume-averaged equivalent plastic strain reaches 1e-4, after 5 % ramping and a 5-sub-step refinement.

**Analytic families** (exact, noise-free; all convex and antipodal):

| Family | Definition | Used in |
|---|---|---|
| p-norm of norms | g(z) = (q(z)^p + λ Σₖ (wₖ\|aₖ·z\|)^p)^{1/p}, q = √(zᵀMz) | complexity study, A1–A6 |
| LSF-v2 | the same form with one weak orientation; the critical resolved stress is set from the base support function | LSF-v2 |
| **Yld91** (H1 holdout) | g(σ) = ((\|S₁−S₂\|^m + \|S₂−S₃\|^m + \|S₃−S₁\|^m)/2)^{1/m}, where Sᵢ are the principal values of s′ = L(a, b, c, h)·σ in a rotated material frame | H1 |

In every family r = R0/g. Yld91 is convex for m ≥ 1.

## 2. Surrogate: frozen GPv2 (V2-A20), used by every method and by every evaluator

**Data and model.**

- Data: (**d**ᵢ, yᵢ = log rᵢ), i = 1…n.
- Mean: constant, equal to median(y).
- Augmentation: each real observation is accompanied by its antipodal pseudo-observation (−**d**ᵢ, yᵢ), with extra noise (1e-4)².
- Kernel: k(**x**, **x**′) = a·m_ℓ(‖**x** − **x**′‖) + σ_w² δ, where m_ℓ is Matérn ν = 3/2: m_ℓ(ρ) = (1 + √3ρ/ℓ)·exp(−√3ρ/ℓ).

**Bounds:** ℓ ∈ [0.05, 5], a ∈ [1e-5, 10], σ_w² ∈ [1e-6, 0.1]. Jitter is 1e-8.

**Schedule A20:**

- **n < 20 real observations:** fixed ℓ = 0.5, a = var(y_aug), σ_w² = 1e-3.
- **n ≥ 20:** ML-II with 2 restarts, refitted every step.

**Posterior.** Let K = k(X, X) + Λ, where Λ holds the noise, jitter and pseudo-noise terms. Then:

  μₙ(**x**) = m₀ + k(**x**, X) K⁻¹ (y − m₀),
  σₙ²(**x**) = a − a²·m(**x**, X) K⁻¹ m(X, **x**)   (latent variance).

**Key property.** σₙ² does not contain y. At fixed (a, ℓ, σ_w²) the variance, and any score built only from it, is response-independent. With A20, y enters only through var(y) (n < 20) and θ̂(y) (n ≥ 20).

## 3. Selection methods (candidates = 4,096-axis Fibonacci half-sphere, analytic; 207 axes from the V1 library, B1/B2)

**Shared initial design.** All methods start from the same 8 points: axial maximin (max-min axial angle min(θ, π − θ)) from the first axis of start s. For s ≥ 1, that first axis is `default_rng(7000 + s)`.

**Angular (axial) maximin.** Greedy: the next axis maximises its minimum axial angle to the current design. It is response-independent.

**Random.** The initial 8, then unique axes in `default_rng(9000 + s)` order.

**V2 BO+MES.** **x**_{n+1} = argmax over untested axes of σₙ²(**x**) (maximum entropy sampling). Ties are deterministic: the lowest index among scores within rel 1e-10 / abs 1e-12 of the maximum.

**V3 LRV (v3.0, frozen).** **x**_{n+1} = argmax σₙ²(**x**) · ρₙ(**x**). The risk is built as follows.

1. **Block leave-one-out residual for real observation i**, removing the block B = {i, its pseudo-image}, with Q = K⁻¹:
   - e_B = (Q_BB)⁻¹ (Q y)_B and C_B = (Q_BB)⁻¹;
   - uᵢ = e_B[i] / √C_B[i, i].
2. **Interpolated risk:**
   - ρₙ(**x**) = Σᵢ cᵢ(**x**) uᵢ² / Σᵢ cᵢ(**x**) / meanᵢ(uᵢ²);
   - cᵢ(**x**) = max(m₀.₅(**x**, **x**ᵢ), m₀.₅(**x**, −**x**ᵢ)).
   - The interpolation length is fixed at 0.5, V2's pre-threshold value, because ML-II ℓ saturates at 5.
3. **Safeguards:** ρ ≡ 1 if the uᵢ are non-finite or mean(u²) ≤ 1e-300; ρ ≥ 1e-12. If the uᵢ are spatially constant, ρ ≡ 1 and V3 = V2 (verified bitwise).

## 4. Evaluation (identical for every method within a study)

- **Fit and prediction:** the frozen GPv2 (A20) is refitted from scratch on each design's first n observations and predicts at independent references.

| Study | References |
|---|---|
| B1/B2 | 60 area-stratified MOOSE references (≥ 3° from the library) + 12 sweep references (B2: 11) |
| Analytic studies | 20,000 equal-area Fibonacci points |
| H1 | a new rotated 20,000-point set |

- **Error:** |r̂/r − 1|.
- **Global error:** the area-weighted mean (for equal-area points, the plain mean).
- **Regional metrics:** feature region F (LSF-v2: weakening ≥ 1 %); high-curvature region (top 10 % of κ₁r); maximum error; 90 %/95 % coverage of the predictive intervals using latent + white variance.

## 5. Statistics

**Paired statistic:** L_j = log(err_other,j / err_method,j).

**Tests and estimators:**

- the exact two-sided Wilcoxon signed-rank test;
- the Hodges–Lehmann estimator exp(median of Walsh averages);
- the exact 95 % CI from the signed-rank null distribution (dynamic programming), [W₍c₎, W₍N−c₋₁₎] with c = max{t : P(T ≤ t) ≤ 0.025};
- wins = #(L > 0).

**Multiplicity.**

- **H1:** Holm over H1–H3.
- **Complexity study:** Holm over the 7 budgets n ≥ 256.
- **B1/B2:** a pre-specified claim rule (≥ 8/9 pairs, median improvement > δ_num, robust under FIXED, mean supports), with Wilcoxon reported. Mann–Whitney is used vs random.

**Confirmed practical superiority (LSF-v2 Outcome A and H1):** p (Holm where applicable) < 0.05, **and** HL ≥ 1.10, **and** the lower CI > 1, **and** ≥ 16/24 wins.

## 6. Diagnostic constructions (no method claims)

| Construction | Definition | Status |
|---|---|---|
| Oracle error-greedy design ORC-E | next = argmax over the library of the *true* \|r̂ₙ/r − 1\| | a reference design, not implementable |
| Response ablations | the same controller observing a smooth surrogate surface: DBU base R0/q, or isotropic J2 for H1 V3_ISO. The resulting design is re-observed on the actual surface | |
| Hill-quadratic fit | 1/r² = **z**ᵀM**z**, by linear least squares on r² | |
