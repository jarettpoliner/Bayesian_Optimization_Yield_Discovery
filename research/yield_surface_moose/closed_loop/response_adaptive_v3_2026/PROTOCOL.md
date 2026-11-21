# V3 response-adaptive acquisition pilot: frozen protocol

**Written:** 2026-11-20.

**Freeze:** hashed in `PROTOCOL.sha256` **before any V3 campaign was run on the gate or benchmark populations** and before any V3 reconstruction error was computed.

**Basis:** `closed_loop/methodology_diagnosis_2026/NEXT_EXPERIMENT_PLAN.md`.

## 1. Scope and honesty statement

- **One change.** This pilot tests one change to the frozen V2 BO+MES controller: the acquisition score. Nothing else changes, including:
  - the GP, Matérn-3/2 kernel, A20 schedule and hyper-parameter bounds;
  - the antipodal symmetry, the 4,096-axis library and the initial designs;
  - the tie rule, the oracles, the references and the evaluator.
- **Imports are read-only.** The frozen V1/V2 sources, `cs_common.py`, `surface_lsf.py`, the DBU population and every release are imported read-only and hash-checked.
- **No MOOSE, no 1,000-point campaign, no new microstructure, no acquisition search.** Exactly one acquisition is implemented and evaluated, and it is not tuned after any result.
- **Status of the evidence:** the LSF-v2 population and its oracle-headroom analysis informed the design of V3 (methodology diagnosis). The Stage 2 result is therefore **development-stage evidence, not independent confirmation**.

## 2. The V3 method (`v3_lrv.py`, VERSION `v3.0-LRV`)

**Score maximised over untested library axes:** s(x) = σ²ₙ(x) · ρₙ(x).

**Fit.** Identical to `V2Controller.fit`:
- A20: fixed ℓ = 0.5, amp = var(y), noise 1e-3 below 20 real observations; ML-II (2 restarts) from 20;
- antipodal pseudo-observations;
- held hyper-parameters.

σ²ₙ is the frozen latent posterior variance, and only past observations are used.

**LOO.** For each real observation i, the block {xᵢ, −xᵢ (its pseudo-observation)} is left out in closed form, at the current hyper-parameters:
- Q = K⁻¹, from the fitted Cholesky factor (including white noise, jitter and pseudo-noise);
- e_B = (Q_BB)⁻¹ (Q y)_B and C_B = (Q_BB)⁻¹;
- uᵢ = e_B[real] / √C_B[real, real].

**Interpolation.** ρₙ(x) = Σᵢ cᵢ(x) uᵢ² / Σᵢ cᵢ(x), with cᵢ(x) = max(m(x, xᵢ), m(x, −xᵢ)), where m is Matérn ν = 1.5 with length **0.5**.

**Normalisation.** ρₙ ← ρₙ / meanᵢ(uᵢ²). This sets the score scale only; the argmax is unchanged.

**Floors and safeguards** (numerical only):
- if any uᵢ or mean(u²) is non-finite, or mean(u²) ≤ 1e-300, then ρₙ ≡ 1 (the V2 behaviour, flagged `fallback`);
- LOO variance floor 1e-300;
- ρₙ floor 1e-12.

**Ties and exclusions.** The frozen deterministic `select` (rel 1e-10, abs 1e-12); tested and unavailable axes are set to −∞.

**Reduction.** If the uᵢ are spatially constant, ρₙ ≡ 1 and V3 = V2 BO+MES (tested in Stage 0).

**Free parameters:** none.

**Documented deviation from NEXT_EXPERIMENT_PLAN.md, made before any test.** The plan wrote cᵢ with "the GP's own Matérn correlation at the current ℓ".
- **Problem:** in the ML-II phase, ℓ sits at its upper bound 5.0 in most steps (methodology diagnosis §3.2). With ℓ = 5, m ≥ 0.846 between *any* two axes (maximum chordal distance 2), so the weights would be almost uniform and V3 would collapse to MES by construction.
- **Change:** the interpolation therefore uses V2's own frozen pre-threshold length scale, 0.5 (`GPv2().length_scale`). This keeps the method parameter-free.
- **Basis:** this follows from the kernel form and the known bound alone. No V3 output was used.

## 3. Stage 0: implementation checks (`stage0_tests.py`; no reconstruction performance computed)

- **T1, reduction:**
  - (a) with `uniform_risk=True`, the sequences and scores equal the saved V2 BO+MES campaigns bitwise (LSF-v2 r00 to 64; DBU r00 and r01 to 33);
  - (b) with spatially constant LOO residuals passed through the full risk path, ρ is within 1e-14 of 1 and the sequences are the same.
- **T2, LOO:** the closed-form block LOO equals a brute-force refit (rel ≤ 1e-8) at a FIXED and an ML-II state.
- **T3, determinism:** repeated runs are identical.
- **T4, information:** perturbing observations after the 20th leaves selections 1–21 unchanged; the controller makes exactly `budget` oracle calls; it receives only the library and the oracle.
- **T5, hashes:** the frozen-source hashes are unchanged and `PROTOCOL.sha256` verifies.

**Any failure stops the study.**

## 4. Stage 1: response-adaptivity gate (`gate_v3.py`; frozen DBU rule, unchanged)

**Realisations:** DBU j = 0–3 (`adaptivity_dbu_2026/population.json`), with their starts and initial 8. Budget 33.

**Arms:**

| Arm | Status |
|---|---|
| A — V3 on the actual radii | new |
| B — V3 on the smooth base R0/q (response ablation) | new |
| C — angular maximin | saved DBU campaign, reused |
| Frozen V2 BO+MES | saved, reused; reported, not in the rule |

**Measures:**

- D = |A selections 9–33 not in B selections 9–33|.
- N = number of selections 9–33 with weakening ≥ 0.1 %.
- **Rule per realisation:** D ≥ 5 **and** N_A ≥ N_B + 2 **and** N_A ≥ N_C + 2.
- **GATE PASS iff ≥ 3 of 4.**

**Descriptive (not in the rule):**
- shoulder, core and F counts;
- enrichment;
- the feature-region and global error trajectories (frozen evaluator);
- the share of steps whose choice differs from MES.

**If the gate fails: STOP.** Stage 2 is not run, and V3 is not modified.

## 5. Stage 2: LSF-v2 benchmark (`run_v3.py lsf`, `evaluate_v3.py lsf`, `analyze_v3.py`; only if the gate passes)

**Population:** the 24 accepted LSF-v2 realisations, with starts j and the same initial 8 as the saved campaigns. Budget 64 (primary n = 33, secondary n = 64).

**New campaigns:** V3 and V3_ABL (response ablation, diagnostic arm).

**Comparators:** the saved, unchanged V2 BO+MES, angular-maximin and random campaigns. Their frozen evaluator rows are reused.

**Evaluator:** a verbatim copy of `evaluate_stress.task` and `analyze_stress.paired`. The copy must reproduce the frozen rows to CSV round-off (≤ 1e-14) and the frozen LSF-v2 primary result exactly (`evaluation/selfcheck.json`). Run it with the `ysd` interpreter that ran the frozen studies.

**Primary endpoint:** n = 33, feature-region mean |error|, V3 vs angular maximin, L = log(err_AXMM / err_V3).
- **OUTCOME A iff** p < 0.05 **and** HL ratio ≥ 1.10 **and** lower 95 % CI > 1 **and** ≥ 16/24 wins (the frozen LSF-v2 criteria).

**Guard** (pre-declared in the plan): n = 33 global mean |error|, V3 vs AXMM; the lower 95 % CI of the HL ratio must be ≥ 0.90.

**Practically meaningful** (interpretive, pre-declared): OUTCOME A **and** the guard is met **and** the median err_AXMM/err_V3 ≥ 1.25 at n = 33, i.e. a ≥ 20 % error reduction.

**Secondary (descriptive, nominal p-values, no claim):**
- n = 64 (Outcome-A criteria applied descriptively);
- feature max, global mean and global max;
- comparisons vs V2 BO+MES, random and V3_ABL;
- sample efficiency (n to feature error ≤ 1 / 0.5 / 0.25 %);
- concentration in F (frozen nearest-reference rule);
- response ablation (selection differences; error ratio);
- the share of steps differing from MES;
- calibration;
- computational overhead;
- the capture fraction of the oracle headroom (realisations 0–7).

**If OUTCOME A holds**, an independent, pre-registered analytic holdout study is *prepared* (protocol only) and not executed.

## 6. Integrity

- **Freeze:** `PROTOCOL.sha256` freezes this file, `v3_lrv.py`, `run_v3.py`, `evaluate_v3.py`, `gate_v3.py`, `analyze_v3.py` and `stage0_tests.py`.
- **After the freeze:** any change to these files is a documented amendment, with the original kept.
- **Checkpointing:** campaigns are checkpointed (a completed `run_info.json` is skipped; incomplete directories are moved to `runs_incomplete/`, never deleted). Every computation is deterministic and single-threaded BLAS.
- **Final audit:** `audit_v3.py` re-verifies the frozen sources, every earlier manifest and protocol freeze, this freeze, a bitwise replay of one V3 campaign, and that nothing outside this directory changed.
