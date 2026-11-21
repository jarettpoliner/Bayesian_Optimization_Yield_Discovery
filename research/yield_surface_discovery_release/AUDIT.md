# Release audit

Audit date: 2026-11-15. Saved result files were treated as the source of truth; earlier summaries were checked against them. All paths below are relative to the research directory `exp3/fuhg_benchmark/` unless stated otherwise.

## 1. Frozen algorithms

**Frozen practical method (physics-informed):**
- 8 angular-maximin initial tests from the candidate library (geometry only);
- direct Matérn-3/2 GP on log yield radius (constant median mean; no analytical template). Hyper-parameters are fixed below 20 training points (amplitude = var(y), length scale 0.5, white noise 1e-3), with ML-II from 20 points;
- maximum predictive-uncertainty (predictive-variance) sampling from the same GP;
- known shear-mirror symmetry (σxy → −σxy) as pseudo-observations;
- expanded physical loading library (Pool 2).

Registered as **same-surrogate Arm C** (`protocol/same_surrogate_ablation_protocol.json`) and confirmed against the symmetry-aware cleanup (`protocol/symmetry_aware_cleanup_protocol.json`, decision: keep Arm C).

**Deployment details** (recommended by the registered cleanup "keep" outcome; bookkeeping only):
- deduplicated candidate list (416 labels → 414 directions);
- deterministic tie handling (relative tolerance 1e-10, lowest candidate index).

This exact combination was never run as a registered arm. In the registered runs it would have changed only choices among tied candidates (16 of 100 Arm C selections).

**Frozen generic blind baseline:** the same algorithm without the shear-mirror symmetry (same-surrogate **Arm A**).

## 2. Authoritative files

| Item | Authoritative file(s) |
|---|---|
| Pool 1 / 2 / 3 definitions | `protocol/frozen_sets.npz` (Pool 1 + held-out rays + regions); `protocol/pool_experiment_sets.npz` (Pools 1–3); `protocol/pool2_candidates.csv`; `protocol/pool2_deduplication.csv` |
| Pool 1 / 2 / 3 results (template-informed + mirror) | `results/pool_experiment_per_run.csv`, `results/pool_experiment/per_ray/` |
| Same-points Fuhg comparison | `results/same_points_per_run.csv`; `results/blind_gp_mes/blind_same_points_fuhg.csv` |
| Blind GP study | `results/blind_gp_mes/{blind_gp_mes_per_run.csv, blind_decomposition.csv, blind_decisions.json, BLIND_GP_MES_REPORT.md}`; `protocol/blind_gp_mes_protocol.json` |
| Same-surrogate ablation (frozen method = Arm C; blind = Arm A) | `results/same_surrogate_ablation/{same_surrogate_per_run.csv, trajectories/, same_surrogate_decisions.json, response_dependence_*.csv, SAME_SURROGATE_REPORT.md}`; `protocol/same_surrogate_ablation_protocol.json` |
| Symmetry-aware cleanup | `results/symmetry_aware_cleanup/{cleanup_per_run.csv, cleanup_decisions.json, SYMMETRY_AWARE_CLEANUP_REPORT.md}`; `protocol/symmetry_aware_cleanup_protocol.json` |
| Direct Matérn-3/2 GP | `blind_lib.fit_blind` (registered); `ss_lib.fit_surface` (verified replica) |
| Angular maximin | `../analytical_library/ana_lib.maximin` |
| Predictive-uncertainty acquisition (same GP) | `ss_lib.entropy_scores` / `ss_lib.run_arm` |
| Shear-mirror augmentation | `blind_lib.fit_blind(sym=True)` |
| Metrics | `bench_lib.metrics`; `blind_lib.shear_bands` |
| Mesh / export | `bench_lib.mesh_eval` → `../edge_followup/edge_lib.export` (pinned icosphere level 5) |
| Earlier figure scripts | `f4_report.py`, `f5_init_experiment.py`, `f7_pool_experiment.py`, `f8_final_figures.py` (template-informed framing; superseded by `release/`), `b4_figures.py`, `s3_figures.py`, `sc4_figures.py` |

The portable package re-implements the GP, maximin, acquisition, symmetry, controller and metrics, and is verified against these files:
- bit-exact re-run of all 10 registered Arm A / Arm C sequences in the pinned environment;
- GP predictions identical to `blind_lib.fit_blind`;
- region masks identical to the frozen ones;
- continuous global p95 identical to ≤ 2e-16.

The pinned mesh exporter was **not** ported; it depends on the research geometry stack. The package reports continuous-prediction metrics instead. On the frozen same-surrogate runs, the difference between mesh-scored and continuous global p95 is ≤ 0.040 percentage points (`data/frozen_results/PROVENANCE.json`).

## 3. Known caveats

1. **Frozen practical method registered on Pool 2 only.** Pool 1 / Pool 3 comparisons with the same reconstructor and symmetry come from the registered decomposition arm: the direct GP + mirror fitted to the *frozen symmetry-assisted test sequences* (Hill-calibration initial tests + symmetric RBF selector). Figure 3b and the main table label this explicitly. The frozen practical method itself was **not** run on Pools 1 and 3, and no new runs were made for the release.
2. **Arm C's 5 starts are 4 distinct designs.** Starts 2 and 3 reach the same initial set, and their results are identical. Arm A's starts 2 and 3 share the initial set but diverge after hyper-parameter learning.
3. **Response-awareness.** Selections respond to measured yield values only through GP hyper-parameter learning: from 20 training points, i.e. test 21 without symmetry and tests 11–16 with it. Before that they are geometry-driven (registered diagnostic).
4. **Exact reproduction needs the pinned environment.** The registered runs used Python 3.12.9 with numpy 1.26.4, scipy 1.15.2, scikit-learn 1.6.1, pandas 2.2.3 and matplotlib 3.10.1.
   - On the workstation these are the base-environment versions with `PYTHONNOUSERSITE=1`. A user-site directory on the same machine provides numpy 2.5.4 / scipy 1.18.1 / scikit-learn 1.9.1, which shadow them when the user site is enabled.
   - Under those versions the registered legacy round-off ties resolve differently (ulp-level differences in radii), so legacy sequences diverge. Deterministic mode, replay checks and metrics agree.
5. **Single target.** All primary numbers are for one target (Mg–4%Li, CB2004).
   - The cross-surface runs (blind study) are secondary.
   - Tsai-Wu remains representation-limited (needle-shaped surface; a smooth stationary GP).
6. **Pool 3 is not an upper bound.** The template-informed Pool 2 beats it at 24–28 tests.
7. **MOOSE has never been run.** The adapter is tested against a fake executable only. The input template is an unverified draft.

## 4. Errata (earlier summaries that the saved files contradict)

| Where | Stated | Correct (from saved files) |
|---|---|---|
| `SAME_SURROGATE_REPORT.md` (note 2, §9) and `same_surrogate_protocol.json → implementation_notes` | Pool 2 contains 36 duplicated pure-shear directions | **2** exact duplicates (`shear+160` ≡ `shear-70`, `shear-160` ≡ `shear+70`). An erratum is appended to the report; the JSON is left unchanged as a record. |
| `BLIND_GP_MES_REPORT.md` (first draft) | 1,713 valid meshes | 2,883 new meshes (1,623 registered + 1,260 post hoc); corrected in place on 2026-11-15 |
| My session's first version check | numpy 2.5.4 / scikit-learn 1.9.1 environment | That was the user-site environment; the registered environment is numpy 1.26.4 / scikit-learn 1.6.1 (caveat 4) |
| README §10, `pool_experiment_protocol.json` | Pool 3 = "upper bound" | Arbitrary-direction **reference**, not an upper bound |
| README §§3–5, early protocols | "MES" | Maximum predictive-uncertainty (entropy) sampling, then on a separate RBF GP. It never meant max-value entropy search. |

## 5. Exploratory only (never use for primary conclusions)

- `posthoc_*` arms in `results/blind_gp_mes/` (`b2y_cross_decomposition_posthoc.py`) and `tsai_wu_same_points_posthoc_sym/` (`b2x_tsaiwu_sym_posthoc.py`);
- `f6_exploratory_reachability.py`, `results/exploratory_reachability.csv`;
- README §6, the arbitrary stress-direction experiment (separate);
- `../edge_followup/` (branch model, flank probes);
- `moose_validation/results/controller_tolerance_check.csv` (analytical oracle, not MOOSE);
- the cross-surface blind/template runs (secondary; post-hoc decomposition).

## 6. Must not be used for primary conclusions

- `figures/final_paper/`, `tables/final_paper/` and `f8_*`: they present the template-informed method as primary, which is superseded.
- Any "upper bound" wording for Pool 3.
- Any description of the physics-informed method as "blind".
- `moose_validation/MOOSE_ORACLE_SPEC.md` Appendix A as a validated deck (it is untested).
- `../explainer/`, `../analytical_library/`, `../paper_package/`: earlier template-bank studies with a different setup. They are used here only as context for the Tsai-Wu limitation.
