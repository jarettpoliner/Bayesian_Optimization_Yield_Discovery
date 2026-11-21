# Figure captions (final showcase)

Figures A–N, TL1_* and TL1V3_* are copied unchanged from final_integrated_release_v3_20261120 (captions below are verbatim from its captions.md). S1 is new, drawn from the frozen H1 median-error file.

**S1 — H1 sample efficiency.** Median (over 24 independent Barlat Yld91 surfaces) global relative radial error versus number of experiments for V3 (error-weighted predictive variance), V2 BO+MES, angular maximin and random sampling; analytic surfaces with exact references (no MOOSE). Annotation: paired Hodges–Lehmann maximin/V3 error ratios with exact 95 % CIs (pre-registered H1 and H3). Source: closed_loop/response_adaptive_v3_2026/holdout_h1_2026/analysis/error_vs_n_median.csv and comparisons.csv.

## A_workflow

Figure A. End-to-end closed-loop MOOSE Bayesian experimental-design workflow: candidate axes, shared initial design, MOOSE virtual yield test, GP surrogate, V2/V3 acquisition, budget loop, and independent-reference evaluation (schematic).

Files: `figures/A_workflow.svg`, `figures/A_workflow.pdf`, `figures/A_workflow.png`

Sources (sha256):


## B_microstructures

Figure B. Microstructures: composite RVE mesh and phases (a), MOOSE equivalent plastic strain at macroscopic yield under uniaxial σxx for B1 (b) and B2 (c), phase properties (d); TL1 periodic truss geometry (e) and MOOSE bar plastic strain at offset yield for three loadings (f–h).

Files: `figures/B_microstructures.svg`, `figures/B_microstructures.pdf`, `figures/B_microstructures.png`

Sources (sha256):

- `release_B1_B2_BO_benchmark_20261119/paraview/all_cases/B1_FB_B1_pool_t01.vtu` `793f7e33059ca25f…`
- `release_B1_B2_BO_benchmark_20261119/paraview/all_cases/B2_FB_B2_pool_t01.vtu` `942e54c4ffdb0cba…`
- `closed_loop/biphasic_benchmark/oracle/materials_frozen.json` `499cde8ef4fd170a…`
- `closed_loop/truss_lattice_benchmark_2026/lattice/TL1_frozen.json` `37cec8d48f349df6…`
- `closed_loop/truss_lattice_benchmark_2026/pilot/tl1_reference_table.csv` `5404f97a7b91c634…`
- `closed_loop/truss_lattice_benchmark_2026/pilot/refs_TREF_v1_20261120T171229/TR410/run_out.csv` `1395f2df9b9bd59d…`
- `closed_loop/truss_lattice_benchmark_2026/pilot/refs_TREF_v1_20261120T171229/TR356/run_out.csv` `9b79874c9061fb29…`
- `closed_loop/truss_lattice_benchmark_2026/pilot/refs_TREF_v1_20261120T171229/TR338/run_out.csv` `1a6e155f18d51658…`

## C_formulation

Figure C. Mathematical formulation: common Matérn-3/2 GP on log r with antipodal augmentation and the A20 schedule; V2 maximum-entropy acquisition; V3 LOO-residual-weighted variance acquisition.

Files: `figures/C_formulation.svg`, `figures/C_formulation.pdf`, `figures/C_formulation.png`

Sources (sha256):


## D_acquisition_step

Figure D. Sequential acquisition in the physical campaigns B1 (start 0, test 13) and B2 (start 0, test 25): GP posterior uncertainty before/after the selected test, recorded MES scores over the 207 candidate axes, tested axes and the selected axis (argmax verified against the recorded campaign).

Files: `figures/D_acquisition_step.svg`, `figures/D_acquisition_step.pdf`, `figures/D_acquisition_step.png`

Sources (sha256):

- `closed_loop/biphasic_benchmark/final_bench_20261119/campaigns/B1/bo_s0/tests.csv` `fe15cd2fc750b0aa…`
- `closed_loop/biphasic_benchmark/final_bench_20261119/campaigns/B1/bo_s0/acquisition_scores.npy` `4606e3b0ea00c753…`
- `closed_loop/biphasic_benchmark/final_bench_20261119/campaigns/B1/bo_s0/steps.csv` `cb539289f563c090…`
- `closed_loop/biphasic_benchmark/final_bench_20261119/campaigns/B2/bo_s0/tests.csv` `d7c50ca3004e409e…`
- `closed_loop/biphasic_benchmark/final_bench_20261119/campaigns/B2/bo_s0/acquisition_scores.npy` `17b8959e296435f4…`
- `closed_loop/biphasic_benchmark/final_bench_20261119/campaigns/B2/bo_s0/steps.csv` `fa470b0a5952216c…`

## E_reconstruction_vs_references

Figure E. GP reconstruction at n = 8, 16, 33 vs independent MOOSE references for B1 and B2 (BO start 0): slice through the continuity sweep (equal biaxial to σxy shear) with 95 % band and in-plane tests; 3-D reconstructed surface with reference errors. n = 64 not available (physical budget 33).

Files: `figures/E_reconstruction_vs_references.svg`, `figures/E_reconstruction_vs_references.pdf`, `figures/E_reconstruction_vs_references.png`

Sources (sha256):

- `closed_loop/biphasic_benchmark/final_bench_20261119/evaluation/references_B1.csv` `e33eae46484efc83…`
- `closed_loop/biphasic_benchmark/final_bench_20261119/campaigns/B1/bo_s0/tests.csv` `fe15cd2fc750b0aa…`
- `closed_loop/biphasic_benchmark/final_bench_20261119/evaluation/references_B2.csv` `7cf5c0b79b453ff4…`
- `closed_loop/biphasic_benchmark/final_bench_20261119/campaigns/B2/bo_s0/tests.csv` `d7c50ca3004e409e…`

## F_convergence

Figure F. Error vs number of MOOSE tests per physical specimen: B1 and B2 (9 BO, 9 maximin, 10 random designs; 60 area-stratified MOOSE references; n ≤ 33) and TL1 (9 matched starts per method, including frozen V3 and V3-Signed from the V3 transfer validation; 338 area-uniform MOOSE references; n = 8, 16, 33, 64). Frozen A20 GP for every reconstruction. Panels are not comparable in absolute error.

Files: `figures/F_convergence.svg`, `figures/F_convergence.pdf`, `figures/F_convergence.png`

Sources (sha256):

- `closed_loop/biphasic_benchmark/final_bench_20261119/evaluation/metrics_all.csv` `4d654913802a37d5…`
- `closed_loop/truss_lattice_benchmark_2026/final_comparison_2026/v3transfer_20261120/evaluation/metrics.csv` `d23ad20b4b7076c0…`

## G_sample_efficiency_ratios

Figure G. Paired (or, for B1/B2 vs random, unpaired) Hodges–Lehmann error ratios with exact CIs and adjusted p-values within each study: physical B1/B2 and TL1 (V2), analytic A1–A6 and LSF-v2 (V2), independent analytic H1 holdout (V3, with the V3_ISO ablation) and the physical TL1 V3 transfer test (frozen V3 confirmatory; V3-Signed exploratory). Filled markers: pre-specified verdict met. TL1 n = 33 V2 vs maximin: Holm p = 0.039 but 7/9 < 8/9 wins, so no supported difference; no V3 comparison on TL1 is supported.

Files: `figures/G_sample_efficiency_ratios.svg`, `figures/G_sample_efficiency_ratios.pdf`, `figures/G_sample_efficiency_ratios.png`

Sources (sha256):

- `closed_loop/biphasic_benchmark/final_bench_20261119/evaluation/metrics_all.csv` `4d654913802a37d5…`
- `closed_loop/truss_lattice_benchmark_2026/final_comparison_2026/campaign3_20261120/evaluation/paired.csv` `4a58e0eacebeea50…`
- `closed_loop/truss_lattice_benchmark_2026/final_comparison_2026/campaign3_20261120/evaluation/claims.json` `db64ea44b02bee4c…`
- `closed_loop/truss_lattice_benchmark_2026/final_comparison_2026/v3transfer_20261120/evaluation/paired.csv` `f44df2234788ecef…`
- `closed_loop/truss_lattice_benchmark_2026/final_comparison_2026/v3transfer_20261120/evaluation/claims.json` `53b6dfcf6833fc3b…`
- `closed_loop/complexity_scaling_2026/evaluation_r3/metrics_all.csv` `7d1b8f28ab13d0be…`
- `closed_loop/adaptivity_stress_test_2026/analysis/comparisons.csv` `0e876d8f4476ec70…`
- `closed_loop/response_adaptive_v3_2026/holdout_h1_2026/analysis/comparisons.csv` `32eb10c5dccb638b…`
- `closed_loop/response_adaptive_v3_2026/holdout_h1_2026/analysis/results.json` `072fb613d6a27537…`
- `project_closeout_v1_20261120/tables/T2_complexity_scaling.csv` `75176081bb537708…`

## H_H1_V3_and_ablation

Figure H. H1 holdout: error vs experiments for all arms (a), paired ratios with HL CIs and Holm p (b), and the shared equal-biaxial sampling bias of V3 and V3_ISO (c).

Files: `figures/H_H1_V3_and_ablation.svg`, `figures/H_H1_V3_and_ablation.pdf`, `figures/H_H1_V3_and_ablation.png`

Sources (sha256):

- `closed_loop/response_adaptive_v3_2026/holdout_h1_2026/evaluation/metrics_all.csv` `0249f96791e9673a…`
- `closed_loop/response_adaptive_v3_2026/holdout_h1_2026/analysis/results.json` `072fb613d6a27537…`
- `closed_loop/response_adaptive_v3_2026/holdout_h1_2026/analysis/post_geometry.json` `b7eee16d067f1b98…`

## I_TL1_complexity_mechanisms

Figure I. TL1 complexity and mechanisms: dominant yielding infill family over 338 area-uniform references (a); measured radius along the equal-biaxial→σxy-shear circle for TL1 (72 points) and B1/B2 (12/11 points) (b); measured complexity metrics computed with identical code (c).

Files: `figures/I_TL1_complexity_mechanisms.svg`, `figures/I_TL1_complexity_mechanisms.pdf`, `figures/I_TL1_complexity_mechanisms.png`

Sources (sha256):

- `closed_loop/truss_lattice_benchmark_2026/pilot/tl1_reference_table.csv` `5404f97a7b91c634…`
- `closed_loop/truss_lattice_benchmark_2026/pilot/complexity_TL1_vs_B1B2.json` `152a4c3fd1fb3e3a…`
- `closed_loop/biphasic_benchmark/final_bench_20261119/evaluation/references_B1.csv` `e33eae46484efc83…`
- `closed_loop/biphasic_benchmark/final_bench_20261119/evaluation/references_B2.csv` `7cf5c0b79b453ff4…`

## J_cost_vs_error

Figure J. Cost vs error: B1/B2 error vs MOOSE CPU-hours per design (median single-thread solve time from the campaign logs); H1 analytic controller wall time vs error at n = 64.

Files: `figures/J_cost_vs_error.svg`, `figures/J_cost_vs_error.pdf`, `figures/J_cost_vs_error.png`

Sources (sha256):

- `closed_loop/biphasic_benchmark/final_bench_20261119/evaluation/metrics_all.csv` `4d654913802a37d5…`
- `closed_loop/biphasic_benchmark/final_bench_20261119/campaigns/B1/nonadaptive/trials.jsonl` `f54535895a4b8cbd…`
- `closed_loop/biphasic_benchmark/final_bench_20261119/evaluation/bo_steps_all.csv` `ea1c30b0471c621a…`
- `closed_loop/biphasic_benchmark/final_bench_20261119/campaigns/B2/nonadaptive/trials.jsonl` `c97bee5f17bace1c…`
- `closed_loop/response_adaptive_v3_2026/holdout_h1_2026/analysis/results.json` `072fb613d6a27537…`
- `closed_loop/response_adaptive_v3_2026/holdout_h1_2026/analysis/error_vs_n_median.csv` `9e0196482968dd0c…`

## K_calibration_parity

Figure K. Actual-vs-predicted yield radii at n = 33 (B1, B2; median over the 9 BO designs, 95 % predictive interval) and empirical coverage vs nominal level at n = 12 and 33 for the three methods.

Files: `figures/K_calibration_parity.svg`, `figures/K_calibration_parity.pdf`, `figures/K_calibration_parity.png`

Sources (sha256):

- `closed_loop/biphasic_benchmark/final_bench_20261119/evaluation/predictions_at_references.csv` `abc45673b3c54dcc…`

## L_evidence_and_limits

Figure L. Evidence map: which claims are demonstrated, not supported, failed or untested, per method and specimen; generalization boundaries.

Files: `figures/L_evidence_and_limits.svg`, `figures/L_evidence_and_limits.pdf`, `figures/L_evidence_and_limits.png`

Sources (sha256):


## M_H1_surface_normals

Figure M (supplement). H1 analytic holdout, n = 33: median and 95th-percentile surface-normal angle errors per method with paired V3-vs-maximin statistics (descriptive); Chamfer distance computed but not reported as evidence because its discretization floor (1.04 % of R0) exceeds the reconstruction error.

Files: `figures/M_H1_surface_normals.svg`, `figures/M_H1_surface_normals.pdf`, `figures/M_H1_surface_normals.png`

Sources (sha256):

- `closed_loop/response_adaptive_v3_2026/holdout_h1_2026/population.json` `92aff54255eb70c1…`
- `closed_loop/response_adaptive_v3_2026/holdout_h1_2026/references_h1.npy` `a3b2b3713da84103…`

## N_TL1_cusp_transition_errors

Figure N. TL1 regional reconstruction errors at n = 33 and 64 (per-start points, median bars) for V2 BO+MES, angular maximin and random, on all area references, the ±equal-biaxial cusp region, mechanism-transition references and the remaining area references; paired maximin/V2 HL ratios with exact 96.1 % CIs and Holm p (over 3 comparisons) from the verified TL1 evaluation.

Files: `figures/N_TL1_cusp_transition_errors.svg`, `figures/N_TL1_cusp_transition_errors.pdf`, `figures/N_TL1_cusp_transition_errors.png`

Sources (sha256):

- `closed_loop/truss_lattice_benchmark_2026/final_comparison_2026/campaign3_20261120/evaluation/metrics.csv` `43d33ed8d2cfd5b7…`
- `closed_loop/truss_lattice_benchmark_2026/final_comparison_2026/campaign3_20261120/evaluation/paired.csv` `4a58e0eacebeea50…`
- `closed_loop/truss_lattice_benchmark_2026/final_comparison_2026/campaign3_20261120/evaluation/regions.json` `af276f85bc731745…`

## TL1V3_X01_convergence_all_methods

Figure TL1 V3 transfer: convergence of frozen V3, V3-Signed, V2, angular maximin and random (primary, max, cusp, transition errors; 9 matched starts). (Imported unchanged from the audited TL1 package; PNG rendered at 300 dpi from its PDF.)

Files: `figures/TL1V3_X01_convergence_all_methods.svg`, `figures/TL1V3_X01_convergence_all_methods.pdf`, `figures/TL1V3_X01_convergence_all_methods.png`

Sources (sha256):

- `closed_loop/truss_lattice_benchmark_2026/final_comparison_2026/v3transfer_20261120/figures/X01_convergence_all_methods.svg` `500a721924bb65c8…`
- `closed_loop/truss_lattice_benchmark_2026/final_comparison_2026/v3transfer_20261120/figures/X01_convergence_all_methods.pdf` `d60b696e03afdfca…`
- `closed_loop/truss_lattice_benchmark_2026/final_comparison_2026/v3transfer_20261120/figures/X01_convergence_all_methods.png` `26c04914a34e3ac4…`

## TL1V3_X04_acquisition_maps_before_after

Figure TL1 V3 transfer: GP variance, LOO risk and weighted acquisition before and after the real 33rd MOOSE experiment, frozen V3 and V3-Signed. (Imported unchanged from the audited TL1 package; PNG rendered at 300 dpi from its PDF.)

Files: `figures/TL1V3_X04_acquisition_maps_before_after.svg`, `figures/TL1V3_X04_acquisition_maps_before_after.pdf`, `figures/TL1V3_X04_acquisition_maps_before_after.png`

Sources (sha256):

- `closed_loop/truss_lattice_benchmark_2026/final_comparison_2026/v3transfer_20261120/figures/X04_acquisition_maps_before_after.svg` `8968e402d964d6c6…`
- `closed_loop/truss_lattice_benchmark_2026/final_comparison_2026/v3transfer_20261120/figures/X04_acquisition_maps_before_after.pdf` `46247afbb2dc02bd…`
- `closed_loop/truss_lattice_benchmark_2026/final_comparison_2026/v3transfer_20261120/figures/X04_acquisition_maps_before_after.png` `59a26491fa89cbc3…`

## TL1V3_X05_paired_ratios

Figure TL1 V3 transfer: paired Hodges–Lehmann ratios with exact 96.1 % CIs at n = 33 and 64. (Imported unchanged from the audited TL1 package; PNG rendered at 300 dpi from its PDF.)

Files: `figures/TL1V3_X05_paired_ratios.svg`, `figures/TL1V3_X05_paired_ratios.pdf`, `figures/TL1V3_X05_paired_ratios.png`

Sources (sha256):

- `closed_loop/truss_lattice_benchmark_2026/final_comparison_2026/v3transfer_20261120/figures/X05_paired_ratios.svg` `2ea50b1110a53f25…`
- `closed_loop/truss_lattice_benchmark_2026/final_comparison_2026/v3transfer_20261120/figures/X05_paired_ratios.pdf` `f6e9646972b27f31…`
- `closed_loop/truss_lattice_benchmark_2026/final_comparison_2026/v3transfer_20261120/figures/X05_paired_ratios.png` `cb48f2d3f76e5937…`

## TL1V3_X06_regional_cusp_transition_errors

Figure TL1 V3 transfer: cusp, transition and other regional errors for all five designs. (Imported unchanged from the audited TL1 package; PNG rendered at 300 dpi from its PDF.)

Files: `figures/TL1V3_X06_regional_cusp_transition_errors.svg`, `figures/TL1V3_X06_regional_cusp_transition_errors.pdf`, `figures/TL1V3_X06_regional_cusp_transition_errors.png`

Sources (sha256):

- `closed_loop/truss_lattice_benchmark_2026/final_comparison_2026/v3transfer_20261120/figures/X06_regional_cusp_transition_errors.svg` `e654a7606201a5c1…`
- `closed_loop/truss_lattice_benchmark_2026/final_comparison_2026/v3transfer_20261120/figures/X06_regional_cusp_transition_errors.pdf` `4737db9e744d6d36…`
- `closed_loop/truss_lattice_benchmark_2026/final_comparison_2026/v3transfer_20261120/figures/X06_regional_cusp_transition_errors.png` `48878eafc5ccd29d…`

## TL1V3_X08_frozen_vs_signed_risk

Figure TL1 V3 transfer: frozen (antipodally symmetric) vs signed LOO-risk weighting evaluated on the real GP states. (Imported unchanged from the audited TL1 package; PNG rendered at 300 dpi from its PDF.)

Files: `figures/TL1V3_X08_frozen_vs_signed_risk.svg`, `figures/TL1V3_X08_frozen_vs_signed_risk.pdf`, `figures/TL1V3_X08_frozen_vs_signed_risk.png`

Sources (sha256):

- `closed_loop/truss_lattice_benchmark_2026/final_comparison_2026/v3transfer_20261120/figures/X08_frozen_vs_signed_risk.svg` `b22a5c0714850d60…`
- `closed_loop/truss_lattice_benchmark_2026/final_comparison_2026/v3transfer_20261120/figures/X08_frozen_vs_signed_risk.pdf` `8b16d2ae7052c57f…`
- `closed_loop/truss_lattice_benchmark_2026/final_comparison_2026/v3transfer_20261120/figures/X08_frozen_vs_signed_risk.png` `0df322c2f6e79526…`

## TL1_C04_v2_acquisition_before_after

Figure TL1: V2-A20 acquisition (GP latent variance) on all 8,192 signed candidates before and after a real MOOSE experiment (start 0, test 33). (Imported unchanged from the audited TL1 package; PNG rendered at 300 dpi from its PDF.)

Files: `figures/TL1_C04_v2_acquisition_before_after.svg`, `figures/TL1_C04_v2_acquisition_before_after.pdf`, `figures/TL1_C04_v2_acquisition_before_after.png`

Sources (sha256):

- `closed_loop/truss_lattice_benchmark_2026/final_comparison_2026/campaign3_20261120/figures/C04_v2_acquisition_before_after.svg` `a26f06d8816da2e9…`
- `closed_loop/truss_lattice_benchmark_2026/final_comparison_2026/campaign3_20261120/figures/C04_v2_acquisition_before_after.pdf` `55a068ac1fd40db2…`
- `closed_loop/truss_lattice_benchmark_2026/final_comparison_2026/campaign3_20261120/figures/C04_v2_acquisition_before_after.png` `bdbd36f596b28a4b…`

## TL1_C05_selections

Figure TL1: actual sequential selections of V2, angular maximin and random (start 0 in order; all starts pooled). (Imported unchanged from the audited TL1 package; PNG rendered at 300 dpi from its PDF.)

Files: `figures/TL1_C05_selections.svg`, `figures/TL1_C05_selections.pdf`, `figures/TL1_C05_selections.png`

Sources (sha256):

- `closed_loop/truss_lattice_benchmark_2026/final_comparison_2026/campaign3_20261120/figures/C05_selections.svg` `87044434ecc420ff…`
- `closed_loop/truss_lattice_benchmark_2026/final_comparison_2026/campaign3_20261120/figures/C05_selections.pdf` `f572374b8b5f3e1e…`
- `closed_loop/truss_lattice_benchmark_2026/final_comparison_2026/campaign3_20261120/figures/C05_selections.png` `21f649000fca1951…`

## TL1_C06_reconstruction_sweep_eb_ps

Figure TL1: GP reconstructions (predictions, 95 % band) at n = 8, 16, 33, 64 along the equal-biaxial→shear great circle vs held-out MOOSE references, per method (start 0). (Imported unchanged from the audited TL1 package; PNG rendered at 300 dpi from its PDF.)

Files: `figures/TL1_C06_reconstruction_sweep_eb_ps.svg`, `figures/TL1_C06_reconstruction_sweep_eb_ps.pdf`, `figures/TL1_C06_reconstruction_sweep_eb_ps.png`

Sources (sha256):

- `closed_loop/truss_lattice_benchmark_2026/final_comparison_2026/campaign3_20261120/figures/C06_reconstruction_sweep_eb_ps.svg` `804063377a4ab971…`
- `closed_loop/truss_lattice_benchmark_2026/final_comparison_2026/campaign3_20261120/figures/C06_reconstruction_sweep_eb_ps.pdf` `7c696ea05d336790…`
- `closed_loop/truss_lattice_benchmark_2026/final_comparison_2026/campaign3_20261120/figures/C06_reconstruction_sweep_eb_ps.png` `06562b0780198c99…`

## TL1_C07_error_vs_n

Figure TL1: reconstruction error vs number of MOOSE experiments (primary area-weighted, maximum and cusp errors), median and spread over 9 starts. (Imported unchanged from the audited TL1 package; PNG rendered at 300 dpi from its PDF.)

Files: `figures/TL1_C07_error_vs_n.svg`, `figures/TL1_C07_error_vs_n.pdf`, `figures/TL1_C07_error_vs_n.png`

Sources (sha256):

- `closed_loop/truss_lattice_benchmark_2026/final_comparison_2026/campaign3_20261120/figures/C07_error_vs_n.svg` `f8ded69de5f204b7…`
- `closed_loop/truss_lattice_benchmark_2026/final_comparison_2026/campaign3_20261120/figures/C07_error_vs_n.pdf` `4e5fca362fc92ed5…`
- `closed_loop/truss_lattice_benchmark_2026/final_comparison_2026/campaign3_20261120/figures/C07_error_vs_n.png` `5c5fc8361560ff97…`

## TL1_C08_error_vs_cpu

Figure TL1: error vs cumulative CPU-hours (MOOSE + controller). (Imported unchanged from the audited TL1 package; PNG rendered at 300 dpi from its PDF.)

Files: `figures/TL1_C08_error_vs_cpu.svg`, `figures/TL1_C08_error_vs_cpu.pdf`, `figures/TL1_C08_error_vs_cpu.png`

Sources (sha256):

- `closed_loop/truss_lattice_benchmark_2026/final_comparison_2026/campaign3_20261120/figures/C08_error_vs_cpu.svg` `b54f14a638ace655…`
- `closed_loop/truss_lattice_benchmark_2026/final_comparison_2026/campaign3_20261120/figures/C08_error_vs_cpu.pdf` `b5f4eaac614e05bc…`
- `closed_loop/truss_lattice_benchmark_2026/final_comparison_2026/campaign3_20261120/figures/C08_error_vs_cpu.png` `e36167cd15407bc0…`

## TL1_C09_paired_ratios

Figure TL1: paired Hodges–Lehmann error ratios with exact 96.1 % CIs and Holm-adjusted p (over 3 comparisons), and per-start ratios. (Imported unchanged from the audited TL1 package; PNG rendered at 300 dpi from its PDF.)

Files: `figures/TL1_C09_paired_ratios.svg`, `figures/TL1_C09_paired_ratios.pdf`, `figures/TL1_C09_paired_ratios.png`

Sources (sha256):

- `closed_loop/truss_lattice_benchmark_2026/final_comparison_2026/campaign3_20261120/figures/C09_paired_ratios.svg` `ee4aa1c9c1e43be1…`
- `closed_loop/truss_lattice_benchmark_2026/final_comparison_2026/campaign3_20261120/figures/C09_paired_ratios.pdf` `fb17d4cbac11731a…`
- `closed_loop/truss_lattice_benchmark_2026/final_comparison_2026/campaign3_20261120/figures/C09_paired_ratios.png` `129c7dd126d083b1…`

## TL1_C10b_regional_errors

Figure TL1: regional errors (cusp, transition, non-feature) per method and budget. (Imported unchanged from the audited TL1 package; PNG rendered at 300 dpi from its PDF.)

Files: `figures/TL1_C10b_regional_errors.svg`, `figures/TL1_C10b_regional_errors.pdf`, `figures/TL1_C10b_regional_errors.png`

Sources (sha256):

- `closed_loop/truss_lattice_benchmark_2026/final_comparison_2026/campaign3_20261120/figures/C10b_regional_errors.svg` `93da79a5d4d5489b…`
- `closed_loop/truss_lattice_benchmark_2026/final_comparison_2026/campaign3_20261120/figures/C10b_regional_errors.pdf` `72b06214fc0410d2…`
- `closed_loop/truss_lattice_benchmark_2026/final_comparison_2026/campaign3_20261120/figures/C10b_regional_errors.png` `27669d6451b3e820…`

## TL1_C12_predicted_vs_measured

Figure TL1: GP-predicted vs MOOSE-measured yield radii at the references. (Imported unchanged from the audited TL1 package; PNG rendered at 300 dpi from its PDF.)

Files: `figures/TL1_C12_predicted_vs_measured.svg`, `figures/TL1_C12_predicted_vs_measured.pdf`, `figures/TL1_C12_predicted_vs_measured.png`

Sources (sha256):

- `closed_loop/truss_lattice_benchmark_2026/final_comparison_2026/campaign3_20261120/figures/C12_predicted_vs_measured.svg` `21384e66003b8ff3…`
- `closed_loop/truss_lattice_benchmark_2026/final_comparison_2026/campaign3_20261120/figures/C12_predicted_vs_measured.pdf` `0f5f5f3eeef71ac2…`
- `closed_loop/truss_lattice_benchmark_2026/final_comparison_2026/campaign3_20261120/figures/C12_predicted_vs_measured.png` `945a698a90f841ce…`

## TL1_P_F2_bar_fields_deformed

Figure TL1: deformed lattice (MOOSE nodal displacements, magnified), bar axial stress and plastic strain at six loadings. (Imported unchanged from the audited TL1 package; PNG rendered at 300 dpi from its PDF.)

Files: `figures/TL1_P_F2_bar_fields_deformed.svg`, `figures/TL1_P_F2_bar_fields_deformed.pdf`, `figures/TL1_P_F2_bar_fields_deformed.png`

Sources (sha256):

- `closed_loop/truss_lattice_benchmark_2026/final_comparison_2026/campaign3_20261120/figures/P_F2_bar_fields_deformed.svg` `e75ee0fc565fdc15…`
- `closed_loop/truss_lattice_benchmark_2026/final_comparison_2026/campaign3_20261120/figures/P_F2_bar_fields_deformed.pdf` `148bb411e6a90913…`
- `closed_loop/truss_lattice_benchmark_2026/final_comparison_2026/campaign3_20261120/figures/P_F2_bar_fields_deformed.png` `317ff4f59b9d447d…`

## Animations

`animations/B1_BO_s0_sequence.gif`, `animations/B2_BO_s0_sequence.gif` (+ per-step PNG frames n = 8–33): GP uncertainty over the loading-axis map and errors at the 72 MOOSE references after each actual MOOSE test of the BO start-0 campaigns.
