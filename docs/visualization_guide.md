# Visualization guide

- **Coordinates:** x = σxx, y = σyy, z = τxy (plane stress, σzz = 0), with σ = r(d̂) d̂ in Mandel form d̂ = (σxx, σyy, √2 τxy)/‖·‖.
- **Colour maps** keep the quantities apart: predictive SD of log r (plasma), actual error (red), acquisition score (viridis).
- **Format:** GIFs here are 960×540 copies of the 1280×720 atlas GIFs, with the original per-frame timing.
- **Provenance:** frames are recorded states only, with no interpolation. GP frames are exact replays of the recorded selection sequence through the frozen GP, with selections and score vectors verified bitwise against the records. Composite deformation is decoded from the original MOOSE checkpoints, and truss deformation comes from recorded node displacements.

## Animations

| File | What it shows | Frames | Length (s) | Size (MB) |
|---|---|---|---|---|
| [`gp_posterior/B1_composite_GP_reconstruction_vs_MOOSE_references.gif`](../animations/gp_posterior/B1_composite_GP_reconstruction_vs_MOOSE_references.gif) | B1 composite, V2 start 0, n = 8…33: GP posterior-mean surface, 72 independent MOOSE references, parity plot (±2 SD), actual error vs predictive SD. | 27 | 27.5 | 1.7 |
| [`gp_posterior/B1_composite_GP_predictive_uncertainty.gif`](../animations/gp_posterior/B1_composite_GP_predictive_uncertainty.gif) | Latent predictive SD of log r on the surface and on an equal-area direction map; uncertainty collapses near observations; ML-II starts at n = 20. | 27 | 27.5 | 2.0 |
| [`gp_posterior/B1_composite_GP_posterior_mean.gif`](../animations/gp_posterior/B1_composite_GP_posterior_mean.gif) | Posterior-mean radius and the Matérn hyper-parameters (ℓ, amplitude) at every observation. | 27 | 27.5 | 1.7 |
| [`gp_posterior/B1_composite_GP_actual_error.gif`](../animations/gp_posterior/B1_composite_GP_actual_error.gif) | Actual error |r̂/r − 1| at the 72 MOOSE references (known only there). | 27 | 27.5 | 1.9 |
| [`gp_posterior/TL1_truss_GP_reconstruction_vs_MOOSE_references.gif`](../animations/gp_posterior/TL1_truss_GP_reconstruction_vs_MOOSE_references.gif) | TL1 truss, V2 start 0, n = 8…64, against 482 MOOSE references; the GP smooths the truss surface's cusp-like corners. | 58 | 39.8 | 4.0 |
| [`gp_posterior/H1_Yld91_V3_GP_predictive_uncertainty.gif`](../animations/gp_posterior/H1_Yld91_V3_GP_predictive_uncertainty.gif) | H1 Yld91 realisation 17, V3, n = 8…64: predictive uncertainty. | 58 | 39.8 | 4.5 |
| [`gp_posterior/H1_Yld91_V3_GP_actual_error.gif`](../animations/gp_posterior/H1_Yld91_V3_GP_actual_error.gif) | Same campaign: actual error against the exact analytical surface (20,000 directions). | 58 | 39.8 | 4.6 |
| [`bayesian_optimization/B1_composite_V2_max_entropy_selection.gif`](../animations/bayesian_optimization/B1_composite_V2_max_entropy_selection.gif) | Two frames per BO step: (a) the acquisition score over 207 candidate axes and the selected direction with the controller's recorded prediction ±2 SD; (b) the recorded MOOSE radius and the updated GP. One recorded failed MOOSE test is shown as a failure. | 52 | 57.0 | 3.3 |
| [`bayesian_optimization/B2_composite_V2_max_entropy_selection.gif`](../animations/bayesian_optimization/B2_composite_V2_max_entropy_selection.gif) | Same for B2. | 50 | 55.0 | 3.0 |
| [`bayesian_optimization/H1_Yld91_V3_weighted_variance_selection.gif`](../animations/bayesian_optimization/H1_Yld91_V3_weighted_variance_selection.gif) | V3 (σ²·ρ) selection over 4,096 candidate axes on an analytical Yld91 surface. | 112 | 117.0 | 9.8 |
| [`bayesian_optimization/H1_Yld91_V3_vs_angular_maximin.gif`](../animations/bayesian_optimization/H1_Yld91_V3_vs_angular_maximin.gif) | V3 vs angular maximin from the same 8 initial directions; realisation 17 was chosen because its n = 33 ratio (1.26×) is closest to the 24-surface median. | 57 | 34.5 | 3.8 |
| [`bayesian_optimization/TL1_truss_V3_selection_method_demonstration.gif`](../animations/bayesian_optimization/TL1_truss_V3_selection_method_demonstration.gif) | V3 on the TL1 truss (8,192 signed candidates). Method demonstration only; no superiority claim. | 112 | 117.0 | 8.8 |
| [`yield_surface_3d/A6_1000_observation_BO_vs_exact_surface.gif`](../animations/yield_surface_3d/A6_1000_observation_BO_vs_exact_surface.gif) | A genuine recorded 1,000-step V2 trajectory on analytical surface A6 (oracle evaluations, not MOOSE). Frames: every n to 64, every 8th to 256, every 16th to 1,000. | 129 | 51.1 | 6.9 |
| [`yield_surface_3d/H1_Yld91_V3_reconstruction_vs_exact_surface.gif`](../animations/yield_surface_3d/H1_Yld91_V3_reconstruction_vs_exact_surface.gif) | V3 observations accumulating on the exact Yld91 surface, and the error-coloured GP reconstruction. | 58 | 35.2 | 3.3 |
| [`yield_surface_3d/B1_composite_yield_surface_turntable.gif`](../animations/yield_surface_3d/B1_composite_yield_surface_turntable.gif) | Final B1 GP reconstruction (n = 33) with 72 MOOSE references; the camera orbits. | 60 | 7.2 | 2.4 |
| [`yield_surface_3d/B2_composite_yield_surface_turntable.gif`](../animations/yield_surface_3d/B2_composite_yield_surface_turntable.gif) | Final B2 GP reconstruction (n = 33) with 71 MOOSE references. | 60 | 7.2 | 2.3 |
| [`yield_surface_3d/TL1_truss_yield_surface_turntable.gif`](../animations/yield_surface_3d/TL1_truss_yield_surface_turntable.gif) | Final TL1 GP reconstruction (n = 64) with 482 MOOSE references. | 60 | 7.2 | 2.5 |
| [`yield_surface_3d/H1_Yld91_yield_surface_turntable.gif`](../animations/yield_surface_3d/H1_Yld91_yield_surface_turntable.gif) | Exact Yld91 surface (realisation 17) with the 64 V3 observations. | 60 | 7.2 | 2.2 |
| [`microstructure/B1_composite_tension_dominated_deformation.gif`](../animations/microstructure/B1_composite_tension_dominated_deformation.gif) | Recorded MOOSE ramp steps decoded from the original checkpoints (no MOOSE re-run): total displacement ×25, von Mises, max principal stress, eqps; direction d̂ = (0.104, 0.442, 0.891). | 11 | 15.9 | 1.3 |
| [`microstructure/B1_composite_compression_dominated_deformation.gif`](../animations/microstructure/B1_composite_compression_dominated_deformation.gif) | Same camera; d̂ = (−0.923, 0.154, 0.352). | 8 | 12.6 | 1.0 |
| [`microstructure/B2_composite_tension_dominated_deformation.gif`](../animations/microstructure/B2_composite_tension_dominated_deformation.gif) | B2 (weak inclusion), tension-dominated direction. | 12 | 17.0 | 1.5 |
| [`microstructure/B2_composite_compression_dominated_deformation.gif`](../animations/microstructure/B2_composite_compression_dominated_deformation.gif) | B2, compression-dominated direction. | 9 | 13.7 | 1.1 |
| [`microstructure/B1_composite_uniaxial_tension_compression_stress_evolution.gif`](../animations/microstructure/B1_composite_uniaxial_tension_compression_stress_evolution.gif) | Uniaxial ±x: no displacement history was recorded, so this shows the recorded phase-averaged stresses, plastic volume fractions and the yield-state σxx field. | 9 | 14.5 | 0.6 |
| [`microstructure/B2_composite_uniaxial_tension_compression_stress_evolution.gif`](../animations/microstructure/B2_composite_uniaxial_tension_compression_stress_evolution.gif) | Same for B2: the weak inclusion saturates at its 90 MPa yield stress. | 7 | 12.5 | 0.5 |
| [`truss/TL1_truss_uniaxial_tension_x.gif`](../animations/truss/TL1_truss_uniaxial_tension_x.gif) | 3×3 periodic cells from the recorded node displacements (×110), bar axial stress and plastic strain; K0 infill yields. | 64 | 30.4 | 3.7 |
| [`truss/TL1_truss_uniaxial_compression_x.gif`](../animations/truss/TL1_truss_uniaxial_compression_x.gif) | Uniaxial compression −x. | 64 | 30.4 | 3.6 |
| [`truss/TL1_truss_uniaxial_tension_y.gif`](../animations/truss/TL1_truss_uniaxial_tension_y.gif) | Uniaxial tension +y; K60/K120 yield. | 41 | 22.4 | 2.6 |
| [`truss/TL1_truss_uniaxial_compression_y.gif`](../animations/truss/TL1_truss_uniaxial_compression_y.gif) | Uniaxial compression −y. | 41 | 22.4 | 2.4 |

## Figures

| Folder | Content |
|---|---|
| `figures/yield_surface_3d/` | 3-D surfaces with exact GP sections (σxx–σyy and σxx–τxy planes, ±2 SD bands) for B1, B2, TL1 and H1, plus an overview. In the SVGs the 3-D panel is an embedded raster. |
| `figures/results/` | Physical final GP vs MOOSE references (parity, error histograms, coverage); H1 V3 vs maximin at n = 33 and n = 64; A6 at n = 1,000 |
| `figures/physical/` | Final-step composite and truss fields; uniaxial composite yield states |
| `figures/bayesian_optimization/` | Example acquisition step |
| `figures/presentation/` | The verified final-presentation figures (PNG + SVG), with captions in `CAPTIONS_from_final_presentation.md` |
| `figures/animations_contact_sheet.png` | Thumbnails of all atlas animations |
