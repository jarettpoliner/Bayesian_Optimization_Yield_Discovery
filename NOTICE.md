# Notice: licensing scope, attribution and provenance

## 1. What the MIT License covers

The [`LICENSE`](LICENSE) (MIT, © 2026 Jarett Poliner) applies to the author's **original** work in this repository:
- the Bayesian-optimization software (controllers, Matérn‑3/2 GP, acquisition functions, candidate libraries, oracles, MOOSE adapter, public API, tests, examples);
- the MOOSE input decks and meshes written for this project;
- the recorded results and benchmark summaries;
- the figures, animations and presentation;
- the documentation.

The author has authorized this public redistribution (2026-11-21).

## 2. Earlier notices in archived files

Some frozen files predate this licensing decision and are kept **byte-for-byte unchanged**, so that their recorded SHA-256 hashes remain verifiable. In particular, `research/yield_surface_discovery_release/software/yield_surface_discovery/LICENSE_OR_NOTICE.md` says "All rights reserved pending the choice of an open-source licence" and "Contact the authors before redistribution". **For this repository, those statements are superseded by the MIT License above**, except for the third-party material listed in section 3, to which they never granted rights.

## 3. Third-party material: not owned, not relicensed

| Material | Where | Status |
|---|---|---|
| Derived numerical values (published fixed designs, error metrics, Cazacu–Barlat 2004 Mg–4%Li benchmark parameters) from **Fuhg, Fau, Bouklas & Marino (2023), *Eur. J. Mech. A/Solids* 99:104925**, and their public repository `FuhgJan/modelDataDriven_YieldFunction` (commit a4fb9ae) | `research/yield_surface_discovery_release/software/yield_surface_discovery/data/frozen_results/{fuhg_same_points.csv, runs.csv, PROVENANCE.json}`, `.../src/yield_surface_discovery/{analytical.py, figures.py, tables.py}`, `research/yield_surface_discovery_release/tables/{fuhg_same_points.*, main_results.csv}`, and the release reports that describe them | Attributed to the original authors. **No ownership is claimed, and they are not relicensed under MIT.** No code or model files from that work are included. Cite the original publication when using these values. |
| Published yield-function models implemented from the literature: Barlat, Lege & Brem (1991) Yld91; Hill (1948) quadratic; von Mises; Cazacu & Barlat (2004) | analytical oracles (`surfaces_h1.py`, `cs_common.py`, V1 `analytical.py`) | The implementations are original code. The models themselves belong to their authors and should be cited. |

## 4. External software (not included)

| Software | Status |
|---|---|
| **MOOSE** (idaholab/moose) | Not included. Install it separately under its own licence. Only this project's input decks and the adapter code that calls an external MOOSE executable are distributed here. |
| Python dependencies (numpy, scipy, pandas, scikit-learn, matplotlib, pyyaml) | Installed from PyPI under their own licences (`requirements-lock.txt`) |

## 5. Not included

- The unpublished manuscript on conditional-diffusion design of experiments, which has co-authors, and its LaTeX source.
- The full archival research tree, raw MOOSE checkpoints and trial directories (except the minimal recorded records needed by the replay tests), MP4s, ParaView data and full-resolution atlas media.

## 6. Provenance

- **Frozen modules:** the 43 frozen modules used by `software/yield_surface_bo` are verified by SHA-256 before import (`software/yield_surface_bo/FROZEN_SOURCES.json`, `./ysbo verify-frozen`).
- **Workstation paths:** a few frozen records and MOOSE input decks contain absolute paths from the original workstation. These are provenance only; no code here depends on them.
