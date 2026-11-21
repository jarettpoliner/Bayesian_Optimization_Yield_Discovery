# Reproduction

## Without MOOSE

```bash
./setup.sh                    # ./.venv from requirements-lock.txt (PyPI); Python 3.12 for bitwise reproduction
./verify_installation.sh      # 16 tests
./run_analytic_demo.sh
./run_replay_demo.sh          # replays the recorded B1/B2/TL1 MOOSE campaigns through the frozen controllers
```

`verify_installation.sh` covers:
- imports and frozen-source SHA-256;
- library and symmetry checks;
- V2/V3/maximin analytic smoke tests;
- the frozen H1 V3 regression;
- bitwise replays of recorded B1/B2 V2 and TL1 V2/V3/V3-Signed campaigns;
- checkpoint resume;
- a truss input dry run;
- MOOSE-adapter detection;
- external-command plumbing;
- independence from the original workstation.

**Analytical campaigns:** `./ysbo run examples/analytic_{v2_mes,v3_weighted,angular_maximin}.yaml`. To run on an H1 surface, change `specimen.analytic.realisation` (0–23) and `start` (2000 + i).

**Frozen research runners.** These are included as code; the full studies are long. Some of them resolve paths below `~/projects/...` from the original workstation:
- H1: `research/yield_surface_moose/plane_stress/closed_loop/response_adaptive_v3_2026/holdout_h1_2026/run_h1.py`, then `evaluate_h1.py` and `analyze_h1.py`;
- A1–A6: `.../complexity_scaling_2026/run_analytic*.py`, then `evaluate_analytic_r3.py`.

The recorded inputs for their full re-evaluation are not part of this repository; only the results are.

## With MOOSE (installed separately; not included)

**Archived build:**
- idaholab/moose commit `08315af82ed0b71e3edf82c04161b2da6a1aff19` (2026-11-13);
- conda `moose-dev 2026.11.09-mpich`, `moose-libmesh 2026.10.21`, `moose-petsc 3.26.0`, `mpich 5.0.2`;
- application `modules/combined/combined-opt`.

Other versions are untested.

```bash
mamba create -n moose moose-dev && mamba activate moose
git clone https://github.com/idaholab/moose.git && cd moose && git checkout 08315af82ed0b71e3edf82c04161b2da6a1aff19
cd modules/combined && make -j 8
export YSBO_MOOSE_EXE=$PWD/combined-opt
./ysbo moose-check --exe "$YSBO_MOOSE_EXE" --check-input     # no simulation
./run_moose_demo.sh --live                                    # 33 TL1 truss experiments with V2 (CPU cap 0.5 h)
```

- **Live TL1 use** was validated only up to `--check-input` on the original build; a live run on a fresh system has not been executed.
- **B1/B2 composite campaigns** are not wrapped as a live adapter. The frozen scripts are `closed_loop/biphasic_benchmark/final_bench_20261119/run_bo.py` → `bench_oracle.py` → `oracle/biphasic_worker_r2.py`, and they use workstation paths. Other MOOSE models can be connected through `specimen.kind: external_command`.
