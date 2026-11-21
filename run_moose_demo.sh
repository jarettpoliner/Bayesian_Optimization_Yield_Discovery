#!/usr/bin/env bash
# Mode C (LIVE MOOSE; requires a compatible MOOSE build, see INSTALL.md).
#   ./run_moose_demo.sh                 detect + validate the executable, generate a TL1 input and run `--check-input` (no simulation)
#   ./run_moose_demo.sh --live          run examples/moose_tl1_v2.yaml (V2-A20 on TL1, 33 MOOSE runs, CPU cap 0.5 h; ~5 min on one core)
# Executable: export YSBO_MOOSE_EXE=/path/to/combined-opt   (optional YSBO_MOOSE_ENV_PREFIX=<conda env of the MOOSE build>)
# Stop safely: touch outputs/moose_tl1_v2/STOP      Resume: run the same command again
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
[ -z "${YSBO_MOOSE_EXE:-}" ] && { echo "set YSBO_MOOSE_EXE to a MOOSE executable with the solid_mechanics module (e.g. modules/combined/combined-opt)"; exit 2; }
mkdir -p "$ROOT/outputs"; "$ROOT/ysbo" moose-check --exe "$YSBO_MOOSE_EXE" --check-input --workspace "$ROOT/outputs/moose_check_workspace"
if [ "${1:-}" = "--live" ]; then "$ROOT/ysbo" run "$ROOT/examples/moose_tl1_v2.yaml"; else echo "dry run only; add --live to start MOOSE experiments"; fi
