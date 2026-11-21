#!/usr/bin/env bash
# Mode B (no MOOSE): replay archived REAL MOOSE campaigns through the frozen controllers and verify every selection against the frozen records:
# B1 and B2 composites (V2), TL1 truss (V2, frozen V3, V3-Signed).  Outputs in test_outputs/replay_<stamp>/.
set -uo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; OUT="$ROOT/test_outputs/replay_$(date +%Y%m%dT%H%M%S)"
export YSBO_ROOT="$ROOT" PYTHONPATH="$ROOT/software" PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLBACKEND=Agg
PY="${YSBO_PYTHON:-}"; [ -z "$PY" ] && [ -x "$ROOT/.venv/bin/python" ] && PY="$ROOT/.venv/bin/python"; [ -z "$PY" ] && PY="$(command -v python3)"
"$PY" -B "$ROOT/software/tests/run_tests.py" --outdir "$OUT" --only T05,T06,T07,T08,T09,T10
