#!/usr/bin/env bash
# Verify the installation from THIS directory only: Python dependencies, frozen-source hashes and the full test suite
# (analytic smoke tests, frozen H1 replay, recorded MOOSE replays B1/B2/TL1, checkpoint/resume, library and symmetry checks).
# HOME is redirected to an empty temporary directory so nothing on the original workstation can be used.  No MOOSE run.
#   ./verify_installation.sh [--quick] [--manifest] [--with-moose]
set -uo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; STAMP="$(date +%Y%m%dT%H%M%S)"; OUT="$ROOT/test_outputs/verify_$STAMP"; mkdir -p "$OUT"
ARGS=(); MAN=0; for a in "$@"; do case "$a" in --quick) ARGS+=(--quick);; --with-moose) ARGS+=(--with-moose);; --manifest) MAN=1;; esac; done
if [ "$MAN" = 1 ]; then echo "== verifying MANIFEST.sha256 (all files; may take minutes)"; (cd "$ROOT" && sha256sum --quiet -c MANIFEST.sha256) && echo "MANIFEST OK" || { echo "MANIFEST FAILED"; exit 1; }; fi
FAKEHOME="$(mktemp -d)"; export MPLCONFIGDIR="$FAKEHOME/.mpl"
echo "== HOME redirected to $FAKEHOME"
HOME="$FAKEHOME" "$ROOT/ysbo" verify-frozen > "$OUT/frozen.json" && echo "frozen sources OK" || { echo "frozen sources FAILED (see $OUT/frozen.json)"; exit 1; }
export YSBO_ROOT="$ROOT" PYTHONPATH="$ROOT/software" PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLBACKEND=Agg
PY="${YSBO_PYTHON:-}"; [ -z "$PY" ] && [ -x "$ROOT/.venv/bin/python" ] && PY="$ROOT/.venv/bin/python"; [ -z "$PY" ] && PY="$(command -v python3)"
HOME="$FAKEHOME" "$PY" -B "$ROOT/software/tests/run_tests.py" --outdir "$OUT" "${ARGS[@]}"; RC=$?
echo "report: $OUT/tests_report.json"; rm -rf "$FAKEHOME"; exit $RC
