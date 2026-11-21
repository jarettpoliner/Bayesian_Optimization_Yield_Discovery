#!/usr/bin/env bash
# Create the Python environment for Yield Surface BO inside this directory (.venv).  Nothing outside the release directory is modified.
#   ./setup.sh                 venv + pinned PyPI wheels (requirements-lock.txt; needs Python 3.10+, ideally 3.12 for exact reproduction)
#   ./setup.sh --conda         conda/micromamba environment from environment.yml (prefix ./.conda-env)
#   PYTHON=/path/to/python3.12 ./setup.sh
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [ "${1:-}" = "--conda" ]; then
  CM="$(command -v micromamba || command -v mamba || command -v conda || true)"; [ -z "$CM" ] && { echo "no conda/mamba/micromamba found"; exit 1; }
  "$CM" create -y -p "$ROOT/.conda-env" -f "$ROOT/environment.yml"; ln -sfn "$ROOT/.conda-env" "$ROOT/.venv"; echo "environment: $ROOT/.conda-env"; exit 0
fi
PY="${PYTHON:-}"; if [ -z "$PY" ]; then for c in python3.12 python3.11 python3.10 python3; do if command -v "$c" >/dev/null; then PY="$(command -v $c)"; break; fi; done; fi
[ -z "$PY" ] && { echo "Python 3.10+ not found; set PYTHON=/path/to/python3.12"; exit 1; }
"$PY" -c 'import sys; v=sys.version_info; assert v >= (3,10), "Python >= 3.10 required"; print("using Python", sys.version.split()[0], "(exact bitwise reproduction of the frozen records was verified with 3.12)" if v[:2]!=(3,12) else "")'
"$PY" -m venv "$ROOT/.venv"
"$ROOT/.venv/bin/python" -m pip install --upgrade pip >/dev/null
"$ROOT/.venv/bin/python" -m pip install -r "$ROOT/requirements-lock.txt"
echo "done: $ROOT/.venv  ->  next: ./verify_installation.sh"
