#!/usr/bin/env bash
# Mode A (no MOOSE): V2-A20 BO+MES, V3 weighted and angular maximin on a frozen H1 Yld91 analytic surface; prints the error curves.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
for m in v2_mes v3_weighted angular_maximin; do
  rm -rf "$ROOT/outputs/analytic_$m"; "$ROOT/ysbo" run "$ROOT/examples/analytic_$m.yaml"
  echo "   error curve:"; sed 's/^/     /' "$ROOT/outputs/analytic_$m/error_curve.csv"
done
echo "outputs: $ROOT/outputs/analytic_*  (selections.csv, reconstruction_grid.csv, reconstruction.png, error_curve.csv, steps.csv)"
