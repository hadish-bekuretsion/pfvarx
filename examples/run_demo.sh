#!/usr/bin/env bash
# PfVarX worked example: synthetic data -> call -> QC -> report.
#
# Uses the Pixi-managed Python core; no external bioinformatics binaries are needed.
#
#   bash examples/run_demo.sh [output-directory]
#
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="${1:-$ROOT/examples/output}"
CFG="$ROOT/config/profiles/wgs.yaml"
cd "$ROOT"

pixi install --locked -e default

echo "[PfVarX] version:"; pixi run -e default pfvarx --version
echo "[PfVarX] release status:"; pixi run -e default status | head -n 2

echo "[PfVarX] generating synthetic data and calling variants -> $OUT"
pixi run -e default pfvarx demo --config "$CFG" --out "$OUT" --coverage 60 --seed 11

echo "[PfVarX] unified QC"
pixi run -e default pfvarx qc --config "$CFG" --results "$OUT/results" --out "$OUT/results/run_qc.json"

echo "[PfVarX] exporting reports"
pixi run -e default pfvarx report all --config "$CFG" --results "$OUT/results" --out "$OUT/results/report" || \
    echo "[PfVarX] report export failed; inspect the error and installed report dependencies."

echo "[PfVarX] done. Outputs under: $OUT"
find "$OUT" -maxdepth 3 -type f | sort
