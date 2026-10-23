#!/usr/bin/env bash
# Validate a run configuration and samplesheet, and print the routing plan.
#
#   bash scripts/validate_config.sh config/config.yaml config/samples.tsv
#
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CFG="${1:-$ROOT/config/config.yaml}"
SHEET="${2:-$ROOT/config/samples.tsv}"

echo "[PfVarX] validating $CFG against $SHEET"
pixi run -e default pfvarx validate --config "$CFG" --samplesheet "$SHEET" --check-paths --json
