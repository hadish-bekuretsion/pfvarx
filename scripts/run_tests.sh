#!/usr/bin/env bash
# Run the full PfVarX test suite and lint through the project-managed Pixi env.
#
#   bash scripts/run_tests.sh
#
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

pixi install --locked -e default

echo "[PfVarX] lint"
pixi run -e default ruff check src tests

echo "[PfVarX] tests"
pixi run -e default test

echo "[PfVarX] release status"
pixi run -e default status | sed -n '1,2p'
