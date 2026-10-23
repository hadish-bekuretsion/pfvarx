#!/usr/bin/env bash
# Build the distributable archive for PfVarX.
#
# The archive contains the package, workflow, configuration, container
# definitions, tests, examples, benchmarks and documentation.  Caches, build
# artefacts and the rendered PDF preview are excluded; the PDF is a
# convenience copy and is regenerated from the DOCX by scripts/build_methodology_docx.sh.
set -euo pipefail

cd "$(dirname "$0")/.."
NAME="pfvarx-$(python -c 'import pfvarx; print(pfvarx.__version__)')"
OUT="$(cd .. && pwd)/${NAME}.zip"

rm -f "$OUT"
cd ..

# Archive a single top-level directory so that extraction does not scatter
# files into the caller's working directory.
# The version-control directory is excluded: a release archive carries the
# working tree, not the history.  History is transferred separately with
# `git bundle`, which is the supported way to move a repository between
# machines.
zip -r -q "$OUT" "pfvarx" \
    -x 'pfvarx/.git/*' \
    -x '*/__pycache__/*' \
    -x '*.pyc' \
    -x '*/.pytest_cache/*' \
    -x '*/.ruff_cache/*' \
    -x '*.egg-info/*' \
    -x '*/build/*' \
    -x '*/dist/*' \
    -x '*.docx~' \
    -x 'pfvarx/todo.md' \
    -x '*/docs/PfVarX_Methodology.pdf'

echo "wrote $OUT"
unzip -l "$OUT" | tail -3
