#!/usr/bin/env bash
# Iterated build of docs/PfVarX_Methodology.docx.
#
# A static, page-accurate contents table is used rather than a Word TOC field,
# because field codes are not populated by non-Word renderers and would leave
# the contents page blank.  Because inserting the contents table can shift
# pagination, the build is iterated until the recorded page numbers match the
# rendered document.
set -euo pipefail

cd "$(dirname "$0")/.."

OUT=docs/PfVarX_Methodology.docx
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

TITLES_JSON="$TMP/titles.json"
cat > "$TITLES_JSON" <<'JSON'
["1. Objectives, scope and methodological principles",
 "2. Architecture and workflow",
 "3. Input validation and mode-aware routing",
 "4. Preprocessing, alignment and read-level quality control",
 "5. Candidate discovery and evidence extraction",
 "6. Statistical models and formulations",
 "7. Model integration, calling, filtering and uncertainty",
 "8. Calibration and threshold justification",
 "9. Annotation and downstream analyses",
 "10. Unified quality control, reporting and dashboard",
 "11. Tools, dependencies, configuration and deployment",
 "12. Testing, benchmarking, validation, metrics and acceptance criteria",
 "13. Reproducibility, error handling, limitations and troubleshooting",
 "14. References and software citations",
 "15. Validation-status table"]
JSON

render() {  # $1 = docx path, $2 = outdir
    soffice --headless --convert-to pdf --outdir "$2" "$1" >/dev/null 2>&1
}

locate() {  # $1 = pdf, $2 = out json
    python - "$1" "$2" "$TITLES_JSON" <<'PY'
import json, re, subprocess, sys

pdf, out, titles_path = sys.argv[1], sys.argv[2], sys.argv[3]
titles = json.load(open(titles_path))

info = subprocess.run(["pdfinfo", pdf], capture_output=True, text=True).stdout
n_pages = int(re.search(r"Pages:\s+(\d+)", info).group(1))

pages = {}
for p in range(1, n_pages + 1):
    txt = subprocess.run(
        ["pdftotext", "-f", str(p), "-l", str(p), "-layout", pdf, "-"],
        capture_output=True, text=True).stdout
    flat = " ".join(txt.split())
    # The contents page lists every section title; skip it so that the table
    # of contents does not match itself.
    if "Section Page" in flat or flat.startswith("Contents"):
        continue
    for t in titles:
        if t not in pages and t in flat:
            pages[t] = p

missing = [t for t in titles if t not in pages]
if missing:
    sys.exit("could not locate: " + "; ".join(missing))
json.dump(pages, open(out, "w"), indent=1)
print(f"      {len(pages)} sections across {n_pages} pages")
PY
}

echo "[pass 1] build without contents"
python scripts/build_methodology_docx.py --out "$TMP/pass1.docx" >/dev/null
render "$TMP/pass1.docx" "$TMP"
locate "$TMP/pass1.pdf" "$TMP/toc.json"

for i in 2 3 4 5; do
    echo "[pass $i] build with contents"
    python scripts/build_methodology_docx.py --toc "$TMP/toc.json" --out "$TMP/pass$i.docx" >/dev/null
    render "$TMP/pass$i.docx" "$TMP"
    locate "$TMP/pass$i.pdf" "$TMP/toc_new.json"

    if cmp -s "$TMP/toc.json" "$TMP/toc_new.json"; then
        echo "      page numbers stable"
        cp "$TMP/pass$i.docx" "$OUT"
        echo "done: $OUT"
        exit 0
    fi
    echo "      pagination shifted; re-iterating"
    mv "$TMP/toc_new.json" "$TMP/toc.json"
done

echo "ERROR: contents page numbers did not stabilise" >&2
exit 1
