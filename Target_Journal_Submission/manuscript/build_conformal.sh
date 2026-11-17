#!/usr/bin/env bash
# Build the conformal paper: fill tokens -> resolve crossrefs/figures -> DOCX (+PDF).
# Requires: conformal_metrics.json (analyze.py) and figures (figures.py) already produced.
set -euo pipefail
cd "$(dirname "$0")"

echo "[1/5] fill tokens from conformal_metrics.json"
python3 fill_conformal.py

echo "[2/5] resolve crossrefs + embed figures"
python3 build_conformal.py

echo "[3/5] pandoc -> DOCX (ACS style)"
pandoc conformal_manuscript_build.md \
  --reference-doc=assets/acs_template.docx \
  --citeproc --csl=assets/acs.csl --bibliography=references.bib \
  --resource-path=.:../code/figures \
  -o Conformal_Manuscript.docx

echo "[4/5] DOCX -> PDF (LibreOffice)"
if command -v libreoffice >/dev/null 2>&1; then
  libreoffice --headless --convert-to pdf --outdir . Conformal_Manuscript.docx >/dev/null 2>&1 || echo "  (PDF conversion skipped)"
else
  echo "  (libreoffice unavailable; skipping PDF)"
fi

echo "[5/5] done"
ls -la Conformal_Manuscript.* 2>/dev/null || true
