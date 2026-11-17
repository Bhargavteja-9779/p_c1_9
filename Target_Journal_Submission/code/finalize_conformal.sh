#!/usr/bin/env bash
# One-shot downstream pipeline for the conformal paper, run after the full study + alpha sweep finish.
# analyze -> figures -> (manuscript) fill -> crossref/figures -> DOCX/PDF -> SI
set -euo pipefail
CODE="$(cd "$(dirname "$0")" && pwd)"
MAN="$CODE/../manuscript"

echo "[analyze] primary metrics + tables"
python3 "$CODE/src/conformal/analyze.py" | tail -2

echo "[figures]"
python3 "$CODE/src/conformal/figures.py"

echo "[manuscript] fill tokens"
cd "$MAN" && python3 fill_conformal.py

echo "[manuscript] crossrefs + embed figures"
python3 build_conformal.py

echo "[SI] generate supporting information"
python3 si_conformal.py | tail -1

echo "[DOCX] pandoc"
pandoc conformal_manuscript_build.md \
  --reference-doc=assets/acs_template.docx \
  --citeproc --csl=assets/acs.csl --bibliography=references.bib \
  --resource-path=.:../code/figures \
  -o Conformal_Manuscript.docx
pandoc conformal_SI.md \
  --reference-doc=assets/acs_template.docx \
  --citeproc --csl=assets/acs.csl --bibliography=references.bib \
  -o Conformal_SI.docx

echo "[PDF] libreoffice"
libreoffice --headless --convert-to pdf --outdir . Conformal_Manuscript.docx >/dev/null 2>&1 || echo "  (manuscript PDF skipped)"
libreoffice --headless --convert-to pdf --outdir . Conformal_SI.docx >/dev/null 2>&1 || echo "  (SI PDF skipped)"

echo "[done]"
ls -la Conformal_Manuscript.* Conformal_SI.* 2>/dev/null || true
