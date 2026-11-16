# Preliminary premise check (2026-11-16, before full experiments)
Sampled up to 300 molecules from 11 datasets; applied the six named standardization pipelines.
- 59.2% of molecules yield a DIFFERENT canonical SMILES across pipelines (per-dataset 13%-86%).
- 9.9% yield a DIFFERENT ECFP4(2048) across pipelines (per-dataset 0%-31%).
- Ionizable/tautomerizable-rich datasets (ChEMBL potency, Caco2, PPBR, Lipophilicity) vary most;
  simple physchem sets (ESOL, FreeSolv) vary least -> supports RQ3 mechanism.
Script: research_notes/scripts/premise_probe.py (reproducible, seed 0).
This is a cheap descriptive probe, NOT a modeling result.
