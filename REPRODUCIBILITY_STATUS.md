# REPRODUCIBILITY_STATUS

- One-command reproduction: `cd Target_Journal_Submission/code && python run_all.py --mode full` (CPU-only, ~2-3 h, checkpointed/resumable), then `python -m src.visualization.{analyze,tables,figures}` and `python scripts/leakage_audit.py` and `python scripts/revision_experiments.py`.
- Determinism: all randomness seeded (seeds 0-2); results cached; figures/tables regenerated only from machine-readable result files.
- Environment: requirements.txt pinned (Python 3.13, numpy 2.5, scikit-learn 1.9, lightgbm 4.7, rdkit 2026.09, dimorphite-dl, torch 2.15 CPU optional). environment.yml provided.
- Data: all public, download instructions in code/README.md (TDC Dataverse 4426004, MoleculeACE GitHub, MoleculeNet S3).
- Tests: 11 unit tests pass (`pytest tests/`).
- Manuscript numbers: 45 tokens auto-filled from results/metrics.json via manuscript/fill_manuscript.py; values.json is the audit trail; final consistency check passed (spot-checks match).
- Hardware used: 4 CPU cores, 15 GB RAM, no GPU.
- Status: REPRODUCIBLE.
