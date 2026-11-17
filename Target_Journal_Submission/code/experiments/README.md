# Experiments
Experiment drivers live in `../scripts/` and are orchestrated by `../run_all.py`:
- `run_all.py --mode {quick,full}` — sensitivity (RQ1/RQ2), ablation (RQ3), SORB+invariance (RQ4/RQ5), shortcut (RQ5b).
- `scripts/revision_experiments.py` — cross-pipeline disagreement, canonical-tautomer fix, random-split replication.
- `scripts/leakage_audit.py` — data-integrity audit.
Analysis/figures/tables: `src/visualization/{analyze,figures,tables}.py`.
