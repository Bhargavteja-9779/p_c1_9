# PROJECT_STATE

Last updated: 2026-11-17 (session continues across container restarts)

## TOPIC (LOCKED)
Standardization-invariance violations in molecular property prediction + SORB method. Target: JCIM (ACS).
Thesis: chemically-equivalent standardization choices barely change aggregate benchmark metrics but
materially change individual predictions (representation-dependent); SORB (orbit augmentation + inference
averaging) fixes it. See research_notes/research_plan.md (+ amendment).

## JOURNAL (SELECTED)
Journal of Chemical Information and Modeling (ACS). JIF ~6.4 (2025), Q1. Census median sub->acc 84 d (<=100).
journal_selection/*.md complete.

## STATUS
- Full experiment running (background task; relaunch `cd Target_Journal_Submission/code && python3 run_all.py --mode full --nproc 4` on restart; it RESUMES from results/_cells/ checkpoints + data_cache/std_cache/).
- Phases: sensitivity grid (RQ1/RQ2), ablation (RQ3), method+invariance (RQ4/RQ5), semi-synthetic shortcut (RQ5b).
- Grid: 55 datasets; sens/abl models [lgbm,rf,linear] x {ecfp,desc or ecfp} x 3 seeds; method/shortcut [lgbm,linear].
- Pilot results (consistent, documented in EXPERIMENT_LOG): aggregate spread ~0.5% median; orbit per-prediction
  SD 0.10 (tree) to 0.21-0.70 (linear-fp) of label SD; SORB train-aug reduces invariance 2-3x; test-orbit safe.

## CODE (complete, tested; 11 unit tests pass)
- src/preprocessing/standardize.py (orbit, pipelines, bounded tautomer enumerator, Dimorphite pH7.4)
- src/preprocessing/{featurize,cache}.py ; src/data/{datasets,splits}.py ; src/models/{zoo,sorb}.py
- src/training/{harness,invariance,shortcut,grid}.py ; src/stat_tests/compare.py
- src/visualization/{analyze,figures,tables}.py ; run_all.py ; scripts/leakage_audit.py ; tests/

## MANUSCRIPT
- manuscript/manuscript.md: full draft with {{TOKENS}} (23) for all numbers.
- manuscript/fill_manuscript.py substitutes tokens from results/metrics.json + tables + leakage_audit.
- references.bib (22 verified entries); DOCX pipeline (pandoc + assets/acs.csl + acs_template.docx) tested.
- cover_letter.md, code/LICENSE done.

## NEXT EXACT STEPS (when grid DONE)
1. python scripts/leakage_audit.py
2. python -m src.visualization.analyze ; python -m src.visualization.tables ; python -m src.visualization.figures
3. cd ../manuscript && python fill_manuscript.py  -> verify 0 tokens missing; sanity-check numbers vs pilot
4. VERIFY manuscript prose matches actual magnitudes (adjust hedges/claims to data); claim-evidence audit (Phase 22)
5. pandoc build Final_Manuscript.docx ; Phase 24-33 ten adversarial review rounds; Phase 34-42 audits + risk report
