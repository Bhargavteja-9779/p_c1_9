# PROJECT_STATE

Last updated: 2026-11-15 (session 1)

## Status
- Session 1 hit the 5-hour usage limit mid-workflow (2026-11-15/16); 3 of 28 research agents completed. Resumed 2026-11-16.
- Journal timeline evidence now gathered by script census (Crossref assertions + PubMed history) -> Target_Journal_Submission/journal_selection/timeline_census_table.md
- Leading journal candidate: Journal of Chemical Information and Modeling (ACS): census median received->accepted 84 d (n=347, 64% <=100 d). Verification + 20-paper benchmark + topic scan running as workflow wf_364cfa10-0ca.
- Data cached (gitignored): code/data_cache/admet_group (TDC ADMET group, 22 datasets, Dataverse file 4426004), code/data_cache/moleculeace_repo (MoleculeACE, 30 datasets).
- Hardware: 4 CPU cores, 15 GB RAM, no GPU (CPU-only experiments).
- Repo branch: claude/busy-knuth-2w5kpv

## Decisions so far
- Target journal: NOT YET SELECTED
- Topic: NOT YET SELECTED

## Next exact step
- Await broad-novelty-triage (wf_4928b6ce-c2a, 8 lower-pre-emption directions). Pick winner by: SURVIVES novelty + method-with-ablations + cannot-collapse-to-null + CPU + no heavy curation + JCIM-fit>=4 + low/medium pre-emption. Then Phase 5 (formulate contribution), implement, run experiments.

## Abandoned (Phase-4 novelty defense, see research_notes/preemption_findings.json)
- C1 conformal-under-censored-calibration: pre-empted by Sesia & Svetnik arXiv:2607.20736 (Jul 2026).
- D MoleculeACE pairwise-cliff/leakage: pre-empted by Cuccarese arXiv:2604.07560, ACNet, GraphCliff.
- delta assay-heterogeneity combination: pre-empted by Chan-Verdonk-Poelking 2308.09086, AssayMatch 2511.16087.
- P2 class-conditional conformal imbalance: pre-empted by Tursunbadalov arXiv:2607.06605.
- Lesson: molecular UQ/conformal/cliff/assay spaces are saturated 2025-26; target lower-pre-emption corners.

## Journal (leading, to finalize)
- Journal of Chemical Information and Modeling (ACS). JIF ~6.4 (2025, per aggregators; verify), Q1. ISSN 1549-9596/1549-960X.
- Census median received->accepted 84 days (n=347), 64% <=100 d -> satisfies <=100-day target.
- Mandatory Data & Software Availability statement (ACS editorial 10.1021/acs.jcim.0c01389). Exact author guidelines to confirm from accessible mirror at manuscript stage (ACS pages block automated fetch).
- Scope fit CONFIRMED: JCIM published many 2026 papers on conformal prediction, UQ, applicability domain, imbalanced regression, evaluation methodology, activity cliffs.

## Code foundation (done, tested)
- code/src/preprocessing/featurize.py (ECFP + 25 RDKit descriptors + standardization + Tanimoto)
- code/src/data/datasets.py (TDC 22, MoleculeACE 30, MoleculeNet 3 = up to 55 datasets)
- code/src/data/splits.py (random, Bemis-Murcko scaffold, repeated, three-way)
- code/src/evaluation/metrics.py (regression, classification, calibration/ECE, interval coverage, error-ranking)
- code/src/stat_tests/compare.py (paired Wilcoxon+rank-biserial, bootstrap CIs, Cohen dz, Friedman+Nemenyi CD, Holm-Bonferroni)
- code/src/utils/seeds.py (deterministic seeding, device detection)
