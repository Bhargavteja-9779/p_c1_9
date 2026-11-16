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
- Collect wf_364cfa10-0ca output (JCIM verification, benchmark papers, triaged topic candidates); write journal_selection/*.md; stress-test + novelty-defend chosen topic; implement.

## Recovery instructions
- All research outputs saved under research_notes/ (raw JSON from workflows) and Target_Journal_Submission/.
- Resume from the "Next exact step" above; do not re-run completed experiments (see EXPERIMENT_LOG.md).
