# PROJECT_STATE

Last updated: 2026-11-15 (session 1)

## Status
- Stage: Phase 0 (journal selection) + Phase 1 (topic landscape) running as background workflows.
- Hardware: 4 CPU cores, 15 GB RAM, no GPU (CPU-only experiments).
- Repo branch: claude/busy-knuth-2w5kpv

## Decisions so far
- Target journal: NOT YET SELECTED
- Topic: NOT YET SELECTED

## Next exact step
- Collect workflow outputs (journal-verify A/B, topic-landscape), write journal_selection/*.md, select journal, run novelty-defense on top topic candidates.

## Recovery instructions
- All research outputs saved under research_notes/ (raw JSON from workflows) and Target_Journal_Submission/.
- Resume from the "Next exact step" above; do not re-run completed experiments (see EXPERIMENT_LOG.md).
