# Submission package — conformal paper (PRIMARY deliverable)

**Title.** A False Sense of Security: Standard Conformal Prediction Under-Covers Novel Chemotypes, and
an Applicability-Domain-Conditional Fix

**Target journal.** Journal of Chemical Information and Modeling (ACS), Q1. See FINAL_DECISION_conformal.md.

**Headline result (scaffold split, 41 datasets × 3 models × 2 featurizations × 3 seeds).**
Standard conformal meets marginal coverage (0.895) but under-covers the most novel stratum
(far-AD coverage 0.820 vs 0.90). Applicability-domain-conditional (Mondrian-AD) conformal raises it to
0.901 (+8.0 points; paired Wilcoxon p = 6×10⁻¹²; improved on 95% of datasets; Friedman over 6 methods
p = 8×10⁻¹⁷). Learned-difficulty normalized (0.840) and scaffold-Mondrian (0.828) do not fix it;
distance-normalized conformal (0.892) helps, confirming the applicability-domain axis is the right one.
Average width cost 4%, concentrated in the novel strata (far-stratum width ~doubles). Every number is
token-substituted from code/results/conformal_metrics.json — no hand-typed values.

## Files to submit
| Item | File |
|---|---|
| Manuscript (ACS Word) | manuscript/Conformal_Manuscript.docx |
| Manuscript (PDF) | manuscript/Conformal_Manuscript.pdf |
| Supporting Information (Word) | manuscript/Conformal_SI.docx |
| Supporting Information (PDF) | manuscript/Conformal_SI.pdf |
| Cover letter (Word) | manuscript/Conformal_Cover_Letter.docx |
| Figures (600 dpi, vector PDF) | code/figures/cfig_{nearfar,profile,tradeoff,perdataset,selective}.pdf |

## Source of truth (reproducibility)
- Manuscript source (token-templated): manuscript/conformal_manuscript.md
- Numbers: code/results/conformal_metrics.json ← code/src/conformal/analyze.py ← conformal_results.csv
- Fill/build/SI: manuscript/{fill_conformal.py, build_conformal.py, si_conformal.py}
- Study: code/run_conformal.py (+ run_strata_sens.py, run_profile.py); methods in code/src/conformal/
- Reproduction command chain: code/README.md §3a
- Verified Crossref citations incl. prior art: references.bib (24 entries; dAD 2023, Hosni 2026,
  Papadopoulos-kNN 2011)

## Process record
- research_notes/research_plan_v2.md — pre-registration of the pivot
- REVIEW_HISTORY_conformal.md — 8-reviewer adversarial panel (round 1) + dispositions; round-2 verification
- REJECTION_RISK_conformal.md — pre-submission risk assessment (post-panel status)
- FINAL_DECISION_conformal.md — journal selection + honest timeline disclosure
- CHANGELOG.md / EXPERIMENT_LOG.md / PROJECT_STATE.md — running record

## Author to complete before submission
Author names/affiliations/ORCIDs, acknowledgments, funding; insert repository URL on acceptance
(currently a policy placeholder under MIT license).

## Honest scope (stated in the paper)
Approximate (not exact) per-stratum coverage; smoothed-AD carries no guarantee; cal/test only
approximately exchangeable under scaffold splits; one AD descriptor; regression only; CPU-scale base
models; MoleculeACE dominates the dataset pool (per-collection breakdown + win-rate reported). Deferred
to future work: temporal splits, GNN/pretrained base models, alternative AD descriptors, CQR/localized-CP.
