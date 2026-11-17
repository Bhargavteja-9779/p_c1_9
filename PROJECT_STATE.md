# PROJECT_STATE

Last updated: 2026-11-17 (session continues across container restarts)

## PIVOT (current primary deliverable)
After the standardization-invariance paper was judged to have **limited novelty** (reviewer-style
rejection on contribution), the project pivoted to a **genuinely novel, strong-positive** finding on the
same benchmarks. Six candidate directions were empirically pre-checked; only one produced a strong,
defensible positive result, and it is now the submission.

### TOPIC (LOCKED, current)
**A False Sense of Security: Standard Conformal Prediction Under-Covers Novel Chemotypes, and an
Applicability-Domain-Conditional Fix.** Standard (marginal) split-conformal meets its marginal coverage
target yet systematically under-covers the most novel chemotypes (far applicability domain); the popular
normalized/heteroscedastic variant retains marginal coverage but does not repair the conditional gap; an
**applicability-domain-conditional (Mondrian-AD) conformal** procedure (+ a smoothed/shrinkage variant)
restores novel-compound coverage to nominal with a finite-sample per-stratum guarantee and no retraining.
Pre-registration: research_notes/research_plan_v2.md.

### JOURNAL (SELECTED): Journal of Chemical Information and Modeling (ACS), Q1.
Strong scope fit (conformal + applicability domain + QSAR; several cited refs are JCIM). Timeline handled
honestly (first-decision vs acceptance distinguished; ≤100 d targeted, not guaranteed). See
Target_Journal_Submission/FINAL_DECISION_conformal.md.

## STATUS (conformal paper)
- Full study RE-RUNNING in background (detached nohup) after a baseline-fairness fix to normalized
  conformal (standard additive-β regularizer |r|/(σ+β); β=median σ). Command to resume on restart:
  `cd Target_Journal_Submission/code && python3 run_conformal.py --mode full --nproc 3 --splits scaffold,random --alphas 0.1`
  (RESUMES from results/_confcells/ checkpoints). Writes results/conformal_results.csv.
- Grid: 41 regression datasets × [lgbm,rf,linear] × [ecfp,desc] × seeds[0,1,2] × {scaffold,random} × α=0.1 = 1476 cells.
- Pending after DONE: run α-robustness subset (α∈{0.05,0.2}); then `code/finalize_conformal.sh`.
- Validated (partial, 19 datasets, scaffold): marginal worst-bin ~0.80 → Mondrian-AD ~0.86 (+5pts,
  paired p<1e-4, ~95% datasets improved); normalized worst-bin low; scaffold-Mondrian no help. Friedman p<1e-11.

## CODE (conformal; tested end-to-end on partial data)
- code/src/conformal/{methods,experiment,analyze,figures}.py ; code/run_conformal.py (checkpointed)
- code/finalize_conformal.sh (analyze→figures→fill→build→SI→DOCX/PDF)
- Reuses: src/data/*, src/preprocessing/*, src/models/*, src/stat_tests/compare.py
- NOTE: normalized() now uses the textbook additive-β form (fairness fix, 2026-11-17).

## MANUSCRIPT (conformal)
- manuscript/conformal_manuscript.md: full draft, {{TOKENS}} (25) for every number; Abstract→Conclusion written.
- manuscript/fill_conformal.py (tokens ← conformal_metrics.json), build_conformal.py (crossrefs+figures),
  si_conformal.py (data-derived SI tables). references.bib (21 verified entries). ACS DOCX pipeline tested.
- manuscript/conformal_cover_letter.md ; Target_Journal_Submission/{FINAL_DECISION,REJECTION_RISK}_conformal.md

## NEXT EXACT STEPS (when study DONE)
1. (optional) α-robustness: `python3 run_conformal.py --mode full --models lgbm,rf,linear --feats ecfp,desc --seeds 0 --splits scaffold --alphas 0.05,0.2 --out conformal_alpha_results.csv`
2. `bash code/finalize_conformal.sh`  → verify 0 tokens missing; sanity-check numbers vs partial
3. VERIFY manuscript prose matches actual magnitudes (normalized narrative especially); claim-evidence audit
4. 10 adversarial review rounds (focus: novelty R1/R2, baseline fairness R3, CQR request) + revise
5. Final package manifest under Target_Journal_Submission/

## ARCHIVED (superseded, kept for provenance — NOT the submission)
Standardization-invariance paper + SORB: src/{preprocessing,models/sorb,training}, manuscript/manuscript.md,
Final_Manuscript.docx/pdf, cover_letter.*, REVIEW/, FINAL_DECISION.md, REJECTION_RISK_REPORT.md,
research_notes/research_plan.md. Rejected for limited novelty; retained as negative-result provenance.
