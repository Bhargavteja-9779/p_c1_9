# CHANGELOG

(created 2026-11-15; entries appended chronologically)
- 2026-11-15: project skeleton and state files created.
- 2026-11-16: resumed after usage limit; journal census scripts; JCIM shortlisted; TDC ADMET + MoleculeACE cached.
- 2026-11-16: Phase-4 novelty defense killed C1,D,delta,P2 (all pre-empted by 2023-2026 work). Pivoting to lower-pre-emption directions.
- 2026-11-16: TOPIC LOCKED = standardization-induced benchmark instability + SORB method, for JCIM. Research plan pre-registered.
- 2026-11-16: analysis+figures pipeline; pre-registration amendment (invariance-centric framing) after blinded pilot.
- 2026-11-16: container restart during precompute (24/55); relaunching checkpointed full run.
- 2026-11-17: full grid executing (sens ~820/990); state snapshot for restart resilience.
- 2026-11-17: full experiment complete; first manuscript assembled (Final_Manuscript.docx), all 26 numeric tokens from machine results.
- 2026-11-17: Revision complete (honest reframing: canonicalization is the fix, SORB negative result). Final audits, rejection-risk report, READY-FOR-SUBMISSION decision. Package complete.
- 2026-11-17: PIVOT to AD-conditional conformal prediction (validated strong positive) after standardization paper judged insufficient novelty; 6 pilot validations recorded.
- 2026-11-17: Conformal study built (src/conformal/*, run_conformal.py, checkpointed); quick-mode validated. Full study launched (1476 cells).
- 2026-11-17: BASELINE FAIRNESS FIX — normalized conformal changed from hard 1e-6 floor to textbook additive-beta score |r|/(sigma+beta), beta=median(sigma); half-width q*(sigma+beta). Wiped _confcells; restarted full study. Rationale: old form caused width blow-up (mean rel-width ~1e5) and under-coverage artifact; new form retains marginal coverage (as theory requires) so our claim narrows honestly to "normalized does not repair the CONDITIONAL gap". Manuscript/cover-letter normalized narrative revised accordingly.
- 2026-11-17: Conformal manuscript drafted (Abstract->Conclusion, 25 tokens), fill_conformal.py / build_conformal.py / si_conformal.py / finalize_conformal.sh; cover letter, FINAL_DECISION_conformal, REJECTION_RISK_conformal written. Downstream pipeline smoke-tested end-to-end on partial data -> clean DOCX, 0 missing tokens.
- 2026-11-17: BASELINE CORRECTNESS FIX #2 — normalized conformal's difficulty model sigma(x) was being fit on calibration residuals and evaluated in-sample on calibration, breaking score exchangeability (full-run normalized marginal coverage 0.769, invalid). Fixed to fit sigma on OUT-OF-FOLD training residuals (2-fold cross-fit), disjoint from cal/test -> valid near-nominal marginal coverage AND genuine heteroscedastic sigma. Smoke-test confirmed (normalized marginal ~0.88, modest far-coverage help, mondrian_ad still best). Wiped _confcells, relaunched full study. Manuscript/SI/risk-report normalized narrative updated to the honest "modest help at best, not a fix" finding.
