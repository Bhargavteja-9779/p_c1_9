# Phase 41 — Final Rejection-Risk Analysis

Honest assessment of residual risks for submission to the Journal of Chemical Information and Modeling.
No claim of zero risk or guaranteed acceptance is made.

| # | Risk | Evidence | Severity | Mitigation | Status |
|---|---|---|---|---|---|
| 1 | **Contribution judged insufficiently impactful**: the realistic per-prediction effect is modest (~0.07× label SD cross-pipeline), and the method (SORB) is a largely-negative result. | §6.4, §7.4; revision_results.csv | **High** | Reframed as a measurement + best-practice contribution with a concrete, quantified recommendation (consistent canonicalization) and a cautionary shortcut result; this genre is published in JCIM (Ash 2025; curation reviews). | Partially mitigated; genuine residual risk, disclosed. |
| 2 | Novelty vs Kikuchi 2025 (confounder concept) and the augmentation/TTA lineage. | §1, §2 | Medium | Kikuchi credited for the concept; SORB explicitly framed as the established paradigm on the graph axis; the graph axis (affecting fingerprints/descriptors) is the defensible delta. | Mitigated. |
| 3 | Headline could read as cherry-picked (0.47× is linear+fingerprint, full orbit). | §6.4 | Medium | Lead with realistic pipeline SD; orbit SD explicitly labeled a worst-case upper bound; tree-model values given. | Mitigated. |
| 4 | Endpoint-invariance validity (logD/permeability are protonation-coupled). | §11(ii) | Medium | Caveat added; tautomer-dominant effect is unaffected; protomer-defined endpoints excluded. | Mitigated (caveat, not re-run). |
| 5 | Benchmarks pre-standardized → effect is a lower bound; raw-data effect unquantified. | §11(i) | Medium | Disclosed as a limitation; framed as lower bound. | Disclosed, not eliminated. |
| 6 | Statistical power: 3 seeds (pre-registered 5); method on 48/54 datasets. | §5, §7 | Low–Medium | Dataset is the unit of analysis (48 datasets); deviation disclosed and justified by CPU budget; Friedman/Nemenyi + Holm reported. | Mitigated/disclosed. |
| 7 | Model scope: no GNN or chemical language models. | §11(iv) | Low–Medium | Disclosed; CPU-only constraint stated; future work. | Disclosed. |
| 8 | Reproducibility / data-software compliance. | README, run_all.py, code release | Low | Full open code, pinned versions, one-command reproduction, machine-generated figures/tables; public data. | Mitigated. |
| 9 | Semi-synthetic shortcut is artificial. | §7.5 | Low | Labeled controlled/semi-synthetic; complements, not replaces, observational results. | Mitigated. |
| 10 | Citation accuracy. | references.bib (26 entries) | Low | All DOIs/arXiv IDs verified via Crossref/arXiv; no fabricated references. | Mitigated. |
| 11 | Journal scope/formatting (ACS). | DOCX, Data/Software statement | Low | ACS reference style via CSL; Data & Software Availability statement included; within length norms (~5.2k words). Final ACS template/word-limit check recommended before submission. | Mitigated; confirm guidelines at submission. |

## Summary
The dominant risk is editorial judgment of impact (Risk 1): the paper is a rigorous, honest, useful
measurement-and-best-practice contribution rather than a high-impact new method, and it reports a
negative method result. This is disclosed rather than hidden. All identified FATAL and MAJOR
methodological issues from the ten-reviewer round have been resolved or honestly scoped.
