# Selected Target Journal

**TARGET JOURNAL:** Journal of Chemical Information and Modeling (JCIM)
**PUBLISHER:** American Chemical Society (ACS)
**ISSN:** 1549-9596 (print) · 1549-960X (online)
**HOMEPAGE:** https://pubs.acs.org/journal/jcisd8

| Field | Value | Verification status |
|---|---|---|
| JCR indexed (SCIE) | Yes | Secondary (MIAR/Clarivate MJL; ACS metrics pages). Clarivate JCR itself requires login. |
| Journal Impact Factor | ≈ 6.4 (JIF year 2025, June 2026 JCR release); 5.6 (2024) | Secondary aggregators (scijournal, slogix). To re-confirm on the ACS journal page at submission. |
| JCR quartile | Q1 (Computer Science, Information Systems; Chemistry, Multidisciplinary) | Secondary; Q1 consistently reported. Per-category rank not independently read from JCR. |
| JCR subject categories | Chemistry, Multidisciplinary; Computer Science, Information Systems; Computer Science, Interdisciplinary Applications | Secondary |
| APC (gold OA, optional) | ACS Open/AuthorChoice CC-BY ≈ USD 2500–5000 (tiered); subscription (no-fee) route available | To confirm exact figure at submission |
| Submission→acceptance (our census) | **median 84 days** (n = 347 research articles, published 2025-11-15…2026-11-15), Q1 63 d, Q3 112 d; **64% ≤ 100 days**, 22.5% ≤ 60 days | Primary: computed from PubMed `<History>` received→accepted dates (`research_notes/scripts/pubmed_timeline_census.py`) |
| Decision timeline | Satisfies the ≤100-day (and frequently ≤60-day) target | Primary census above |

## Why JCIM

1. **Scope fit (verified).** JCIM is the flagship ACS venue for cheminformatics and molecular-ML methodology. A Crossref scan of 2022–2026 JCIM articles (`research_notes/jcim_methodology_corpus.json`, 235 methodology-adjacent papers) shows the journal actively publishes exactly our genre: uncertainty quantification, applicability domain, probability calibration, evaluation/benchmark methodology, data-quality and curation studies, and molecular representation analysis. Directly relevant 2026 examples include "Uncertainty Quantification in Molecular Machine Learning for Property Prediction", "Toward Imbalanced Molecular Property Regression: A Benchmark Study", a "Three-Level Noise Floor" data-quality paper, "AssayMatch", and "Practically Significant Method Comparison Protocols for Machine Learning in Small Molecule Drug Discovery" (Ash et al. 2025, 10.1021/acs.jcim.5c01609), whose comparison protocol we adopt.
2. **Acceptance timeline (primary evidence).** Our own date census gives a median submission→acceptance of 84 days, comfortably inside the user's ≤100-day requirement and far better than most high-IF AI/Elsevier/Springer venues we measured (e.g., Knowledge-Based Systems / ESWA / Information Sciences medians >140 days; see `timeline_census_table.md`).
3. **Q1 JCR, reputable publisher.** ACS is a legitimate society publisher; JCIM is long-established with rigorous editorial standards and a strict Data & Software Availability policy (ACS editorial 10.1021/acs.jcim.0c01389), which our fully-reproducible, open-code study is built to satisfy.
4. **Methodological + reproducibility emphasis.** The selected topic (standardization-induced benchmark instability and a standardization-robust method) is a cheminformatics best-practice + method contribution, precisely the kind of work JCIM's editors and reviewers value.

## Honest caveats
- JIF value/quartile taken from aggregators because Clarivate JCR requires login; to be re-verified on the ACS page at submission (does not affect topic suitability).
- Acceptance-time census uses publisher-deposited received dates; if a journal resets the "received" date on major revision, true first-submission-to-acceptance could be longer. The median (84 d) and the ≤100-day fraction (64%) are robust to this caveat relative to the 100-day target.
