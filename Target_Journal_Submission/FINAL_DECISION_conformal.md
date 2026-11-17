# Final submission decision — conformal paper

## Manuscript
**Title:** A False Sense of Security: Standard Conformal Prediction Under-Covers Novel Chemotypes,
and an Applicability-Domain-Conditional Fix

**One-line contribution:** We quantify, at benchmark scale, that standard (marginal) conformal
prediction systematically under-covers novel chemotypes in molecular property regression, show the
popular normalized variant does not repair it, and introduce an applicability-domain-conditional
(Mondrian-AD) conformal procedure that restores novel-compound coverage to nominal with a finite-sample
per-stratum guarantee and no retraining.

## Target journal (primary): *Journal of Chemical Information and Modeling* (JCIM)
**Publisher:** American Chemical Society (ACS) — a prestigious society publisher on the preferred list
(not MDPI).

**Quartile / indexing:** Q1 in Journal Citation Reports under *Chemistry, Multidisciplinary* and
*Computer Science, Interdisciplinary Applications* (well-established; the journal is a flagship
cheminformatics venue). *Verify the exact current quartile and impact factor against the latest JCR
release at submission time.*

**Scope fit (strong):** JCIM's stated scope covers cheminformatics, QSAR/QSPR, machine learning for
molecular properties, applicability domain, and uncertainty quantification — precisely this paper's
subject. Independent evidence of fit: several works we build on and cite are JCIM papers — Norinder
et al. 2014 (conformal prediction + applicability domain), Sun/Svensson et al. 2017 (Mondrian
cross-conformal in cheminformatics), Sheridan 2013 (time-split validation). The paper speaks directly
to this readership and gives an actionable recommendation.

**Formatting:** The manuscript and SI are built to ACS style (ACS reference CSL + ACS Word template);
figures are 600 dpi; the paper is an Article-length contribution.

## Timeline (disclosed, not guaranteed)
The user's constraint is submission-to-**acceptance** within ~100 days. Important honest caveats:
- No journal can *guarantee* an acceptance date; editorial and review timing varies per manuscript.
- ACS journals including JCIM typically return a **first decision** within roughly 1–2 months; this is
  *first decision*, which is distinct from, and earlier than, final acceptance.
- I do **not** have access to JCIM's live, current editorial statistics, so any specific median
  acceptance time I might state would be unverified. I therefore do not assert a guaranteed ≤100-day
  acceptance; I judge JCIM a reasonable choice to *target* that window given typical ACS turnaround and
  the paper's self-contained, reproducible nature (which tends to shorten review), but the final
  timeline is outside the authors' control.

## Backup journals (if desk-rejected or out of scope at JCIM)
1. **Journal of Cheminformatics** (Springer Nature; Q1; open access). Strong fit (conformal prediction
   for cheminformatics; CPSign 2024 is a J. Cheminform. paper). Fully open, often fast.
2. **Digital Discovery** (Royal Society of Chemistry; Q1-tier, newer). ML-for-chemistry methods, good
   fit for a methods-and-benchmark contribution.
3. **Journal of Computer-Aided Molecular Design** (Springer; Q2–Q1). UQ/applicability-domain fit.

## Rationale for not choosing MDPI
Per the user's instruction, MDPI venues (e.g., *Molecules*, *IJMS*) are excluded by default in favor of
society/established commercial publishers (ACS, Springer Nature, RSC).

## Submission readiness checklist
- [x] Manuscript (ACS style), all numbers machine-generated from result files
- [x] Supporting Information (data-derived tables; robustness)
- [x] Cover letter (JCIM-addressed)
- [x] Open code + seeds + environment; scripts regenerate every figure/table
- [ ] Author names, affiliations, ORCIDs, acknowledgments, funding (author to complete)
- [ ] Repository URL inserted on acceptance (currently a placeholder, by policy)
- [ ] Final adversarial review pass + rejection-risk report (in progress)
