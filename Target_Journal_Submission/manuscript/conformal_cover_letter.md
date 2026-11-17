**[Author Name(s)]**
[Affiliation]
[Address]
[corresponding@author.email]

[Date]

**To the Editors,**
*Journal of Chemical Information and Modeling*

Dear Editors,

We are pleased to submit our manuscript, **"A False Sense of Security: Standard Conformal
Prediction Under-Covers Novel Chemotypes, and an Applicability-Domain-Conditional Fix,"** for
consideration as an Article in the *Journal of Chemical Information and Modeling*.

**The problem.** Conformal prediction is being adopted across cheminformatics because it promises a
finite-sample, distribution-free coverage guarantee for molecular property models — a principled
alternative to heuristic applicability-domain scores. That guarantee, however, is *marginal*: it holds
on average over the test set. We show, at scale, that this averaging conceals a failure that strikes
exactly where uncertainty matters most. Across a large, public benchmark suite of molecular-property
regression datasets, three model families, and two featurizations, standard split-conformal intervals
meet their nominal marginal coverage yet **systematically under-cover the most novel compounds** — the
far-applicability-domain chemotypes a prospective campaign most wants to score. The nominal 90%
interval a chemist would trust for a genuinely new candidate covers the truth markedly less often than
advertised: a quantified "false sense of security."

**Why it matters, and what is new.** It is textbook that marginal coverage does not imply conditional
coverage, and recent cheminformatics work has begun to address applicability-domain-aware conformal
calibration (notably the dynamic-calibration method of Oršolić and Šmuc, *Bioinformatics* 2023, and,
concurrently, Hosni et al., *ACS Omega* 2026). Our contribution is the *systematic, benchmark-scale*
characterization these do not provide — 41 single-task molecular-regression datasets, three model
families, two featurizations, with coverage resolved continuously along the applicability domain — and,
with the obvious remedies ruled out as negative controls (a learned-difficulty normalized score, a
distance-normalized score, and a scaffold-conditional variant), a demonstration that a *simple,
interpretable, retraining-free* applicability-domain-conditional wrapper (Mondrian-AD, with a shrinkage
variant) raises novel-compound coverage to nominal. We are explicit about scope: the restored per-stratum
coverage is approximate (the strata are data-dependent and the most novel stratum is open-ended), the
gain is in *honest* uncertainty rather than sharper predictions, and we verify the fix improves coverage
even on an independent (scaffold) partition it was not calibrated on. The method is model- and
featurization-agnostic and needs no retraining.

**Rigor and reproducibility.** Every quantitative claim in the manuscript is generated directly from
machine-readable result files by a token-substitution pipeline — no number is hand-transcribed. We use
scaffold (out-of-distribution) splits as the primary protocol, report paired nonparametric tests with
bootstrap confidence intervals and a Friedman/Nemenyi omnibus across methods, and validate the finding
and the fix across both featurizations, all three model families, the random split, and multiple
significance levels. We also report, honestly, the method's cost (wider intervals) and its limits
(per-stratum rather than object-conditional validity). All code, environment specifications, seeds, and
scripts that regenerate every figure and table from the raw results are released under an open license.

**Fit for JCIM.** The work sits squarely within the journal's scope on predictive modeling,
applicability domain, and uncertainty quantification for chemical data, and it offers the readership an
immediately actionable recommendation: report and calibrate conformal coverage *conditional on the
applicability domain*, so that the uncertainty attached to a novel candidate can actually be trusted.

This manuscript is original, is not under consideration elsewhere, and all authors have approved the
submission. We have no conflicts of interest to declare. We thank the editors and reviewers in advance
for their time.

Sincerely,

[Author Name(s)], on behalf of all authors
