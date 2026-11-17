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
coverage; what has not been done, to our knowledge, is to *quantify this specific failure for molecular
models along the applicability-domain axis at benchmark scale*, to show that the community's go-to
remedy — normalized/heteroscedastic conformal — does not address the problem (it preserves marginal
coverage but only reallocates interval width by a learned difficulty estimate that is itself unreliable
out of domain, leaving novel-chemotype coverage low), and to provide a simple, guaranteed fix. We introduce **applicability-domain-conditional (Mondrian-AD)
conformal**, which calibrates separately across strata defined by distance to the training set, with a
shrinkage variant for sparse strata. It restores novel-compound coverage to nominal with a
finite-sample per-stratum guarantee, at a modest and well-placed width cost, with no retraining and no
dependence on the base model or featurizer.

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
