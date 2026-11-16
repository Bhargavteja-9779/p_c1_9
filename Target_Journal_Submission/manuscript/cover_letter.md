---
title: Cover Letter
---

[Date]

Dear Editors of the *Journal of Chemical Information and Modeling*,

We are pleased to submit our manuscript, **"Chemically Equivalent, Computationally Different:
Standardization-Invariance Violations in Molecular Property Prediction and a Representation-Agnostic
Remedy,"** for consideration as an Article in the *Journal of Chemical Information and Modeling*.

Every molecular machine-learning pipeline begins by standardizing chemical structures — choosing a
tautomer, a protonation and salt state, a stereochemistry convention, and an aromaticity model.
These choices are individually defensible, seldom reported, and differ between toolkits and groups.
Our manuscript asks, and answers with evidence, whether they affect what models predict.

**What we find.** Across a large, public, multi-dataset benchmark (regression and classification),
we show that standardization choice has only a modest effect on *aggregate* benchmark metrics — which
is precisely why the problem has gone unnoticed — but a substantial effect on *individual*
predictions: the same molecule, prepared under different but chemically-equivalent standardizations,
receives materially different predictions, by an amount comparable to a fraction of experimental
assay noise. We characterize this per-prediction standardization-invariance violation, attribute it
to specific operations, and show it is strongly representation-dependent (largest for fingerprint
linear models, smallest for tree ensembles).

**What we contribute.** We introduce SORB (Standardization-Orbit Robust learning), a
representation-agnostic method that augments training over, and averages inference across, each
molecule's standardization orbit. We evaluate it with multi-seed, multi-dataset statistics and a
controlled semi-synthetic experiment with a known invariance ground truth, and report honestly where
it helps and where it does not.

**Why JCIM.** This work is a cheminformatics methodology and best-practice contribution squarely
within the journal's scope, building directly on recent JCIM work on rigorous method comparison and
data quality. It adopts the practically-significant comparison protocol recently published in this
journal, and we provide complete open-source code and public-data instructions to meet the journal's
data and software availability requirements.

**Declarations.** This manuscript is original, has not been published previously, and is not under
consideration elsewhere. All data are public and require no credentials; all code is released openly.
The authors declare no competing financial interest. All authors have approved the manuscript and
agree to its submission.

We believe the molecular-modeling community will find this a useful and actionable contribution to
the reliability and reproducibility of property-prediction models. We thank you for your
consideration and look forward to your response.

Sincerely,

[Corresponding Author Name], on behalf of all authors
[Affiliation] · [Email]
