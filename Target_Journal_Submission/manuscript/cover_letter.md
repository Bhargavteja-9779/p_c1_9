---
title: Cover Letter
---

[Date]

Dear Editors of the *Journal of Chemical Information and Modeling*,

We are pleased to submit our manuscript, **"Chemically Equivalent, Computationally Different:
Standardization-Induced Prediction Variability in Molecular Property Models, and a Simple Fix,"** for
consideration as an Article.

Every molecular machine-learning pipeline begins by standardizing chemical structures — choosing a
tautomer, a protonation and salt state, a stereochemistry convention, and an aromaticity model. These
choices are individually defensible, seldom reported, and differ between toolkits and groups. Our
manuscript asks, and answers with evidence across 54 public benchmarks and three model families,
whether and how much they affect what models predict.

**What we find.** Standardization choice has only a small effect on *aggregate* benchmark metrics
(median relative change ~0.5%), which is why it is easily overlooked, but it changes *individual*
predictions more than the aggregate numbers suggest — for a single molecule, re-standardization shifts
a prediction by a median of ~0.12 of the label standard deviation over its full set of equivalent
forms (up to ~0.47 for high-variance fingerprint linear models), driven chiefly by tautomer choice. We
then show the remedy is simple: applying one consistent canonical-tautomer standardization to both
training and inference inputs collapses this per-prediction sensitivity to essentially zero. We
evaluate a more elaborate orbit-robust training method (SORB) honestly and report, as a negative
result, that it does not reduce cross-tool disagreement and carries a real accuracy cost, though it
does resist a standardization "shortcut" that we demonstrate with a controlled experiment.

**Why JCIM.** This is a cheminformatics methodology and best-practice contribution squarely within the
journal's scope, extending the curation literature with a specific, quantified recommendation and
adopting the practically-significant method-comparison protocol recently published in this journal. We
provide complete open-source code and public-data instructions, and all figures and tables are
regenerated from machine-readable result files, to meet the journal's data and software availability
requirements.

**Declarations.** This manuscript is original, has not been published previously, and is not under
consideration elsewhere. All data are public and require no credentials; all code is released openly.
The authors declare no competing financial interest. All authors have approved the submission. We note
that our central message is a cautionary and corrective one — quantifying an overlooked confounder and
showing that a simple, well-understood fix suffices — which we believe is of direct practical value to
the molecular-modeling community.

We thank you for your consideration.

Sincerely,

[Corresponding Author Name], on behalf of all authors
[Affiliation] · [Email]
