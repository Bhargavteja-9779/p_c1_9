---
title: "A False Sense of Security: Standard Conformal Prediction Under-Covers Novel Chemotypes, and an Applicability-Domain-Conditional Fix"
author:
  - name: "[Author Name(s)]"
    affiliation: "[Affiliation]"
    email: "[corresponding@author.email]"
keywords: "conformal prediction; uncertainty quantification; applicability domain; QSAR; molecular property prediction; distribution shift; prediction intervals; drug discovery"
bibliography: references.bib
csl: assets/acs.csl
---

# Abstract

Reliable uncertainty is most valuable exactly where molecular property models are least reliable: on
novel chemotypes far from the training set, the compounds a prospective drug-discovery campaign most
wants to score. Conformal prediction is increasingly adopted for this purpose because it provides a
finite-sample, distribution-free coverage guarantee. That guarantee, however, is **marginal** —
averaged over the test distribution — and we show this averaging hides a systematic and
decision-relevant failure. Across {{N_DATASETS}} public molecular-property regression datasets, three
model families, and two featurizations, standard split-conformal intervals attain their nominal
marginal coverage ({{MARG_COV}} at a 90% target) yet **under-cover the most novel compounds**: for
test molecules in the farthest applicability-domain stratum, coverage falls to a median of
{{MARG_FAR_COV}}, and the worst-stratum coverage is {{MARG_WORST_COV}} — a false sense of security for
the extrapolation regime. The widely-recommended normalized (heteroscedastic) conformal variant does
not repair this: it retains marginal coverage but reallocates interval width by a learned difficulty
estimate that is itself unreliable out of domain, leaving novel-compound coverage low (worst-stratum
{{NORM_WORST_COV}}). We introduce an **applicability-domain-conditional
conformal** procedure that calibrates separately across strata defined by distance to the training set
(a Mondrian taxonomy on the applicability domain), with a shrinkage variant for sparse strata. It
restores novel-compound coverage to nominal (worst-stratum coverage {{MOND_WORST_COV}}, a
{{WORST_GAIN}}-point gain; paired $p={{WORST_P}}$) at a modest width cost ({{WIDTH_COST}}%), while
preserving marginal coverage. The method is model- and featurization-agnostic, requires no retraining,
and is validated across splits and significance levels. We provide an open, reproducible implementation
and a simple recommendation: report and calibrate conformal coverage conditional on the applicability
domain, so that the uncertainty attached to a novel candidate can actually be trusted.

# 1. Introduction

Machine-learning models of molecular properties — potency, ADMET endpoints, physicochemical
quantities — are now routine in drug discovery and materials design. Their point predictions are
useful, but decisions (synthesize, assay, prioritize) require knowing *how much to trust* a given
prediction, and that trust is most consequential for **novel** molecules that are structurally
dissimilar from anything in the training set. These are precisely the compounds a model is least
equipped to predict, and precisely the ones a prospective campaign is most interested in.

Conformal prediction [@vovk2005; @papadopoulos2002; @angelopoulos2023] has become a popular way to
attach calibrated uncertainty to molecular models [@norinder2014; @svensson2018; @cpsign2024] because,
unlike heuristic applicability-domain scores, it offers a finite-sample, distribution-free guarantee:
a conformal prediction interval at level $1-\alpha$ contains the true value with probability at least
$1-\alpha$, under only the assumption that calibration and test data are exchangeable. This guarantee
is attractive and, importantly, it is the *marginal* coverage — the probability averaged over the whole
test distribution. It is well known in the conformal literature that marginal coverage does **not**
imply conditional coverage [@vovk2012conditional]: a predictor can hit its target on average while
systematically under-covering some regions and over-covering others.

For molecular property prediction this theoretical caveat has a sharp, practical consequence that, to
our knowledge, has not been quantified at scale. If standard conformal intervals under-cover exactly
the novel-chemotype region — where uncertainty matters most — then the nominal 90% interval a chemist
reads off for a new candidate may in truth cover far less than 90% of the time, a false sense of
security in the highest-stakes case. Conversely, over-covering familiar compounds wastes nothing but
attention. The question is empirical: how large is this effect, and can it be fixed cheaply and
reliably?

**Contributions.**

1. **A quantified "false security" finding.** Across {{N_DATASETS}} public molecular-property
   regression datasets, three model families, and two featurizations, we show that standard
   split-conformal intervals meet their marginal target yet systematically under-cover novel
   chemotypes (coverage degrading monotonically with distance to the training set), and that the
   popular normalized/heteroscedastic variant retains marginal coverage but does not repair the
   novel-chemotype gap (Section 6).
2. **An applicability-domain-conditional conformal method.** We calibrate conformal quantiles
   separately across strata of a molecule's distance to the training set — a Mondrian taxonomy on the
   applicability domain — with a shrinkage variant for sparsely-populated strata. The procedure keeps
   the finite-sample guarantee per stratum and restores novel-compound coverage to nominal at a modest
   width cost, without retraining and independent of model and featurization (Section 4, 7).
3. **A comprehensive, reproducible benchmark and a practical recommendation.** We compare marginal,
   normalized, scaffold-conditional, and applicability-domain-conditional conformal across datasets,
   models, splits, and significance levels, with multi-dataset statistics, and recommend reporting and
   calibrating conformal coverage conditional on the applicability domain (Section 7, 8, 10).

# 2. Related Work

**Conformal prediction and its guarantees.** Inductive (split) conformal prediction [@papadopoulos2002;
@vovk2005] attaches distribution-free prediction sets to any point predictor with a finite-sample
*marginal* coverage guarantee; see [@angelopoulos2023] for a modern introduction. That the guarantee is
marginal — not conditional — is classical [@vovk2012conditional]: exact object-conditional coverage is
impossible in finite samples without distributional assumptions, motivating approximate schemes.

**Approaches to conditional coverage.** Two broad strategies improve conditional behavior:
normalized/heteroscedastic nonconformity scores that scale the residual by an estimated difficulty
$\sigma(x)$ [@papadopoulos2002], including conformalized quantile regression [@romano2019]; and Mondrian
(group-conditional) conformal prediction, which calibrates separately within a partition of object
space and thereby guarantees coverage *within each group* [@vovk2005]. Mondrian conformal has been used
in cheminformatics chiefly to handle class imbalance in bioactivity classification [@svensson2018].
Localized conformal prediction [@guan2023] reweights calibration points near the test object. We find,
empirically, that the normalized variant does not fix the novel-chemotype coverage problem for
molecular regression — it preserves marginal coverage but reallocates width by a difficulty estimate
that is unreliable out of domain — whereas a Mondrian taxonomy defined directly on the *applicability
domain* does.

**Conformal prediction and applicability domain in cheminformatics.** Norinder et al. [@norinder2014]
proposed conformal prediction as a statistically-grounded alternative to heuristic applicability-domain
scoring for QSAR, and tools such as CPSign [@cpsign2024] make conformal modeling routine. The
applicability domain itself — how far a query is from the training set — is a long-standing concept for
flagging unreliable predictions, often via similarity or distance measures, and is tied to the harder
generalization seen under scaffold and time splits [@sheridan2013]. Our contribution connects these
threads: we use the applicability domain not as a separate heuristic flag but as the *conditioning
variable* for conformal calibration, and we quantify, at scale, the conditional-coverage failure that
makes this necessary.

# 3. Problem Formulation

Let $f$ be a trained regression model and let a calibration set $\{(x_i, y_i)\}_{i=1}^{n}$, exchangeable
with the test data, give nonconformity scores (absolute residuals) $s_i = |y_i - f(x_i)|$. Split
conformal forms the interval $f(x) \pm \hat q$ with $\hat q$ the $\lceil (n+1)(1-\alpha)\rceil/n$
empirical quantile of $\{s_i\}$; this guarantees $\Pr(y \in \text{interval}) \ge 1-\alpha$ marginally.
We are interested in **applicability-domain-conditional coverage**: writing $d(x)$ for a molecule's
distance to the training set (one minus its mean Tanimoto similarity to the $k$ nearest training
molecules), we ask whether $\Pr(y \in \text{interval}\mid d(x) \in B) \approx 1-\alpha$ for strata $B$
ranging from familiar ($d$ small) to novel ($d$ large). A clinically or synthetically useful guarantee
is one that holds in the novel stratum, not merely on average.

# 4. Method: Applicability-Domain-Conditional Conformal

**Mondrian-AD.** We partition molecules into $m$ strata by quantiles of the calibration-set distance
$d(\cdot)$, and compute a *separate* conformal quantile $\hat q_B$ within each stratum $B$. A test
molecule is assigned the quantile of the stratum its distance falls into, giving the interval
$f(x) \pm \hat q_{B(x)}$. Because within each stratum the calibration and test molecules are
exchangeable, Mondrian conformal retains the finite-sample coverage guarantee *conditional on the
stratum* [@vovk2005]. Intuitively, novel test molecules are calibrated against other novel calibration
molecules — which share their larger residuals — rather than against the easy bulk, so their intervals
widen to the level they actually require.

**Smoothed-AD (shrinkage).** The weakness of a hard taxonomy is that the farthest stratum may be
sparsely populated, making its quantile noisy. We therefore shrink each stratum quantile toward the
pooled quantile, $\tilde q_B = (1-\lambda)\hat q_B + \lambda\,\hat q_{\text{pooled}}$, trading a little
conditional accuracy for stability in the sparse-novel regime. We report both variants.

Both methods are post-hoc wrappers: they require only calibration residuals and the cheap
distance-to-training descriptor, no retraining, and apply to any model and featurization. We compare
against three baselines: **marginal** split conformal (the standard); **normalized** conformal with the
textbook additive-regularized score $|r|/(\sigma(x)+\beta)$ and interval half-width
$\hat q\,(\sigma(x)+\beta)$ [@papadopoulos2002; @romano2019], where the difficulty model $\sigma(x)$ is a
gradient-boosted regressor of the absolute residual, $\beta$ is the median calibration difficulty, and —
crucially for validity — $\sigma$ is fit on *out-of-fold training* residuals (a two-fold cross-fit),
disjoint from the calibration and test sets, so that normalized conformal keeps its marginal guarantee
rather than degenerating; and **Mondrian-scaffold** (the same Mondrian idea conditioned on Bemis–Murcko
scaffold cluster rather than on the applicability domain). We verified that the normalized baseline
attains near-nominal marginal coverage, confirming it is correctly specified and not handicapped.

# 5. Experimental Setup

**Datasets.** {{N_DATASETS}} public molecular-property regression benchmarks: the MoleculeACE
collection (30 ChEMBL potency sets) [@vantilborg2022], the Therapeutics Data Commons ADMET regression
tasks [@tdc2021], and the MoleculeNet physicochemical sets ESOL and FreeSolv [@moleculenet2018]
(MoleculeNet Lipophilicity duplicates the TDC AstraZeneca logD benchmark and is removed). All public,
no credentials.

**Models and featurization.** Three model families spanning the linear-to-nonlinear range — ridge
regression, random forest, and gradient-boosted trees (LightGBM) [@lightgbm2017] — on two
featurizations: 2048-bit ECFP4 fingerprints [@ecfp2010] and 2D RDKit descriptors. The applicability-
domain descriptor is one minus the mean top-5 Tanimoto similarity of a molecule to the training set.

**Protocol.** Each dataset is split three ways (train / calibration / test). The primary split is a
Bemis–Murcko scaffold split [@bemismurcko1996] (so that test molecules are genuinely out-of-scaffold,
making the applicability-domain gradient meaningful); a random split is reported as a robustness check.
Conditional coverage is evaluated by binning the test set into {{N_BINS}} equal-frequency
applicability-domain strata (defined on the calibration distribution) and measuring coverage, worst-
stratum coverage, the maximum stratum-wise deviation from nominal (the conditional-coverage gap), and
interval width per stratum. The primary target is $1-\alpha=0.9$ ($\alpha=0.2$ and $0.05$ as
robustness). Each configuration uses three seeds; the dataset is the unit of analysis.

**Statistics.** Per-method values are averaged over seeds to one value per dataset and compared across
datasets by paired Wilcoxon signed-rank tests with bootstrap 95% confidence intervals, and by a
Friedman omnibus with Nemenyi critical-difference analysis over the methods [@demsar2006]. All code is
released; every figure and table is regenerated from machine-readable result files.

# 6. Results: Standard Conformal Gives a False Sense of Security

## 6.1 Marginal coverage is met — on average
At the 90% target, split-conformal intervals attain near-nominal marginal coverage across the
{{N_DATASETS}} datasets (median {{MARG_COV}}; Table @tbl-summary), confirming the familiar guarantee.
Naive model error bars, by contrast, are badly miscalibrated under scaffold splits (median coverage
far below nominal), underscoring why a conformal procedure is needed in the first place.

## 6.2 But novel chemotypes are systematically under-covered
Coverage is not uniform across the applicability domain. Binning the test set into {{N_BINS}} strata by
distance to the training set, marginal conformal over-covers familiar compounds (nearest stratum
coverage {{MARG_NEAR_COV}}) and under-covers novel ones: coverage in the farthest stratum is
{{MARG_FAR_COV}} and the worst-stratum coverage is {{MARG_WORST_COV}}, a {{MARG_GAP}}-point maximum
deviation from nominal (Figure @fig-nearfar; Table @tbl-summary). The degradation is monotone in
distance and consistent across datasets, models, and featurizations. In decision terms, the nominal
90% interval a practitioner would trust for a genuinely novel candidate covers the truth only about
{{MARG_FAR_PCT}}% of the time.

## 6.3 Normalized conformal does not fix the conditional gap
The remedy most often recommended for conditional coverage — a normalized (heteroscedastic)
nonconformity score, $|r|/(\sigma(x)+\beta)$ with a learned difficulty model $\sigma(x)$ — does not
solve the problem. Implemented correctly — with the difficulty model fit on out-of-fold training
residuals so that the procedure remains valid — it retains marginal coverage (median {{NORM_COV}}) and
reallocates interval width by predicted difficulty, which yields at most a modest improvement in
novel-compound coverage: its worst-stratum coverage ({{NORM_WORST_COV}}) stays well short of nominal,
far from the restoration achieved by conditioning on the applicability domain. The difficulty axis is
not the applicability-domain axis, and the difficulty estimate is itself unreliable for the compounds
furthest from training. Conditioning instead on Bemis–Murcko scaffold cluster
(Mondrian-scaffold) likewise leaves the applicability-domain gradient essentially unchanged
(worst-stratum coverage {{SCAF_WORST_COV}}), because scaffold identity is not the axis along which
coverage fails.

# 7. Results: Applicability-Domain-Conditional Conformal Restores Coverage

## 7.1 The fix
Calibrating conformal quantiles conditional on the applicability domain (Mondrian-AD) restores
novel-compound coverage to nominal: the worst-stratum coverage rises from {{MARG_WORST_COV}} (marginal)
to {{MOND_WORST_COV}} (Mondrian-AD), a median gain of {{WORST_GAIN}} points (paired Wilcoxon
$p={{WORST_P}}$; 95% CI [{{WORST_CI_LO}}, {{WORST_CI_HI}}]; improved on {{WORST_FRAC}}% of datasets),
and the conditional-coverage gap falls from {{MARG_GAP}} to {{MOND_GAP}} (Figure @fig-nearfar,
@fig-perdataset). Marginal coverage is preserved ({{MOND_COV}}). The shrinkage variant (smoothed-AD)
achieves nearly the same conditional coverage (worst-stratum {{SMOOTH_WORST_COV}}) with slightly
tighter intervals, and is the more robust default when the novel stratum is sparse.

## 7.2 The cost, and the trade-off
Reliability for novel compounds is not free: Mondrian-AD widens intervals by {{WIDTH_COST}}% on average
(relative width {{MOND_RELWIDTH}} vs {{MARG_RELWIDTH}} for marginal), concentrated in the novel strata
where the extra width is warranted (Figure @fig-tradeoff). This is the correct behavior — the method
reports larger, honest uncertainty exactly where the model is extrapolating — and it is the price of a
coverage guarantee that actually holds out of domain. A Friedman omnibus over the four conformal
strategies on worst-stratum coverage (closeness to nominal) rejects equality ($p={{FRIEDMAN_P}}$; mean
ranks place Mondrian-AD and smoothed-AD best, marginal and normalized worst; Nemenyi critical
difference {{NEMENYI_CD}}).

## 7.3 Robustness
The finding and the fix are stable across conditions (Supporting Information): both featurizations, all
three model families, the random split (where the applicability-domain gradient is milder but marginal
conformal still under-covers the far stratum, and Mondrian-AD still restores it), and the $\alpha=0.2$
and $\alpha=0.05$ targets. We fix the stratum count at {{N_BINS}} and the shrinkage at $\lambda=0.3$
throughout; the hard-binned (Mondrian-AD, $\lambda=0$) and shrunk (smoothed-AD, $\lambda=0.3$) variants
already bracket the shrinkage choice, and a sensitivity analysis over the number of strata (Supporting
Information) shows the qualitative conclusion is unchanged, with more strata giving finer conditional
coverage at the cost of noisier sparse-stratum quantiles.

# 8. Downstream: Reliable Selection of Novel Candidates
The practical payoff is in prospective selection. When a campaign uses prediction intervals to screen
novel candidates — e.g., to flag molecules whose interval clears a potency threshold, or to abstain
where uncertainty is high — the guarantee must hold for the novel compounds under consideration.
Because marginal conformal under-covers that regime, its abstention and selection decisions are
miscalibrated exactly where it matters; applicability-domain-conditional conformal restores the
coverage on which such decisions rely, at the cost of flagging more novel compounds as uncertain — the
honest and desirable behavior (Supporting Information).

# 9. Discussion
Our results sharpen a known theoretical caveat into a concrete, measurable, and fixable problem for
molecular modeling. The marginal coverage guarantee that makes conformal prediction attractive is, by
itself, insufficient for the use case the field most cares about: trustworthy uncertainty on novel
chemotypes. Standard conformal intervals quietly under-cover exactly those compounds, giving a false
sense of security in prospective settings; the popular normalized variant does not help, because it
conditions on a learned difficulty estimate that is itself unreliable out of domain and is not aligned
with the applicability-domain axis along which coverage actually fails. The remedy is simple and
principled — condition the calibration on the applicability domain — and it inherits a finite-sample
per-stratum guarantee while remaining a cheap, model-agnostic, post-hoc wrapper. We recommend that
conformal coverage in molecular property prediction be reported and calibrated conditional on the
applicability domain, not merely marginally, so that the uncertainty attached to a novel candidate can
be trusted at its stated level.

# 10. Limitations
(i) Conditional validity is restored *per applicability-domain stratum*, not for every individual
molecule; exact object-conditional coverage is impossible in finite samples [@vovk2012conditional], and
extremely novel molecules may fall into a sparse stratum where the shrinkage variant is advisable.
(ii) The applicability-domain descriptor (Tanimoto distance to the training set) is one reasonable
choice; other descriptors (embedding distance, density) could define the strata and are left to future
work. (iii) We study regression; the same conditioning idea applies to classification but is not
evaluated here. (iv) Our models are CPU-feasible (fingerprint/descriptor learners); large pretrained
graph or language models, and conformalized quantile regression with strong base conditional models,
are natural extensions. (v) The benchmarks are curated and exchangeability may be only approximate
under strong temporal shift; the guarantee weakens gracefully as exchangeability is violated.

# 11. Conclusion
Standard conformal prediction gives molecular property models a coverage guarantee that holds on
average but fails where it matters most — on novel chemotypes — producing a false sense of security in
exactly the prospective, out-of-domain setting that drug discovery depends on. Across {{N_DATASETS}}
datasets we quantified this systematic under-coverage, showed that the common normalized variant retains
marginal coverage but does not repair it, and introduced an applicability-domain-conditional conformal procedure that
restores novel-compound coverage to nominal at a modest, well-placed width cost, with a finite-sample
per-stratum guarantee and no retraining. We release an open implementation and recommend conditioning
conformal calibration on the applicability domain as standard practice for trustworthy molecular
uncertainty.

# Data and Software Availability
All datasets are public and require no credentials (MoleculeACE [@vantilborg2022], Therapeutics Data
Commons [@tdc2021], MoleculeNet [@moleculenet2018]); download instructions are in the repository README.
Complete source code, environment specifications, deterministic seeds, and scripts that regenerate every
figure and table from machine-readable result files are released under the MIT License at [repository
URL — to be inserted on acceptance]. Models use scikit-learn [@sklearn2011] and LightGBM [@lightgbm2017];
fingerprints and descriptors use RDKit [@rdkit].

# Acknowledgments
[To be completed.]

# Author Contributions
[To be completed.]

# Competing Interests
The authors declare no competing financial interest.
