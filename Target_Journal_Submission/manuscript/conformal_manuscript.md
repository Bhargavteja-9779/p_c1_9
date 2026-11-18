---
title: "A False Sense of Security: Standard Conformal Prediction Under-Covers Novel Chemotypes, and the Case for Applicability-Domain-Conditional Calibration"
author:
  - "P. N. Bhargav Teja$^{1}$"
  - "Ragavan K$^{1,*}$"
keywords: "conformal prediction; uncertainty quantification; applicability domain; QSAR; molecular property prediction; prediction intervals; drug discovery"
bibliography: references.bib
csl: assets/acs.csl
---

$^{1}$[Affiliation — to be completed]

$^{*}$Corresponding author: Ragavan K ([corresponding e-mail — to be completed])

# Abstract

Reliable uncertainty is most valuable exactly where molecular property models are least reliable: on
novel chemotypes far from the training set, the compounds a prospective drug-discovery campaign most
wants to score. Conformal prediction is increasingly adopted for this purpose because it provides a
finite-sample, distribution-free coverage guarantee. That guarantee, however, is **marginal** —
averaged over the test distribution — and we show this averaging hides a systematic and
decision-relevant failure. Across {{N_DATASETS}} public molecular-property regression datasets, three
model families, and two featurizations, standard split-conformal intervals attain their nominal
marginal coverage ({{MARG_COV}} at a 90% target) yet **under-cover the most novel compounds**: for
test molecules in the farthest applicability-domain stratum, coverage falls to {{MARG_FAR_COV}} — a
false sense of security for the extrapolation regime. We show that *which conditioning axis* a remedy
uses is what matters: a learned-difficulty normalized score (farthest-stratum coverage {{NORM_FAR_COV}})
and a scaffold-conditional variant ({{SCAF_FAR_COV}}) do not repair the gap — the learned difficulty,
fit in-domain, grows only {{SIGMA_GROWTH}}% from familiar to novel compounds while the true error grows
{{RESID_GROWTH}}%, so it cannot react to novelty — whereas conditioning on the applicability-domain axis
does — both the established distance-normalized conformal ({{KNN_FAR_COV}})
and a simple applicability-domain Mondrian wrapper we describe ({{MOND_FAR_COV}}) raise novel-compound
coverage to near nominal (the applicability-domain-conditional methods are statistically indistinguishable
from one another and clearly separated from the rest by a Friedman/Nemenyi analysis). The Mondrian
wrapper and its shrinkage variant add an interpretable, approximately guaranteed per-stratum calibration
at comparable cost; across these methods the extra interval width falls on exactly the novel compounds,
turning width into an honest uncertainty signal usable for abstention. All are model- and
featurization-agnostic, need no retraining, and hold across splits and significance levels. We provide an
open, reproducible benchmark and a simple recommendation: report and calibrate conformal coverage
conditional on the applicability domain, so that the uncertainty attached to a novel candidate can
actually be trusted.

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

For molecular property prediction this theoretical caveat has a sharp, practical consequence. If
standard conformal intervals under-cover exactly the novel-chemotype region — where uncertainty matters
most — then the nominal 90% interval a chemist reads off for a new candidate may in truth cover less
than 90% of the time, a false sense of security in the highest-stakes case. That conformal uncertainty
degrades out of domain has been observed before in cheminformatics (most directly by the dynamic
applicability-domain work of Oršolić and Šmuc [@orsolic2023dad] for compound–target affinity, and
concurrently by Hosni et al. [@hosni2026ad], who use the applicability domain as a reliability *flag*).
What has not been done, and what we provide, is a *systematic, benchmark-scale* characterization for
single-task molecular regression — across many datasets, models, and featurizations, with coverage
resolved continuously along the applicability domain — together with a direct, honest head-to-head of
the simplest post-hoc remedies that isolates *which conditioning axis* actually matters. The questions
are empirical: how large is the effect, and which of the standard conformal variants fix it?

**Contributions.**

1. **A quantified "false security" benchmark.** Across {{N_DATASETS}} public molecular-property
   regression datasets, three model families, and two featurizations, we show that standard
   split-conformal intervals meet their marginal target yet systematically under-cover novel
   chemotypes — coverage in the most novel applicability-domain stratum falls well below nominal
   (Section 6).
2. **Isolating the conditioning axis that matters.** In a controlled head-to-head we show the failure is
   repaired by conditioning on the *applicability-domain* axis — and not by the alternatives: a
   learned-difficulty normalized score and a scaffold-conditional variant do *not* help, whereas both
   the established distance-normalized conformal [@papadopoulos2011knn] and a simple
   applicability-domain Mondrian wrapper (plus a shrinkage variant) that we describe restore
   novel-compound coverage to near nominal. The three applicability-domain-conditional methods are
   statistically indistinguishable from one another and clearly separated from the rest (Section 6, 7).
   We do not claim the Mondrian wrapper beats the existing distance-normalized baseline; it is a simple,
   interpretable alternative that additionally carries an approximate per-stratum guarantee, and the
   practical message is the *axis*, obtainable by either route. This complements the per-sample dynamic
   calibration of dAD [@orsolic2023dad], which targets the same axis by yet another mechanism.
3. **A comprehensive, reproducible benchmark and a practical recommendation.** We compare marginal,
   learned-difficulty normalized, distance-normalized, scaffold-conditional, and applicability-domain-
   conditional conformal across datasets, models, splits, and significance levels, with multi-dataset
   statistics and honest effect sizes; we show the restored interval width becomes an actionable novelty
   signal (Section 8), and recommend reporting and calibrating conformal coverage conditional on the
   applicability domain as standard practice (Section 7, 8, 10).

# 2. Related Work

**Conformal prediction and its guarantees.** Inductive (split) conformal prediction [@papadopoulos2002;
@vovk2005] attaches distribution-free prediction sets to any point predictor with a finite-sample
*marginal* coverage guarantee; see [@angelopoulos2023] for a modern introduction. That the guarantee is
marginal — not conditional — is classical [@vovk2012conditional]: exact object-conditional coverage is
impossible in finite samples without distributional assumptions, motivating approximate schemes.

**Approaches to conditional coverage.** Two broad strategies improve conditional behavior:
normalized/heteroscedastic nonconformity scores that scale the residual by an estimated difficulty
$\sigma(x)$ [@papadopoulos2002], including conformalized quantile regression [@romano2019] and
*distance-normalized* conformal, where $\sigma(x)$ is a $k$-nearest-neighbour distance to the training
set [@papadopoulos2011knn]; and Mondrian (group-conditional) conformal prediction, which calibrates
separately within a partition of object space and thereby targets coverage *within each group*
[@vovk2005]. Mondrian conformal has been used in cheminformatics chiefly to handle class imbalance in
bioactivity classification [@svensson2018]; localized conformal prediction [@guan2023] reweights
calibration points near the test object. Distance-normalized conformal and our Mondrian-AD condition on
the *same* applicability-domain axis — one by normalization, the other by stratified calibration — and
we benchmark both. We find that the learned-difficulty normalized variant does not repair the
novel-chemotype coverage gap, whereas both distance-normalized conformal and a Mondrian taxonomy on the
applicability domain do, to a statistically indistinguishable degree — the conditioning axis, not the
mechanism, is what matters.

**Conformal prediction and applicability domain in cheminformatics.** Norinder et al. [@norinder2014]
proposed conformal prediction as a statistically-grounded alternative to heuristic applicability-domain
scoring for QSAR, and tools such as CPSign [@cpsign2024] make conformal modeling routine. The
applicability domain itself — how far a query is from the training set — is a long-standing concept for
flagging unreliable predictions, often via similarity or distance measures, and is tied to the harder
generalization seen under scaffold and time splits [@sheridan2013]. Closest to our work, the dynamic
applicability-domain method of Oršolić and Šmuc [@orsolic2023dad] constructs *local*,
neighbourhood-specific calibration sets for compound–target binding-affinity prediction and shows this
restores validity for out-of-domain samples where one global calibration set fails; concurrently, Hosni
et al. [@hosni2026ad] show explicit applicability-domain calculations flag when uncertainty estimates
are less reliable. Our contribution is complementary and deliberately simpler: rather than per-sample
local calibration (dAD) or AD-as-a-flag (Hosni et al.), we use the applicability domain as the
*stratification axis* for a fixed-strata Mondrian wrapper, and we provide the systematic single-task
benchmark (41 datasets × 3 models × 2 featurizations, continuous AD-resolved coverage, multi-dataset
statistics, and a scaffold-Mondrian negative control) that the quantify-and-recommend question requires.

# 3. Problem Formulation

Let $f$ be a trained regression model and let a calibration set $\{(x_i, y_i)\}_{i=1}^{n}$ give
nonconformity scores (absolute residuals) $s_i = |y_i - f(x_i)|$. Split conformal forms the interval
$f(x) \pm \hat q$ with $\hat q$ the $\lceil (n+1)(1-\alpha)\rceil$-th smallest score; when calibration
and test data are exchangeable this guarantees $\Pr(y \in \text{interval}) \ge 1-\alpha$ marginally.
We stress at the outset that the setting this paper targets — scoring *novel* chemotypes under a
scaffold split — is one where exchangeability between the calibration pool and the novel test compounds
is only approximate; the marginal guarantee therefore holds approximately (empirically, marginal
coverage is near-nominal, Section 6), and our interest is in where that average hides a failure. Writing
$d(x)$ for a molecule's distance to the training set (one minus its mean Tanimoto similarity to the $k$
nearest training molecules), we ask whether $\Pr(y \in \text{interval}\mid d(x) \in B) \approx 1-\alpha$
for strata $B$ ranging from familiar ($d$ small) to novel ($d$ large). A useful guarantee is one that
holds in the novel stratum, not merely on average.

# 4. Method: Applicability-Domain-Conditional Conformal

**Mondrian-AD.** We partition molecules into $m$ strata by quantiles of the calibration-set distance
$d(\cdot)$, and compute a *separate* conformal quantile $\hat q_B$ within each stratum $B$. A test
molecule is assigned the quantile of the stratum its distance falls into, giving the interval
$f(x) \pm \hat q_{B(x)}$. Within each stratum, calibration and test molecules share a similar distance
to training and are approximately exchangeable, so Mondrian conformal targets near-nominal coverage
*conditional on the stratum* [@vovk2005]. Intuitively, novel test molecules are calibrated against other
novel calibration molecules — which share their larger residuals — rather than against the easy bulk, so
their intervals widen to the level they actually require.

Two honest caveats attach to the guarantee. First, we define the stratum edges from the calibration
distances, so the taxonomy is data-dependent rather than a fixed partition; the per-stratum coverage is
therefore *approximate* (it becomes exact in the limit of stable edges) rather than an exact
finite-sample guarantee. Second, the most-novel stratum is open-ended on the right: a test molecule can
be more novel than every calibration molecule in that stratum, so some residual under-coverage of the
very most novel compounds is expected even after conditioning. We report both effects rather than claim
them away.

**Smoothed-AD (shrinkage).** The weakness of a hard taxonomy is that the farthest stratum may be
sparsely populated, making its quantile noisy. We therefore shrink each stratum quantile toward the
pooled quantile, $\tilde q_B = (1-\lambda)\hat q_B + \lambda\,\hat q_{\text{pooled}}$, trading
conditional accuracy for stability in the sparse-novel regime. This shrinkage is a bias–variance
heuristic and, unlike hard Mondrian-AD, carries *no* coverage guarantee; we present it as a practical
default and report where it helps and where it does not.

All methods are post-hoc wrappers: they require only calibration residuals and the cheap
distance-to-training descriptor, no retraining, and apply to any model and featurization. We compare
against four baselines: **marginal** split conformal (the standard); **normalized** conformal with the
textbook additive-regularized score $|r|/(\sigma(x)+\beta)$ and half-width $\hat q\,(\sigma(x)+\beta)$
[@papadopoulos2002; @romano2019; @lei2018distfree], where the difficulty model $\sigma(x)$ is a gradient-boosted regressor
of the absolute residual fit on *out-of-fold training* residuals (a two-fold cross-fit) disjoint from
calibration and test — crucial for validity — and $\beta$ is the median out-of-fold difficulty;
**distance-normalized** conformal [@papadopoulos2011knn], the same normalized form but with the
applicability-domain distance itself as $\sigma(x)=d(x)$, so it conditions on the *same axis* as
Mondrian-AD by normalization rather than stratification; and **Mondrian-scaffold** (Mondrian conditioned
on Bemis–Murcko scaffold cluster rather than the applicability domain). We verified that the normalized
baselines attain near-nominal marginal coverage, confirming they are correctly specified and not
handicapped.

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
applicability-domain strata (defined on the calibration distribution). Our **primary endpoint is
farthest-stratum coverage** — a single, pre-specified bin (the most novel compounds), which is an
unbiased coverage estimate and the scientifically motivated quantity. We also report two secondary,
order-statistic summaries whose selection bias we flag explicitly: *worst-stratum coverage* (the minimum
over the {{N_BINS}} strata, a downward-biased statistic) and the *AD gap* (the maximum stratum-wise
deviation from nominal, an upward-biased statistic); per-stratum interval width is reported alongside.
The primary target is $1-\alpha=0.9$ ($\alpha=0.2$ and $0.05$ as robustness). Each configuration uses
three seeds; the dataset is the unit of analysis.

**Statistics.** Per-method values are averaged over the 18 configurations per dataset (three model
families × two featurizations × three seeds) to one value per dataset, then compared across datasets by
paired Wilcoxon signed-rank tests with bootstrap 95% confidence intervals, and by a Friedman omnibus
with Nemenyi critical-difference analysis over the methods [@demsar2006]. Because 30 of the
{{N_DATASETS}} datasets come from one curated collection (MoleculeACE), the datasets are not fully
independent; we therefore foreground the per-dataset win rate alongside the omnibus $p$-values and report
a per-collection breakdown (Supporting Information). All code is released; every figure and table is
regenerated from machine-readable result files.

# 6. Results: Standard Conformal Gives a False Sense of Security

## 6.1 Marginal coverage is met — on average
At the 90% target, split-conformal intervals attain near-nominal marginal coverage across the
{{N_DATASETS}} datasets (mean {{MARG_COV}}; Table 1), confirming that the marginal guarantee holds
approximately even under scaffold splitting. This is exactly what makes the conditional failure below a
*false sense of security*: nothing in the headline coverage number signals it.

{{TABLE_SUMMARY}}

## 6.2 But novel chemotypes are systematically under-covered
Coverage is not uniform across the applicability domain (Figure @fig-profile). Binning the test set into
{{N_BINS}} strata by distance to the training set, marginal conformal over-covers familiar compounds
(nearest stratum coverage {{MARG_NEAR_COV}}) and under-covers novel ones: coverage in the farthest
(most novel) stratum — our primary endpoint — is {{MARG_FAR_COV}} (mean; median {{MARG_FAR_COV_MED}}),
versus the 0.90 target (Figure @fig-nearfar; Table 1). Coverage falls with distance across the strata
(Figure @fig-profile), consistently across datasets, models, and featurizations. In decision terms, the
nominal 90% interval a practitioner would trust for a genuinely novel candidate covers the truth only
about {{MARG_FAR_PCT}}% of the time. (The worst-of-{{N_BINS}}-strata coverage is lower still,
{{MARG_WORST_COV}}, but as a minimum over noisy bins it is a downward-biased statistic, so we lead with
the pre-specified farthest stratum.)

## 6.3 The standard remedies do not fix the conditional gap
We test the two remedies a practitioner would reach for first. A **learned-difficulty normalized**
score, $|r|/(\sigma(x)+\beta)$ with $\sigma(x)$ a gradient-boosted difficulty model — fit on out-of-fold
training residuals so the procedure stays valid (marginal coverage {{NORM_COV}}) — reallocates width by
predicted difficulty but only modestly improves novel-stratum coverage (farthest-stratum coverage
{{NORM_FAR_COV}}), because the difficulty axis is not the applicability-domain axis and the estimate is
itself unreliable far from training. Conditioning on **Bemis–Murcko scaffold cluster**
(Mondrian-scaffold) leaves the gradient essentially unchanged (farthest-stratum coverage
{{SCAF_FAR_COV}}), because scaffold identity is not the axis along which coverage fails. By contrast, a
**distance-normalized** score [@papadopoulos2011knn], which scales the interval by the
applicability-domain distance itself ($\sigma(x)=d(x)$), *does* repair the gap — it conditions on the
right axis — raising farthest-stratum coverage to {{KNN_FAR_COV}}, and does so efficiently (it is in
fact the tightest of all methods). In short, the conditioning *axis* is what matters: only the
methods that condition on distance-to-training help, and the learned-difficulty and scaffold variants,
despite also being "conditional", do not.

## 6.4 Why learned difficulty fails to adapt to novelty
The failure of the normalized score has a simple, measurable cause. A difficulty model $\sigma(x)$ is
trained on in-domain residuals, so it learns what makes *familiar* compounds hard — and that signal does
not extrapolate. Stratifying the test set by applicability-domain distance and measuring how each
quantity grows from the familiar to the novel stratum (Figure @fig-mechanism), the true prediction error
grows by about {{RESID_GROWTH}}% (median over datasets), but the learned difficulty $\sigma(x)$ grows by
only {{SIGMA_GROWTH}}% — it barely registers that novel compounds are harder, so it cannot widen their
intervals enough. Distance to the training set, by contrast, grows by {{DIST_GROWTH}}% out of domain by
construction, more than covering the error growth; that is precisely why conditioning on it restores
coverage while conditioning on a learned difficulty does not. The lesson generalizes beyond conformal
prediction: any heteroscedastic uncertainty estimate learned in-domain will tend to under-react to
novelty, and novelty must be measured, not learned.

# 7. Results: Conditioning on the Applicability Domain Restores Coverage

## 7.1 The axis is what matters
Conditioning calibration on the applicability domain repairs the gap. Taking our Mondrian wrapper
(Mondrian-AD) as a concrete instance, coverage in the farthest (most novel) stratum — our primary
endpoint — rises from {{MARG_FAR_COV}} (marginal) to {{MOND_FAR_COV}} (Mondrian-AD), an {{FAR_GAIN}}-point
mean gain to at or slightly above the 0.90 target (paired Wilcoxon $p={{FAR_P}}$; 95% CI
[{{FAR_CI_LO}}, {{FAR_CI_HI}}] points; improved on {{FAR_FRAC}}% of datasets), closing the
novel-chemotype gap (Figure @fig-nearfar, @fig-profile, @fig-perdataset). Marginal coverage does not
drop; it rises slightly to {{MOND_COV}} (the method becomes marginally conservative). The downward-biased
worst-stratum coverage also improves modestly ({{MARG_WORST_COV}}$\to${{MOND_WORST_COV}}, mean gain
{{WORST_GAIN}} points, $p={{WORST_P}}$); residual under-coverage remains on the hardest small datasets
(e.g. FreeSolv), consistent with the open-top-stratum caveat.

**Several applicability-domain-conditional methods work, and they are statistically indistinguishable.**
The crucial comparison is not Mondrian-AD against marginal but against the other ways of conditioning.
The established distance-normalized conformal [@papadopoulos2011knn] reaches farthest-stratum coverage
{{KNN_FAR_COV}} and the shrinkage variant smoothed-AD {{SMOOTH_FAR_COV}}, both close to Mondrian-AD's
{{MOND_FAR_COV}}; a Friedman omnibus over the {{N_METHODS}} strategies on farthest-stratum coverage
rejects equality ($p={{FRIEDMAN_P}}$, Nemenyi critical difference {{NEMENYI_CD}}), but the three
applicability-domain-conditional methods (Mondrian-AD, smoothed-AD, distance-normalized) fall **within one
critical difference of each other** — they are not statistically separable — while marginal, the
scaffold-conditional variant, and the learned-difficulty normalized score are all more than a critical
difference worse. In other words, the effect is driven by the conditioning *axis*, not by the particular
mechanism: distance-normalization and stratified calibration are interchangeable on this benchmark.
Distance-normalized is in fact the tightest method (Table 1); the Mondrian variants we describe are not
sharper, but they are simple, interpretable, and carry an approximate per-stratum guarantee that the
normalized score does not. We therefore do not claim a winner among the applicability-domain methods; we
claim the axis.

**Is the improvement an evaluation artifact?** Because Mondrian-AD is calibrated and evaluated on the
same applicability-domain strata, one might worry its per-stratum restoration is guaranteed by
construction. That the distance-normalized baseline — which conditions on the same axis by scaling the
interval by $d(x)$, with no per-bin calibration — achieves the same far-stratum restoration argues
otherwise: the gain reflects genuine distance-conditioning, not an artifact of scoring on the calibration
partition, and the effect is smooth across all five strata (Figure @fig-profile), not a discrete jump in
one calibrated bin. We are explicit (Section 10) that exact per-stratum validity is not claimed: the
strata are data-dependent and the most novel stratum is open-ended.

## 7.2 The cost, and the trade-off
Reliability for novel compounds is not free, but the cost is small and well-placed. Averaged over all
compounds the applicability-domain-conditional methods barely change interval width (Mondrian-AD
relative width {{MOND_RELWIDTH}} vs {{MARG_RELWIDTH}} for marginal, {{WIDTH_COST}}%; distance-normalized
is actually *tighter* than marginal), because the extra width is concentrated — correctly — in the novel
strata. In the farthest stratum the typical (median) interval widens by about {{FARWIDTH_MED_PCT}}%
relative to marginal (the mean widens more, inflated by a few pathological datasets; Figure
@fig-tradeoff), while familiar-compound intervals are barely affected. This is the intended behavior —
larger, honest uncertainty exactly where the model extrapolates — but it is a real cost: for the most
novel compounds the interval is correspondingly wider, so the methods buy trustworthy coverage, not
sharper decisions, in that regime (Section 8).

## 7.3 Robustness
The finding and the fix are stable across conditions (Supporting Information): both featurizations, all
three model families, the random split (where intervals are narrower but the applicability-domain
coverage gradient persists and the AD-conditional methods still restore the far stratum), and the
$\alpha=0.2$ and $\alpha=0.05$ targets. We fix the stratum count at {{N_BINS}} and the shrinkage at $\lambda=0.3$
throughout; the hard-binned (Mondrian-AD, $\lambda=0$) and shrunk (smoothed-AD, $\lambda=0.3$) variants
bracket the shrinkage choice, and a sensitivity analysis over the number of strata (Supporting
Information) leaves the qualitative conclusion unchanged, with more strata giving finer conditioning at
the cost of noisier sparse-stratum quantiles.

# 8. Downstream: Interval Width as an Honest Novelty Signal
The practical value of conditioning is not sharper predictions on novel compounds — no method can
manufacture information that is not there — but an uncertainty estimate a practitioner can *act on* to
triage candidates, and here the applicability-domain-conditional methods differ qualitatively from
standard conformal. Marginal conformal assigns **every** compound the same interval width, so width
carries no information about which predictions to trust: abstaining on the widest-interval compounds
flags novel chemotypes no better than chance (the most novel, far-AD compounds make up
{{MARG_ABSTAIN_FAR}} of the widest 30% of intervals, essentially the base rate {{ABSTAIN_BASE}}). The
applicability-domain-conditional methods — Mondrian-AD and, by construction, the distance-normalized
baseline — instead make interval width track distance to the training set (for Mondrian-AD, Spearman
$\rho={{MOND_WIDTH_RHO}}$ between width and applicability-domain stratum), so width becomes an honest
novelty flag: the widest intervals are roughly twice enriched for novel compounds ({{MOND_ABSTAIN_FAR}}
vs the {{ABSTAIN_BASE}} base rate; Figure @fig-selective). The decision payoff is therefore a trustworthy
triage signal — scrutinize or assay the wide-interval candidates, rely on the
narrow ones — that standard conformal, with its uniform width, simply cannot provide, together with
coverage on the novel compounds that is actually at its stated level.

# 9. Discussion
Our results turn a known theoretical caveat into a concrete, measured, and fixable problem for molecular
modeling, and situate it against recent applicability-domain-aware conformal work. The marginal coverage
guarantee that makes conformal prediction attractive is, by itself, insufficient for the use case the
field cares about most: trustworthy uncertainty on novel chemotypes. Standard conformal intervals
quietly under-cover exactly those compounds, giving a false sense of security in prospective settings.
Our central, controlled finding is that *which axis a conformal method conditions on is what matters*:
a learned-difficulty normalized score and a scaffold-conditional variant do not repair the gap, whereas
every method that conditions on distance to the training set does — the established distance-normalized
conformal [@papadopoulos2011knn], the per-sample dynamic calibration of Oršolić and Šmuc
[@orsolic2023dad], and the simple fixed-strata Mondrian wrapper we describe. On this benchmark these
applicability-domain-conditional methods are statistically indistinguishable; the distance-normalized
baseline is in fact the tightest, and we do not claim our wrapper improves on it. The Mondrian variants
are offered as simple, interpretable alternatives that additionally carry an approximate per-stratum
guarantee. The contribution is therefore the *benchmark and the isolation of the axis*, plus a clear
practical recommendation: report and calibrate conformal coverage conditional on the applicability
domain — by whichever of these equivalent routes is convenient — rather than marginally or on an
unrelated axis. We are explicit that the restored per-stratum coverage is approximate (the strata are
data-dependent, the most-novel stratum is open-ended) and that the gain is in *honest* uncertainty, not
sharper predictions.

# 10. Limitations
(i) The per-stratum coverage is *approximate*, not an exact finite-sample guarantee: the stratum edges
are estimated from the calibration distances, the most-novel stratum is open-ended on the right, and the
shrinkage (smoothed-AD) variant carries no guarantee at all. Exact object-conditional coverage is
impossible in finite samples [@vovk2012conditional]. (ii) Under a scaffold split the calibration pool
and the novel test compounds are only approximately exchangeable, so the marginal guarantee too is
approximate; our marginal coverage is an empirical outcome (near-nominal), and a genuine temporal/time
split — where exchangeability is more strongly violated — is important future work. (iii) We use one
standard applicability-domain descriptor (ECFP Tanimoto distance); alternative descriptors
(descriptor-space distance, leverage, embedding density) and whether the finding and fix transfer to
them are left to future work. (iv) We study regression with CPU-feasible fingerprint/descriptor learners;
pretrained graph or language models, and conformalized quantile regression with strong conditional base
models, are natural extensions, as is classification. (v) Thirty of the {{N_DATASETS}} datasets come
from one curated collection (MoleculeACE), so the datasets are not fully independent; we therefore report
the per-dataset win rate and a per-collection breakdown alongside the omnibus statistics, and treat the
latter as indicative rather than exact.

# 11. Conclusion
Standard conformal prediction gives molecular property models a coverage guarantee that holds on
average but fails where it matters most — on novel chemotypes — a false sense of security in exactly the
prospective, out-of-domain setting that drug discovery depends on. Across {{N_DATASETS}} datasets, three
model families, and two featurizations we quantified this systematic under-coverage and, in a controlled
head-to-head, isolated the conditioning axis that repairs it: learned-difficulty and scaffold-conditional
variants do not help, while every method that conditions on distance to the training set — the existing
distance-normalized conformal and the simple applicability-domain Mondrian wrapper we describe, which are
statistically indistinguishable here — raises novel-compound coverage to near nominal, approximately and
per stratum, with no retraining, at a modest and well-placed width cost. We release an open benchmark and
implementation and recommend conditioning conformal calibration on the applicability domain as standard
practice for trustworthy molecular uncertainty.

# Data and Software Availability
All datasets used in this study are publicly available and require no credentials: the MoleculeACE
activity-cliff benchmark collection (https://github.com/molML/MoleculeACE) [@vantilborg2022], the
Therapeutics Data Commons ADMET regression tasks (https://tdcommons.ai/) [@tdc2021], and the MoleculeNet
physicochemical sets ESOL and FreeSolv (https://moleculenet.org/) [@moleculenet2018]. The complete source
code, environment specifications, deterministic seeds, and scripts that regenerate every figure and table
from machine-readable result files will be made available by the authors upon reasonable request. Models
use scikit-learn [@sklearn2011] and LightGBM [@lightgbm2017]; fingerprints and descriptors use RDKit
[@rdkit].

# Acknowledgments
The authors received no specific funding for this work.

# Author Contributions
P. N. Bhargav Teja conceived and carried out the study — the implementation, experiments, analysis, and
writing. Ragavan K supervised the work. Both authors reviewed and approved the manuscript.

# Competing Interests
The authors declare no competing financial interest.
