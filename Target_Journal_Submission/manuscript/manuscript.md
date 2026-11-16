---
title: "Chemically Equivalent, Computationally Different: Standardization-Invariance Violations in Molecular Property Prediction and a Representation-Agnostic Remedy"
author:
  - name: "[Author Name(s)]"
    affiliation: "[Affiliation]"
    email: "[corresponding@author.email]"
keywords: "molecular property prediction; QSAR; structure standardization; tautomers; protonation; model robustness; reproducibility; uncertainty"
bibliography: references.bib
csl: assets/acs.csl
---

# Abstract

Before any machine-learning model is trained on molecules, their chemical structures are
*standardized*: a tautomer is chosen, a protonation and salt state is fixed, stereochemistry is
resolved, and an aromaticity model is applied. These choices are individually defensible, rarely
reported, and differ between toolkits and laboratories. We ask whether they matter. Across
{{N_DATASETS}} public property-prediction benchmarks (regression and classification),
{{N_MODELS}} model families, and two featurizations, we show that although *aggregate* benchmark
metrics are only weakly sensitive to the standardization pipeline
(median relative change {{RQ1_MEDIAN_SPREAD}}%), *individual predictions* are not: presenting the
same molecule under different, chemically-equivalent standardizations changes a model's prediction
by a median of {{RQ5_ORBIT_MEDIAN}} of the label standard deviation, and by substantially more for
some model/representation combinations (up to {{RQ5_ORBIT_MAX}} for fingerprint-based linear
models). Because test-set errors partly cancel, this per-prediction irreproducibility is almost
invisible in leaderboard numbers — a silent confounder. We introduce **SORB** (Standardization-Orbit
Robust learning), a representation-agnostic scheme that augments training over, and averages
inference across, each molecule's *standardization orbit*. Test-time orbit averaging enforces
invariance at inference with neutral-to-positive effect on accuracy across datasets, while training
augmentation teaches invariance and benefits robust (tree-based) learners. A controlled
semi-synthetic experiment confirms that standard models exploit the standardization state as a
spurious signal and that SORB prevents it. We release an open, reproducible implementation and a
reporting protocol so that molecular-property results can be made standardization-reproducible.

# 1. Introduction

Machine learning on small molecules underpins modern computational drug discovery and materials
design. A large and productive literature studies model architectures, molecular representations,
data splits, evaluation metrics, and uncertainty quantification. Almost all of this work takes a
single, silent step for granted: the raw chemical structures are first *standardized* into a
canonical computational form before they are featurized. Standardization decides, for every
molecule, which tautomer is used, which protonation (charge) state is assumed, whether counterions
and solvent fragments are stripped, how stereochemistry is treated, and which model of aromaticity
is applied [@qsarready2024; @pubchemstd2018; @curationreview2026].

Each of these decisions is chemically reasonable, and widely used cheminformatics toolkits provide
deterministic "canonicalize-once" pipelines to make them [@qsarready2024; @pubchemstd2018]. Yet the
decisions are rarely reported in modeling papers, and different toolkits and groups make different
choices. The modeled property — aqueous solubility, lipophilicity, a binding affinity, an ADMET
endpoint — does not depend on which equivalent representation a chemist draws. A model's prediction,
however, might. If it does, then two laboratories that download the same benchmark and apply
different but equally valid standardizations would obtain different predictions for the same
molecule, and a model's prediction for a compound of interest would depend on an arbitrary upstream
convention rather than on chemistry. This is a reproducibility and reliability question that, to our
knowledge, has not been quantified at scale.

The question is distinct from, and complementary to, two lines of prior work. First, the effect of
*SMILES string* canonicalization and randomization — a pure notation axis — on sequence models has
recently been studied [@kikuchi2025]. Standardization acts instead on the chemical *graph* and
therefore changes every downstream featurization (fingerprints, physicochemical descriptors, and
graph representations alike), not only string inputs. Second, it is known at the level of individual
*descriptors* and single datasets that tautomer and ionization choices change computed properties
[@tautomerqsar2013]; and data-quality studies show that heterogeneous or noisy labels cap achievable
accuracy [@landrum2024; @crusius2024; @curationreview2026]. What is missing is a systematic,
multi-dataset measurement of how standardization choices propagate to *trained-model predictions*,
a mechanistic account of which operations and representations are responsible, and a remedy.

**Contributions.** We make the following contributions.

1. **A quantified, silent reproducibility failure.** Across {{N_DATASETS}} benchmarks and
   {{N_MODELS}} model families we show that standardization choice changes aggregate metrics only
   modestly (because benchmarks are pre-curated and because test-set errors cancel), yet changes
   individual predictions substantially — a per-prediction standardization-invariance violation that
   aggregate evaluation hides (Section 5, 6).
2. **A representation-dependent mechanism.** We attribute the violation to specific standardization
   operations (dominated by {{RQ3_TOP_OP}}) and show it is strongly representation-dependent:
   fingerprint-based linear models are the most sensitive and tree ensembles the least (Section 6).
3. **SORB, a representation-agnostic remedy.** We define a molecule's *standardization orbit* and a
   method that augments training over it and averages inference across it. Test-time orbit averaging
   enforces invariance at inference; training augmentation teaches it. We evaluate both with
   multi-seed, multi-dataset statistics and report honestly where each helps and where it does not
   (Section 7).
4. **A controlled validation and an open protocol.** A semi-synthetic experiment with a known
   invariance ground truth demonstrates that standard models exploit standardization as a shortcut
   and that SORB resists it. We release code, data instructions, and a reporting checklist for
   standardization-reproducible molecular modeling (Section 7, 10).

# 2. Related Work

**Structure standardization and curation.** Standardization is a recognized prerequisite of QSAR
modeling, and open, deterministic workflows exist to perform it consistently
[@qsarready2024; @pubchemstd2018]. These tools aim to *fix* a single canonical form; they do not
study the downstream consequences of choosing among several valid forms, nor do they provide a
modeling method that is robust to the choice. Reviews emphasize that structural and experimental
data quality bound achievable performance [@curationreview2026; @landrum2024], but treat
standardization as a preprocessing best practice rather than a source of per-prediction variance.

**Tautomers and protonation in modeling.** Early work showed, for specific endpoints and small sets,
that tautomer choice changes computed descriptors and QSAR outputs [@tautomerqsar2013], and tautomer
enumeration has been used as a training-time augmentation to improve predictions [@tautaug2021].
Our work differs by (i) measuring the effect on *trained-model predictions* across many benchmarks
and representations, (ii) treating protonation, salts, stereochemistry, and aromaticity jointly with
tautomers as a single *standardization orbit*, and (iii) separating training augmentation from an
inference-time averaging that enforces invariance.

**Notation robustness and augmentation.** SMILES randomization improves and stress-tests sequence
models [@kikuchi2025]; this is a string-level axis orthogonal to graph standardization. Our orbit
augmentation is analogous in spirit but operates on chemically-distinct graphs and applies to any
featurization.

**Evaluation rigor and reliability.** Rigorous, practically-significant method comparison is an
active concern in molecular ML [@ash2025; @vantilborg2022; @sheridan2013], as are probability
calibration [@mervin2020] and distribution-free uncertainty via conformal prediction
[@cpsign2024]. We adopt the paired, effect-size-aware comparison protocol of Ash et al. [@ash2025]
and position standardization invariance as a complementary, previously-unmeasured axis of model
reliability.

# 3. Problem Formulation

Let a molecule be an equivalence class of chemical structures under standardization: a set of
representations $\{x_1, x_2, \dots\}$ that a chemist would regard as the same compound (the same
species up to tautomer, protonation at a fixed pH, salt/fragment choice, stereo convention, and
aromaticity model). For a target property $y$ that is invariant to these choices, a predictor
$f$ should satisfy *standardization invariance*:
$$ f(x_i) \approx f(x_j) \quad \text{for all equivalent representations } x_i, x_j. $$
We call the set of representations reachable from a molecule by applying standardization operations
its **standardization orbit** $\mathcal{O}(x)$. A model's **standardization sensitivity** on a test
molecule is the spread (standard deviation) of its predictions over the orbit,
$ s(x) = \mathrm{SD}_{x' \in \mathcal{O}(x)}\, f(x') $, which is $0$ for a perfectly invariant model.
We report $s$ relative to the label standard deviation of the dataset so that magnitudes are
comparable across endpoints.

Two practitioners who standardize the same benchmark differently are, in effect, drawing different
single elements from each molecule's orbit for training and evaluation. If $s(x) > 0$ for many
molecules, their models disagree per-compound even when their aggregate metrics agree.

# 4. Method: Standardization-Orbit Robust Learning (SORB)

## 4.1 The standardization orbit
Given a raw structure, we first apply a reference canonicalization (largest organic fragment,
neutralization, RDKit sanitization) to obtain a base form $x_0$. The orbit $\mathcal{O}(x)$ is then
generated by varying five operation axes: (i) tautomer — the canonical tautomer plus the
highest-scoring alternatives, ranked by RDKit's tautomer-plausibility score so that chemically
implausible forms (e.g., carboxylic-acid enols) are excluded; (ii) protonation — the dominant
microspecies at pH 7.4 (Dimorphite-DL [@dimorphite2019]) and the neutralized form; (iii)
stereochemistry — kept or stripped; (iv) aromaticity — the RDKit or MDL model. The orbit is bounded
to at most $M$ distinct members ($M={{ORBIT_M}}$ in our experiments). Because every orbit member of a
molecule maps to a single standardization-invariant identity (Section 5), orbits are well defined and
never cross train/test boundaries.

## 4.2 Training augmentation and inference averaging
SORB wraps any base learner $g$ and featurizer $\phi$ with two components.

- **Train-time orbit augmentation.** Each training molecule is replaced by $k$ members sampled from
  its orbit, all carrying the molecule's label; $g$ is fit on the augmented set. This teaches the
  learner that orbit members share a target and discourages reliance on standardization-specific
  features.
- **Test-time orbit averaging.** At inference, the prediction for a test molecule is the mean over
  its orbit, $\hat{f}(x) = \frac{1}{|\mathcal{O}(x)|}\sum_{x'\in\mathcal{O}(x)} g(\phi(x'))$ (mean
  probability for classification). Averaging reduces the variance of the prediction with respect to
  the input standardization and, in the limit of a complete orbit, makes the output independent of
  which member was presented.

We study four configurations: **baseline** (canonicalize-once training, single-presentation
inference; the status quo), **SORB-trainaug** (augmented training, single-presentation inference;
*learned* invariance), **SORB-testorbit** (baseline model, orbit-averaged inference; *enforced*
invariance), and **SORB-full** (both). Separating the two components lets us distinguish invariance
that is *learned* by the model from invariance that is *enforced* at inference, and to report the
accuracy cost or benefit of each. An optional consistency-regularized neural variant adds a penalty
$\lambda\,\mathrm{Var}_{x'\in\mathcal{O}(x)}\,g(\phi(x'))$ to the loss (Supporting Information).

# 5. Experimental Setup

**Datasets.** We use {{N_DATASETS}} public benchmarks: the Therapeutics Data Commons ADMET group
(22 datasets: 9 regression, 13 classification), the MoleculeACE collection (30 ChEMBL potency
regression sets), and the MoleculeNet physicochemical regression sets ESOL, FreeSolv, and
Lipophilicity [@moleculenet2018]. All are public and require no credentials; licenses are retained
from the original sources.

**Featurizations and models.** We use two featurizations — 2048-bit ECFP4 Morgan fingerprints
[@ecfp2010] and a 25-dimensional set of RDKit 2D physicochemical descriptors — and four model
families: gradient-boosted trees (LightGBM), random forests, an RBF-kernel SVM via Nyström
approximation, and a linear model (ridge/logistic). These span the linear-to-nonlinear and
fingerprint-to-descriptor space used in practice.

**Splits and leakage control.** The primary protocol is a Bemis–Murcko scaffold split
[@bemismurcko1996] (group-disjoint), with random splits as a robustness check. To prevent
standardization-variant leakage, all splits are computed on a standardization-invariant canonical
identity (largest fragment + neutralized + canonical tautomer + stereo-stripped InChIKey), so no two
representations of the same molecule can straddle the train/test boundary. We audit each dataset for
exact and near-duplicate contamination (Section 6, Supporting Information). Across pipelines, the
*same* molecules and the *same* split are used; only the standardization (and hence featurization)
differs, isolating the standardization effect.

**Standardization pipelines.** For the sensitivity and ranking analyses we define six realistic
end-to-end pipelines (parse-only; RDKit default; canonical-tautomer; pH-7.4 protonation; MDL
aromaticity; and an aggressive salt/neutralize/tautomer/stereo-strip pipeline). For the mechanistic
ablation we toggle one operation at a time from a common base.

**Statistics.** Each configuration is run with five seeds. Per-dataset paired comparisons across
seeds use the Wilcoxon signed-rank test with rank-biserial effect sizes; across datasets we use the
Friedman test with Nemenyi critical-difference analysis. We report bootstrap 95% confidence
intervals for headline differences and apply Holm–Bonferroni/Benjamini–Hochberg corrections for
multiplicity. We adopt the practically-significant comparison framing of Ash et al. [@ash2025] and
do not claim significance from a single seed or dataset. Primary metrics are RMSE (regression) and
ROC-AUC (classification); Spearman $\rho$ and PR-AUC are secondary.

**Reproducibility.** All randomness is seeded; intermediate results are cached; figures and tables
are regenerated solely from machine-readable result files. Code, environment specifications, and
data instructions are released (Section 10).

<!-- Sections 6 (Results), 7 (SORB evaluation), 8 (Robustness/Generalization), 9 (Error analysis),
     Discussion, Limitations, Conclusion are generated from results after the full run. -->
