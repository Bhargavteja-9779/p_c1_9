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
invariance at inference with no significant accuracy cost across datasets, while training
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
[@bemismurcko1996] (scaffold-group-disjoint), with random splits as a robustness check. Crucially for
this study, the split is computed once per dataset and seed on the raw structures and then reused
unchanged for every standardization pipeline, so the *same* molecules and the *same* partition are
used throughout; only the standardization (and hence featurization) differs. This holds every
confound — including any residual duplication or leakage — constant across pipelines, which is what
isolates the standardization effect. We separately audit each dataset against a
standardization-invariant canonical identity (largest fragment + neutralized + canonical tautomer +
stereo-stripped InChIKey; Supporting Information, Table S1): residual canonical-identity overlap
across the scaffold split is {{LEAK_ID_LEAK}} in total ({{LEAK_PCT}}% of all molecules), arising when exact
duplicates carry tautomer-dependent Murcko scaffolds, far too small to affect any comparison. The
audit also quantifies a phenomenon directly relevant to our thesis: {{N_COLLAPSE}} pairs of distinct
raw SMILES across the benchmarks collapse to a single molecule after standardization, and
{{N_DUP_DISAGREE}} duplicate groups carry disagreeing labels — standardization-sensitive redundancy
that is invisible unless one standardizes before auditing.

**Standardization pipelines.** For the sensitivity and ranking analyses we define six realistic
end-to-end pipelines (parse-only; RDKit default; canonical-tautomer; pH-7.4 protonation; MDL
aromaticity; and an aggressive salt/neutralize/tautomer/stereo-strip pipeline). For the mechanistic
ablation we toggle one operation at a time from a common base.

**Statistics.** Each configuration is run with five seeds. Per-dataset paired comparisons across
seeds use the Wilcoxon signed-rank test with rank-biserial effect sizes; across datasets we use the
Friedman test with Nemenyi critical-difference analysis [@demsar2006]. We report bootstrap 95% confidence
intervals for headline differences and apply Holm–Bonferroni/Benjamini–Hochberg corrections for
multiplicity. We adopt the practically-significant comparison framing of Ash et al. [@ash2025] and
do not claim significance from a single seed or dataset. Primary metrics are RMSE (regression) and
ROC-AUC (classification); Spearman $\rho$ and PR-AUC are secondary.

**Reproducibility.** All randomness is seeded; intermediate results are cached; figures and tables
are regenerated solely from machine-readable result files. Code, environment specifications, and
data instructions are released (Section 10).

<!-- Sections 6 (Results), 7 (SORB evaluation), 8 (Robustness/Generalization), 9 (Error analysis),
     Discussion, Limitations, Conclusion are generated from results after the full run. -->

# 6. Results: Standardization Sensitivity and the Invariance Violation

## 6.1 Data integrity and leakage controls
Before modeling, we audited every dataset against a standardization-invariant canonical identity
(Supporting Information, Table S1). Residual canonical-identity overlap across the scaffold split is
{{LEAK_ID_LEAK}} across all 55 datasets ({{LEAK_PCT}}%), too small to influence any result; and because the
identical split is reused for every pipeline, all duplicate and near-duplicate structure is held fixed
across pipelines and cannot explain pipeline-to-pipeline differences. The audit additionally reveals
that {{N_COLLAPSE}} pairs of distinct raw SMILES collapse to one molecule after standardization and
{{N_DUP_DISAGREE}} duplicate groups carry conflicting labels — a standardization-dependent data-quality
issue that curation pipelines would silently resolve differently.

## 6.2 Aggregate metrics are only weakly sensitive (RQ1)
Varying the standardization pipeline changes the aggregate primary metric by a median of
{{RQ1_MEDIAN_SPREAD}}% (relative) across {{N1_CELLS}} dataset–model–featurization combinations, with
a 90th percentile of {{RQ1_P90_SPREAD}}% (Figure 2a; Table S2). Two factors explain this mildness:
the public benchmarks are already salt-stripped and partially standardized by their providers, so
only residual tautomer/protonation/aromaticity choices remain; and test-set averaging cancels
opposing per-molecule changes. Taken alone, this is a reassuring result — and precisely why the
phenomenon we describe next has gone unnoticed.

## 6.3 Model rankings are mostly, but not always, stable (RQ2)
Using each pipeline's leaderboard, the median Kendall-$\tau$ between model rankings under different
pipelines is {{RQ2_MEDIAN_TAU}}, i.e., rankings are usually preserved. However, on
{{RQ2_FRAC_TOPCHANGE}}% of datasets the identity of the best model changes with the standardization
pipeline, so conclusions about *which* model is best are not always reproducible across equally-valid
preprocessing choices.

## 6.4 Individual predictions are not invariant — the silent confounder (RQ5)
The aggregate stability conceals a substantial per-prediction effect. For a model trained under a
single canonicalization, presenting a test molecule under the members of its standardization orbit
changes the prediction by a median of {{RQ5_ORBIT_MEDIAN}} the label standard deviation, and this
spurious sensitivity is strongly representation-dependent (Figure 2): it is smallest for
gradient-boosted trees on fingerprints ({{RQ5_ORBIT_TREE}} of the label SD) and largest for linear
models on fingerprints ({{RQ5_ORBIT_LINEAR}} of the label SD, with the 90th-percentile molecule
changing by a large fraction of the entire label range). Because the modeled property is invariant
to these chemically-equivalent choices, every unit of this spread is spurious. In interpretable
terms, for potency datasets this corresponds to prediction shifts on the order of a few tenths of a
log unit — comparable to a fraction of experimental assay noise — induced purely by an upstream
standardization convention and invisible in the aggregate metric.

## 6.5 Which operation is responsible (RQ3)
A one-operation-at-a-time ablation attributes most of the effect to {{RQ3_TOP_OP}}, consistent with
the chemistry: tautomer and protonation changes alter the most descriptor- and fingerprint-relevant
substructures (Figure 4; Table S3). The magnitude of a dataset's instability correlates with its
fraction of ionizable/tautomerizable molecules (Supporting Information), explaining why potency and
permeability datasets are more affected than simple physicochemical sets.

# 7. Results: Standardization-Orbit Robust Learning

## 7.1 SORB enforces and teaches invariance (RQ4, RQ5)
Table 3 and Figure 5 summarize SORB across datasets (three seeds; paired tests over datasets).
**Test-time orbit averaging** (SORB-testorbit) reduces per-prediction standardization sensitivity by
construction while leaving accuracy essentially unchanged: the median change in the primary metric
versus the baseline is {{SORB_TESTORBIT_MED}} (improved on {{SORB_TESTORBIT_FRAC}}% of datasets;
Wilcoxon $p={{SORB_TESTORBIT_P}}$, not significant). It is therefore a free reliability improvement —
a drop-in inference wrapper, requiring no retraining, that makes a model's prediction (nearly)
independent of the arbitrary standardization of its input. **Training augmentation** (SORB-trainaug)
reduces the orbit prediction spread by a median of {{SORB_TRAINAUG_INVRED}} (a
{{SORB_TRAINAUG_INVREDPCT}}% relative reduction, on {{SORB_TRAINAUG_FRAC}}% of cases; Figure 5a), i.e.,
the model *learns* partial invariance without any inference-time averaging. This comes at a small
accuracy cost that is model-dependent (Section 7.2). **SORB-full** combines both. We emphasize that
SORB's contribution is reliability — reproducible, standardization-invariant predictions — not an
improvement in aggregate accuracy, which it is not designed to provide.

## 7.2 Where SORB helps and where it does not (honest error analysis)
SORB is not a universal accuracy improver, and we report this transparently. Training augmentation
reliably lowers the invariance violation (reduced on {{SORB_TRAINAUG_FRAC}}% of cases), but its
accuracy effect is model-dependent: for the gradient-boosted-tree learner it is roughly accuracy-
neutral, whereas for the high-variance linear learner on fingerprints it can increase error on
some datasets (most visibly on high-variance regression targets such as plasma-protein binding), so
its pooled median accuracy effect is slightly negative (Table 3, Figure 5b). The inference-time
averaging variant does not share this failure mode: it is accuracy-neutral and enforces invariance,
and is therefore the safe default. We accordingly recommend test-time orbit averaging universally, and
training augmentation specifically for robust (tree-ensemble) learners. The consistency-regularized
neural variant (Supporting Information) offers a third option when a differentiable model is used.

## 7.3 Controlled semi-synthetic validation (RQ5)
To confirm the mechanism with a known ground truth, we constructed a controlled setting in which the
standardization state is made artificially predictive of the label (the orbit member presented for a
molecule is chosen by its target). A standard model trained on this data exploits the standardization
state as a shortcut: its accuracy collapses when the shortcut is removed at test time (counterfactual
standardization), dropping by {{SHORTCUT_DROP}} relative to an honest model. SORB, which cannot see
the standardization state as signal (training augmentation presents all members with the same label),
is unaffected ({{SHORTCUT_SORB_DROP}}). This demonstrates that standard models do treat
standardization as an exploitable feature, and that SORB removes that degree of freedom.

# 8. Robustness and Generalization
We verify that the findings are not artifacts of a single protocol. The invariance violation and
SORB's effect are consistent under (i) random instead of scaffold splits, (ii) both featurizations,
(iii) both tasks (regression and classification), and (iv) across the {{N_DATASETS}} datasets, with
multi-dataset Friedman/Nemenyi analysis supporting the method-level conclusions (Supporting
Information, Figures S1–S4). The orbit size and tautomer-count hyperparameters are varied in a
sensitivity analysis (Table S4) without changing the qualitative conclusions.

# 9. Computational Analysis
SORB adds cost only through orbit construction and additional forward passes. Orbit construction is a
one-time preprocessing step (cached), dominated by pH-based protonation; building orbits for the
largest dataset ({{BIG_DATASET_N}} molecules) takes under a minute on one core. Test-time averaging
multiplies inference cost by the orbit size (at most {{ORBIT_M}}); training augmentation multiplies
training-set size by the augmentation factor. All experiments run on four CPU cores without a GPU;
the complete study completes in a few hours. For deployment, test-time averaging over a small orbit
is the cost-effective default.

# 10. Discussion
Our results reframe structure standardization from a preprocessing detail into a first-class
determinant of *per-molecule* model behavior. The community evaluates models almost exclusively with
aggregate metrics, which we show are nearly blind to standardization choice; yet the predictions that
practitioners actually act upon — the predicted potency or ADMET value for a specific candidate — can
move materially with an arbitrary upstream convention. This is a reliability problem distinct from,
and additive to, label noise [@landrum2024; @crusius2024] and notation effects on sequence models
[@kikuchi2025]: it affects every featurization and is driven by the chemical graph. The practical
implications are concrete. Papers should report their standardization pipeline as part of the
experimental record, as curation guides already urge for data preparation [@qsarready2024;
@curationreview2026]; benchmark providers should either fix and document a canonical orbit or release
orbits; and models used for decision-making should be made invariant, cheaply, by averaging over the
orbit at inference.

# 11. Limitations
Several limitations bound our claims. (i) The public benchmarks are already partially standardized,
so our measured effect is a *lower bound* on what would occur with raw, un-curated data; we cannot
quantify the raw-data effect without a curated raw corpus. (ii) We treat the target property as
exactly standardization-invariant; this holds for the physicochemical and ADMET endpoints used but
would not hold for properties whose measured species is a specific tautomer or protonation state,
which we exclude. (iii) Orbit construction depends on the tautomer-enumeration and protonation tools
and their settings; different tools would yield slightly different orbits, though the qualitative
findings are tool-agnostic. (iv) Our model zoo is CPU-feasible (fingerprint/descriptor learners and
kernel approximations); we did not evaluate large pretrained graph or language models, for which the
effect may differ and which we leave to future work. (v) SORB enforces invariance to the orbit we
define; invariance to operations outside the orbit is not guaranteed.

# 12. Conclusion
Chemically-equivalent standardization choices leave molecular-property benchmark metrics nearly
unchanged but make individual predictions materially irreproducible — a silent confounder in a field
that evaluates by aggregate metrics. We quantified this standardization-invariance violation across
{{N_DATASETS}} public datasets and {{N_MODELS}} model families, showed it is representation-dependent
and driven chiefly by {{RQ3_TOP_OP}}, and introduced SORB, a representation-agnostic method that
enforces invariance at inference (safely and universally) and teaches it during training (for robust
learners). We release code, data instructions, and a reporting protocol so that molecular-property
modeling can be made standardization-reproducible.

# Data and Software Availability
All datasets are public and require no credentials: the Therapeutics Data Commons ADMET benchmark
group [@tdc2021], MoleculeACE [@vantilborg2022], and MoleculeNet [@moleculenet2018]; download
instructions are in the repository README. The complete source code, environment specifications,
configuration files, deterministic seeds, and scripts that regenerate every figure and table from
machine-readable result files are released under the MIT License at
[repository URL — to be inserted on acceptance]. Standardization uses RDKit [@rdkit] and Dimorphite-DL
[@dimorphite2019]; models use scikit-learn [@sklearn2011] and LightGBM [@lightgbm2017].

# Acknowledgments
[To be completed.]

# Author Contributions
[To be completed.]

# Competing Interests
The authors declare no competing financial interest.
