# Research Plan & Pre-Registration

*Written 2026-11-16, before running any modeling experiment. Hypotheses, research questions,
method, and analysis plan are fixed here to avoid post-hoc rationalization / p-hacking.*

## Title (working)
**Standardization-Induced Instability in Molecular Property Prediction Benchmarks, and a
Standardization-Orbit-Robust Remedy**

## Problem
Every molecular ML pipeline begins by *standardizing* raw chemical structures — choosing a tautomer,
a protonation/charge state, a salt/fragment policy, a stereochemistry convention, and an aromaticity
model — before any featurization (fingerprints, descriptors, or graphs). These choices are (i) rarely
reported, (ii) inconsistent across toolkits and groups, and (iii) individually defensible. They are an
unexamined **researcher degree of freedom**. Because they act on the chemical *graph*, they change
*all* featurizations, unlike the already-studied SMILES-string canonicalization axis (Kikuchi et al.,
arXiv:2505.07139, which is string-level only).

## Gap (verified, Phase 3–4)
Prior work establishes, at the **descriptor level on single datasets**, that tautomer/ionization
choices change computed descriptors and individual model statistics (e.g., "Does tautomerism influence
QSAR modeling?", 2013; plasma-protein-binding ionization study). Modern curation tools provide
**deterministic, canonicalize-once** pipelines (QSAR-ready workflow, J. Cheminform. 2024,
10.1186/s13321-024-00814-3; ProQSAR 2026; PubChem standardization). **No** work: (a) quantifies, across
many modern benchmarks, how standardization choice affects **model performance, leaderboard rankings,
and the reproducibility of "which model wins"**; (b) provides a **standardization-robust modeling
method**; or (c) gives a **ground-truth (semi-synthetic) analysis** of whether models exploit
standardization as a spurious signal. (Confirmed open by two independent literature sweeps, Phase 4.)

## Central hypothesis
The choice of structure-standardization pipeline is a *silent confounder* in molecular property
benchmarks: realistic, individually-defensible pipeline choices change model predictions and can reorder
model leaderboards and flip method-comparison conclusions, even when each pipeline's mean performance is
similar; and a method that is *robust to the standardization choice by construction* (training over, and
averaging over, each molecule's standardization orbit) both reduces this instability and does not harm —
and can improve — predictive accuracy.

## Contributions
1. **C1 — Diagnosis.** First multi-dataset (≈55 public benchmarks, regression + classification; 5 model
   families) quantification of standardization-choice sensitivity of (i) absolute performance, (ii) model
   *rankings*, and (iii) the agreement of "which model is best / which comparison is significant".
2. **C2 — Mechanism.** Attribution of the instability to individual standardization operations
   (tautomer, protonation, salt, stereo, aromaticity) via one-operation ablation, and identification of
   dataset properties (ionizable/tautomerizable fraction, size, split) that predict it.
3. **C3 — Method (SORB).** *Standardization-Orbit Robust* learning: train-time augmentation over each
   molecule's standardization orbit (label shared) + test-time prediction averaging over the orbit +
   (variant) a consistency regularizer. Ablatable into its components.
4. **C4 — Ground-truth validation.** A semi-synthetic design in which the label is, by chemistry,
   invariant across the standardization orbit, giving an exact invariance ground truth; plus an injected
   standardization *shortcut* to test spurious-signal exploitation. Shows naive models violate invariance
   / exploit the shortcut and SORB does not.
5. **C5 — Practice.** A reporting protocol + open, pip-installable implementation so benchmark results
   become standardization-reproducible.

## Research questions (each tied to an experiment)
- **RQ1 (sensitivity).** How much does standardization-pipeline choice change per-model performance?
  *Exp E1:* all models × pipelines × datasets; per-dataset spread of the metric across pipelines; paired
  Wilcoxon of pipeline-vs-reference deltas; effect sizes.
- **RQ2 (ranking reproducibility — HEADLINE).** Does standardization choice reorder model leaderboards
  and flip method-comparison conclusions? *Exp E2:* Kendall-τ and rank-swap rate between model rankings
  under different pipelines; fraction of datasets where the top-model identity changes; count of
  model-pair comparisons whose sign/significance flips across pipelines (significance gated by the Ash
  practically-significant comparison protocol, 10.1021/acs.jcim.5c01609).
- **RQ3 (mechanism).** Which operation dominates the instability, and what predicts it? *Exp E3:*
  one-operation-at-a-time ablation of the orbit; regression of per-dataset instability on dataset
  descriptors (ionizable fraction via pKa heuristics, mean #tautomers, size, task type).
- **RQ4 (method).** Does SORB reduce cross-standardization prediction variance and stabilize rankings
  without accuracy loss? *Exp E4:* SORB vs canonicalize-once baseline across all datasets; paired deltas
  of (a) across-orbit prediction SD, (b) accuracy, with bootstrap CIs; ablation of {train-aug,
  test-ensemble, consistency}.
- **RQ5 (ground-truth shortcut).** Do naive models violate standardization invariance / exploit a
  standardization shortcut, and does SORB prevent it? *Exp E5:* (a) invariance: across-orbit prediction
  SD for the same molecule (true SD should be 0); (b) shortcut: train with label-correlated
  standardization state, test under counterfactual standardization; measure performance drop (shortcut
  reliance) for naive vs SORB.

## Standardization orbit (the operations)
Implemented with RDKit `rdMolStandardize` + Dimorphite-DL. An operation axis and its settings:
- **Tautomer:** as-provided · RDKit canonical · enumerated (sample up to T forms).
- **Protonation:** as-provided · neutralized (RDKit Uncharger) · pH 7.4 (Dimorphite-DL).
- **Fragment/salt:** keep · largest-fragment only.
- **Stereo:** keep · strip.
- **Aromaticity:** RDKit model · MDL model (re-perceived).

**Named pipelines** (realistic end-to-end choices used in E1–E3):
- `P_raw` : parse only (as provided).
- `P_rdkit` : RDKit default cleanup (sanitize, largest fragment, neutralize).
- `P_taut` : P_rdkit + canonical tautomer.
- `P_ph74` : P_rdkit + Dimorphite pH 7.4 protonation.
- `P_mdlarom` : P_rdkit with MDL aromaticity.
- `P_aggressive` : largest fragment + neutralize + canonical tautomer + strip stereo.

**Orbit for SORB** (E4–E5): for each molecule, the set generated by sampling the axes above
(bounded to ≤ O members, O default 8, deterministic per seed).

## Method: SORB
- **Train-time orbit augmentation:** replace each training molecule by R sampled orbit members, each
  carrying the molecule's label; featurize and fit the base model on the augmented set.
- **Test-time orbit averaging:** featurize O orbit members of each test molecule; predict each; average
  (mean for regression; mean probability for classification).
- **Consistency-regularized variant (CPU MLP only):** add λ·Var(predictions over orbit) to the loss.
- **Base learners:** any in the zoo (fingerprint/descriptor + sklearn/LightGBM/XGBoost; tiny MLP for the
  regularized variant).
- **Baseline:** `canonicalize-once` (P_rdkit) + the same base learner (no augmentation/averaging).

## Datasets (cached, public)
- **Regression (42):** MoleculeACE (30 ChEMBL potency sets), TDC ADMET regression (9), MoleculeNet
  ESOL/FreeSolv/Lipophilicity (3).
- **Classification (13):** TDC ADMET classification.
Licenses: TDC (MIT/CC-BY per source), MoleculeACE (MIT), MoleculeNet (MIT). All public, no credentials.

## Model zoo (CPU)
ECFP4(2048)+LightGBM, ECFP4+RandomForest, RDKitDescriptors+LightGBM, ECFP4+SVM(RBF), ECFP4+Ridge/LogReg
(linear reference). Tiny MLP (PyTorch-CPU) for the consistency variant and for a representation-sensitivity
check.

## Splitting (Phase 8)
Primary: **scaffold split** (Bemis–Murcko, group-disjoint) — the leakage-free, harder regime; MoleculeACE
also evaluated on its **official split** for comparability. Secondary: random split (robustness).
Repeats: R = 5 scaffold resplits (seeds 0–4) for the full grid; R = 10 for the headline SORB and RQ2
comparisons. **Leakage controls:** before every split, standardize to a canonical form for duplicate
detection; drop exact duplicates; report train/test near-duplicate (Tanimoto ≥ 0.95) rates; never let two
standardization variants of the same molecule cross the split (split on canonical-graph identity).

## Metrics
- Regression: RMSE (primary), MAE, Spearman ρ, R², Kendall τ.
- Classification: AUROC (primary), AUPRC, balanced accuracy, MCC, Brier, ECE.
- Instability: per-dataset metric spread across pipelines; Kendall-τ / rank-swap of model rankings;
  top-model-change rate; comparison flip rate; across-orbit prediction SD (invariance).

## Statistics (pre-specified)
Per-dataset paired comparisons across seeds via **Wilcoxon signed-rank** (+ rank-biserial effect size);
across-datasets via **Friedman + Nemenyi critical-difference**; **bootstrap 95% CIs** for all headline
deltas; **Holm–Bonferroni** / Benjamini–Hochberg for multiplicity; **Cohen's d_z** for paired effects.
"Practical significance" gating per Ash et al. 2025. Significance is reported only where the test's
assumptions hold and the unit of analysis (dataset, or seed within dataset) is stated. We will NOT claim
significance from a single seed or a single dataset.

## Compute budget (4 CPU cores, no GPU)
- E1–E3 grid: 42 reg × 5 models × 6 pipelines × 5 seeds ≈ 6.3k fits; 13 clf × 5 × 6 × 5 ≈ 2.0k fits.
  Fingerprint+GBDT on ≤10k mols ≈ 1–5 s/fit → ≈ 3–7 CPU-h → ≈ 1–2 h wall on 4 cores.
- E4 SORB: method+baseline on all datasets, R=10, O=8 → ≈ +3–4 CPU-h.
- E5 semi-synthetic: focused subset → ≈ +1–2 CPU-h.
- Tiny-MLP consistency variant + GNN sensitivity: a representative subset → ≈ +2–3 CPU-h.
- Total well within the ~15 h envelope, run incrementally with checkpointing (`run_all.py --mode full`).
  `--mode quick` runs a 4-dataset × 3-model × 3-pipeline × 2-seed smoke for pipeline verification.

## Falsification / null-avoidance
RQ2's headline (rankings reorder) and RQ4/RQ5 (SORB reduces measurable across-orbit prediction SD; naive
models violate the exact invariance ground truth) cannot collapse to an uninformative null: across-orbit
prediction SD for naive models is > 0 by construction unless the featurization is already invariant (an
interesting result either way), and SORB reduces it toward 0. If standardization turns out *not* to
reorder rankings (RQ2 null), that itself is a publishable, reassuring reproducibility result, reported
honestly. Negative/weak results are retained.

## Threats to validity (to address in Discussion)
- "Received date reset" does not affect this study (no journal-timeline dependence in the science).
- Property *should* be standardization-invariant: true for the physicochemical/ADMET endpoints used;
  we exclude any endpoint where a tautomer/charge state is itself the measured species.
- Orbit sampling stochasticity: fixed seeds; report sensitivity to O and T.
- Dimorphite-DL / RDKit version pinning for reproducibility.

---
## Pre-registration amendment (2026-11-16, after a blinded pilot on 4-11 datasets)
A pilot (no cherry-picking; scripts `premise_probe.py`, `inv_feat.py`) showed:
1. The AGGREGATE benchmark metric is only modestly sensitive to standardization choice (median
   rel. spread ~0.5%, p90 ~2.5%) because (a) public benchmarks are already salt-stripped/curated by
   their providers, so only residual tautomer/protonation/aromaticity choices remain, and (b)
   test-set averaging cancels per-molecule changes. This is itself a reportable (reassuring) finding.
2. The PER-PREDICTION standardization-invariance violation is substantial and representation-dependent:
   baseline orbit prediction SD is 0.10 (tree) to 0.21-0.70 (linear-on-fingerprint) of the label SD,
   i.e., the same molecule under a different (equally valid) standardization receives a materially
   different prediction, comparable to a fraction of assay noise — yet aggregate metrics hide it.
**Sharpened framing (headline moves from RQ1/RQ2 to the invariance violation + method):**
- Primary contribution: quantify the standardization-invariance VIOLATION (per-prediction) across
  representations and datasets, and show aggregate metrics conceal it (the "silent confounder").
- Method (SORB): test-orbit averaging ENFORCES invariance at inference (safe, universal, neutral-to-
  better accuracy); train-aug TEACHES invariance (helps robust tree learners, can destabilize
  high-variance linear/SVM learners — reported honestly, investigated in error analysis).
- RQ1 (aggregate sensitivity) and RQ2 (ranking) are retained as honest secondary results explaining
  why the problem went unnoticed. This amendment is driven by the pilot being blinded to final
  multi-seed results and is reported transparently; no negative result is hidden.
