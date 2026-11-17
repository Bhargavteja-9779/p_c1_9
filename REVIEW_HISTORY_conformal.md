# REVIEW_HISTORY — conformal paper (adversarial panel, round 1)

8 adversarial reviewers dispatched (7 returned; R3/benchmark hit the 5h limit mid-run — its likely
points, CQR/localized-CP baselines + circularity, are covered by R2/R4/R6). Consensus: the **empirical
finding is real, honestly reported, and reproducible** (R8 verified every number traces exactly to
conformal_metrics.json; R5 reproduced all tests). Verdicts ranged minor→reject, converging on **major
revision**. The exposure is (a) overclaimed theoretical guarantees, (b) novelty vs uncited prior art,
(c) circular/own-axis evaluation, (d) a missing AD-aware baseline, (e) unsupported SI claims, (f) many
labeling/consistency errors.

## Disposition (what I am changing)

### Experiments / code (consolidated rerun + new analyses)
- [DONE-data] α=0.05/0.2 robustness computed (conformal_alpha_results.csv) → add SI S5. [R2,R8 BLOCKING]
- [DONE-data] Strata-count sensitivity computed (conformal_strata_results.csv) → add SI S5b. [R2,R8 BLOCKING]
- [DO] Add **k-NN distance-normalized conformal** baseline (σ = AD distance): the decisive test of
  "difficulty axis ≠ AD axis". [R4 M1 — most important missing experiment]
- [DO] **far_cov = primary endpoint** (single pre-specified unbiased bin; paired Wilcoxon p=1.8e-12,
  +9.8 pts); flag worst-bin (min over 5) and AD-gap (max over 5) as selection-biased order statistics. [R5]
- [DO] Aggregate **scaffold-conditional coverage** (cond_gap_scaffold, already computed) to show the fix
  helps on an axis it was NOT calibrated on — answers the circularity objection. [R2,R6 B2]
- [DO] **5-stratum coverage profile figure** (run_profile.py) — supports/That replaces the unsupported
  "monotone" claim; shows marginal's empirical slope vs Mondrian-AD flat. [R5,R7,R2]
- [DO] **Downstream abstention analysis** for §8 (coverage + width vs retention when abstaining by
  interval width) OR demote §8 to qualitative — §8 currently cites non-existent SI. [R2,R6,R8 BLOCKING]
- [DO-code] Exact conformal_quantile = ⌈(n+1)(1-α)⌉-th order statistic (was np.quantile 'higher',
  off-by-one conservative). [R1 #4]
- [DO-code] Deterministic scaffold hashing (hashlib, not salted hash()). [R8 #4]
- [DO-code] Normalized β from OOF training difficulties, not σ_cal. [R1 #6]

### Theory framing (honest downgrade — R1 BLOCKING)
- Downgrade "exact finite-sample per-stratum guarantee" → **approximate/asymptotic per-stratum coverage**
  (data-dependent quantile bins are not a fixed Mondrian taxonomy).
- State explicitly **smoothed-AD carries no coverage guarantee** (heuristic shrinkage); confine guarantee
  language to hard Mondrian-AD (approximate).
- Disclose **open top bin** (most-novel test points can exceed all in-bin calibrators → residual
  under-coverage) and **cal/test non-exchangeability** under scaffold splits (marginal guarantee is
  approximate; 0.901 is empirical). Frame as the realistic prospective (covariate-shift) setting.

### Novelty repositioning (R4 BLOCKING)
- Cite + differentiate **dAD (Oršolić & Šmuc, Bioinformatics 2023, btad465)** — local/dynamic AD-conditional
  conformal; our delta: fixed-strata Mondrian wrapper, 41-benchmark single-task breadth, scaffold-Mondrian
  negative control, continuous AD-stratified curves + multi-dataset stats.
- Cite **Hosni, Gillet & Marchese Robinson, ACS Omega 2026** (concurrent; AD as reliability flag vs our
  recalibration axis).
- Cite **Papadopoulos & Haralambous 2011 / Johansson kNN-normalized** (distance-normalized conformal) —
  now added as a baseline.
- Reframe finding as "first systematic benchmark-scale quantification + simple post-hoc wrapper +
  practical recommendation", not "discovery". Temper "false sense of security" / drop "distribution shift".

### Text/consistency fixes (R2,R5,R6,R7,R8)
- Fix mean/median labeling (summary values are MEANS; report median far_cov too; gain is median).
- Friedman "four"→"five"; worst methods are **marginal + Mondrian-scaffold** (normalized ranks better).
- Unify the "gap" definition (SI said "mean", code computes MAX); relabel consistently; lead with far_cov.
- "restores to nominal": far-stratum reaches/slightly exceeds nominal (over-covers 0.927); worst-bin still
  0.861 below nominal; report residual under-coverage (freesolv 0.62, esol 0.80, aqsoldb 0.81...).
- "preserves marginal coverage" → over-covers slightly (0.901→0.929).
- Random-split "milder gradient" is false (0.126≈0.128) → reword (narrower intervals, gradient persists).
- Report far-stratum width explicitly (≈ doubles, 18→40) not just the 16% average.
- Soften smoothed-AD "robust default" (worse on the sparse cases it is advertised for).
- Update README with the conformal reproduction command chain; dedup-check ESOL∩AqSolDB. [R8,R2]

### Deferred to honest Limitations (not doing now; cost/scope)
- GNN/pretrained base model; temporal/time-split; additional AD descriptors (descriptor-space/leverage);
  localized-CP (Guan) and CQR baselines; dataset non-independence (MoleculeACE 30/41) → report per-collection
  breakdown + foreground per-dataset win-rate over raw p-value.

## Round-2 plan: after revision, a lighter verification re-review (focus: did the fixes land, any new overclaim).

## Round 2 (verification panel, 3 reviewers) + honest reframe
Reviewers: consistency (clean — all R1 fixes landed, every number traces, 4 cosmetic nits fixed),
novelty/honesty, and fresh-adversarial. The latter two converged on two legitimate, serious issues that
were FIXED honestly:

1. [BLOCKING] Abstract/Contribution-#3 claimed Mondrian-AD "improves coverage on an independent scaffold
   partition it was not calibrated on" — contradicted by our own SI S4c (scaffold-conditional gap p=0.14,
   null on the primary split; significant only on the random split). REMOVED the claim entirely.

2. [MAJOR] The distance-normalized baseline (Papadopoulos 2011, prior art) is competitive-to-BETTER than
   our Mondrian-AD on our own data: far-cov 0.892 vs 0.901 (indistinguishable under Nemenyi CD 1.18),
   cond_gap 0.099 vs 0.106 (MORE uniform), rel_width 2.98 vs 3.19 (TIGHTER), 100% vs 95% win rate. Our
   "distance-normalized helps but less uniformly / Mondrian helps most" claim was FALSE. REFRAMED the
   entire paper honestly: the contribution is the benchmark + isolation of the conditioning AXIS (only
   distance-to-training conditioning helps; learned-difficulty and scaffold do not) + recommendation; the
   three AD-conditional methods (distance-normalized, smoothed-AD, Mondrian-AD) are statistically
   indistinguishable and we explicitly do NOT claim our wrapper beats the existing baseline. Title changed
   from "...and an AD-Conditional Fix" to "...and the Case for AD-Conditional Calibration". Abstract,
   contributions, Related Work, 6.3, Section 7, Section 8, Discussion, Conclusion all reframed.

Also fixed: far-width "doubles" (mean, outlier-skewed) -> report median +~30%; added distance-normalized
to the Section 8 / SI S7 width-as-signal analysis (re-ran run_profile with knn_normalized); added Lei
et al. 2018 citation; aligned REJECTION_RISK wording (no "guaranteed/working method" overclaim).

Honest standing: this is now a benchmark + axis-isolation + practical-recommendation paper, not a
novel-superior-method paper. That is the honest reading of the data. It is a legitimate JCIM contribution
but carries residual novelty risk; the alternative (claiming superiority) is precluded by the data and by
the project's no-fabrication rule.
