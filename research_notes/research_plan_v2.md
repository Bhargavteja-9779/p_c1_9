# Research Plan v2 (PIVOT) — pre-registered 2026-11-17

## Why the pivot
The standardization paper (v1) was correctly judged "limited novelty" (modest effect; trivial fix).
Six empirical validations (scratchpad) confirmed that beating strong GBDT+ECFP baselines on accuracy,
and clever UQ-efficiency tricks, do NOT yield strong positive results in this mature field. One
direction gave a clean, strong, important, validated positive result and is retained:

## Topic (LOCKED, validated)
**Trustworthy uncertainty for molecular property prediction under chemical distribution shift:
applicability-domain-conditional conformal prediction.**

### Validated core finding (8/8 datasets, scaffold split, pilot)
- Naive model error bars are badly broken: ~0.40–0.58 coverage at nominal 0.90.
- Split conformal restores MARGINAL coverage (~0.90) but gives a FALSE SENSE OF SECURITY: it
  systematically UNDER-covers NOVEL chemotypes (test molecules far from training: coverage ~0.85,
  worst ~0.79) while over-covering familiar ones (~0.93–1.00). This is exactly the extrapolation case
  that matters most in prospective drug discovery.
- **Fix (validated):** applicability-domain-conditional (Mondrian) conformal — calibrate separately by
  distance-to-training-set bins — restores novel-chemotype coverage to ~0.90 (mean 0.848 → 0.901) at a
  small width cost (~+6%), preserving marginal coverage. Similarity-localized and quantile-regression
  variants were tested and did NOT work (recorded as negative controls).

## Contributions
1. A systematic, multi-dataset quantification of applicability-domain-conditional MIScoverage of
   standard conformal prediction in molecular property regression (the "false security" finding).
2. An applicability-domain-conditional conformal method (Mondrian by distance-to-training; plus smoothed
   variants) that restores valid coverage for novel chemotypes, with finite-sample per-group validity.
3. A comprehensive benchmark across ~40 regression datasets x 3 models x multiple conformal strategies
   (marginal, normalized/heteroscedastic, Mondrian-AD, Mondrian-scaffold), measuring marginal +
   conditional coverage, worst-group coverage, width/efficiency, and prospective (novel-compound)
   selective prediction utility.
4. Practical guidance: how to obtain trustworthy uncertainty for prospective molecular predictions.

## RQs (each -> experiment)
- RQ1: How badly do naive and marginal-conformal intervals mis-cover as a function of distance-to-training
  (applicability domain) and across scaffold groups? (coverage vs AD bin; worst-group coverage)
- RQ2: Does AD-conditional (Mondrian) conformal restore conditional coverage, and at what width cost?
  (coverage by bin, width, marginal coverage preserved) — vs normalized and scaffold-Mondrian baselines.
- RQ3: Which applicability-domain descriptor and binning choice works best; sensitivity to #bins, alpha,
  calibration size, model family, featurization.
- RQ4: Does it hold under random vs scaffold vs temporal-like splits, both featurizations, GBDT/RF/linear?
- RQ5: Downstream: for prospective selection of novel compounds with a guaranteed error rate, does
  AD-conditional conformal give reliable abstention/selection where marginal conformal fails?

## Datasets, models, stats
Regression: MoleculeACE (30 ChEMBL) + TDC ADMET regression (9) + MoleculeNet ESOL/FreeSolv = ~41 (drop
duplicate logD). Models: LightGBM, RandomForest, Ridge (ECFP + descriptors). Seeds 0–4 (dataset is the
unit). Stats: Wilcoxon, Friedman/Nemenyi, bootstrap CIs, Holm. alpha=0.1 primary (also 0.2, 0.05).
Metrics: marginal coverage, per-AD-bin coverage, worst-group coverage, conditional-coverage gap (max |cov-
(1-alpha)|), mean/median interval width, size-stratified; selective-prediction error curves.

## Honesty
Mondrian conformal is an established technique; the contribution is its principled application to the
applicability-domain axis for molecular regression, the systematic "false-security" finding, the
comprehensive benchmark, and the validated practical recipe. Negative results (localized, CQR variants)
are reported. The width cost of the fix is reported honestly.
