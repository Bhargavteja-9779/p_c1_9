# EXPERIMENT_LOG

(created 2026-11-15; entries appended chronologically)

## 2026-11-16 pilots (pre-full-run, documented for integrity)
- premise_probe: 59% of molecules change canonical SMILES, 10% change ECFP across 6 pipelines (dataset range 13-86% / 0-31%).
- inv_feat / method pilots: baseline orbit prediction SD 0.10-0.70 of label SD (representation-dependent: linear>tree; ecfp-linear worst). SORB train-aug reduces it ~2-3x; test-orbit averaging neutral-to-better on accuracy; train-aug can hurt linear/SVM on some sets (e.g. ppbr_az linear RMSE 28.8->60.7 — honest negative, to investigate).
- Full multi-seed run launched (55 datasets x 4 models x 2 feats x 5 seeds + sensitivity/ablation grids).

## 2026-11-17 conformal study (current primary deliverable)
- Design: 41 regression datasets x [lgbm,rf,linear] x [ecfp,desc] x seeds[0,1,2] x {scaffold,random} x alpha=0.1; 5 conformal methods per cell (marginal, normalized[additive-beta], mondrian_scaffold, smoothed_ad, mondrian_ad). AD = 1 - mean top-5 Tanimoto to train; 5 quintile strata.
- PARTIAL validation (19 datasets, scaffold, pre-fairness-fix normalized): marginal worst-bin-cov median 0.818, far-cov 0.836, marginal-cov 0.890; mondrian_ad worst-bin 0.870 (+5.4 pts, paired Wilcoxon p=1.9e-5, 94.7% datasets improved, bootstrap CI [3.8,7.4] pts); smoothed_ad worst-bin 0.854; mondrian_scaffold worst-bin 0.810 (no help); normalized (OLD) marginal-cov 0.682 + width blow-up -> FIXED to additive-beta, full re-run.
- Friedman over 5 methods on closeness-to-nominal worst-bin (19 ds): p=1.9e-12; mean ranks mondrian_ad=smoothed_ad 1.68 (best), mondrian_scaffold 3.11, marginal 3.53, normalized 5.0 (worst).
- Integrity: every manuscript number is token-substituted from conformal_metrics.json; no hand-typed values. Negative controls (normalized, scaffold-Mondrian) reported honestly. Final numbers pending full re-run completion.
