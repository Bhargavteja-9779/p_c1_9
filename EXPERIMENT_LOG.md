# EXPERIMENT_LOG

(created 2026-11-15; entries appended chronologically)

## 2026-11-16 pilots (pre-full-run, documented for integrity)
- premise_probe: 59% of molecules change canonical SMILES, 10% change ECFP across 6 pipelines (dataset range 13-86% / 0-31%).
- inv_feat / method pilots: baseline orbit prediction SD 0.10-0.70 of label SD (representation-dependent: linear>tree; ecfp-linear worst). SORB train-aug reduces it ~2-3x; test-orbit averaging neutral-to-better on accuracy; train-aug can hurt linear/SVM on some sets (e.g. ppbr_az linear RMSE 28.8->60.7 — honest negative, to investigate).
- Full multi-seed run launched (55 datasets x 4 models x 2 feats x 5 seeds + sensitivity/ablation grids).
