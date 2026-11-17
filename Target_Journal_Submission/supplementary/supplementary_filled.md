# Supporting Information

## Chemically Equivalent, Computationally Different: Standardization-Invariance Violations in Molecular Property Prediction

All tables are generated from machine-readable result files by `supplementary/build_supplementary.py`.

## Table S1. Datasets and data-integrity audit

Columns: n molecules; exact-duplicate groups by standardization-invariant canonical identity (dup_groups) and groups with disagreeing labels (dup_label_disagree); distinct raw SMILES collapsing to one molecule after standardization (std_collapse_pairs); canonical-identity overlap across the scaffold split (scaffold_id_leak); train/test near-duplicate rate at Tanimoto>=0.95.

| dataset                            | task   |     n |   dup_groups |   dup_label_disagree |   std_collapse_pairs |   scaffold_id_leak |   train_test_nn_ge95pct |
|:-----------------------------------|:-------|------:|-------------:|---------------------:|---------------------:|-------------------:|------------------------:|
| tdc:caco2_wang                     | reg    |   910 |           74 |                   39 |                   97 |                  2 |                  0.0055 |
| tdc:lipophilicity_astrazeneca      | reg    |  4200 |           82 |                   14 |                   98 |                  0 |                  0.0143 |
| tdc:solubility_aqsoldb             | reg    |  9982 |          524 |                  232 |                  848 |                 27 |                  0.003  |
| tdc:ppbr_az                        | reg    |  2790 |          691 |                  303 |                   42 |                  0 |                  0.0018 |
| tdc:vdss_lombardo                  | reg    |  1130 |           21 |                    3 |                    2 |                  0 |                  0      |
| tdc:half_life_obach                | reg    |   667 |            9 |                    0 |                    7 |                  0 |                  0      |
| tdc:clearance_hepatocyte_az        | reg    |  1213 |          201 |                  122 |                   20 |                  0 |                  0      |
| tdc:clearance_microsome_az         | reg    |  1102 |           13 |                    5 |                   14 |                  0 |                  0.0045 |
| tdc:ld50_zhu                       | reg    |  7385 |           40 |                   21 |                    3 |                  1 |                  0.005  |
| molnet:esol                        | reg    |  1128 |           13 |                    1 |                   13 |                  0 |                  0.0044 |
| molnet:freesolv                    | reg    |   642 |            3 |                    0 |                    3 |                  0 |                  0      |
| mace:CHEMBL1862_Ki                 | reg    |   794 |            0 |                    0 |                    0 |                  0 |                  0.0252 |
| mace:CHEMBL1871_Ki                 | reg    |   659 |            1 |                    0 |                    1 |                  0 |                  0.0152 |
| mace:CHEMBL2034_Ki                 | reg    |   750 |            0 |                    0 |                    0 |                  0 |                  0.0067 |
| mace:CHEMBL2047_EC50               | reg    |   631 |            0 |                    0 |                    0 |                  0 |                  0.0317 |
| mace:CHEMBL204_Ki                  | reg    |  2754 |            7 |                    0 |                    7 |                  0 |                  0.029  |
| mace:CHEMBL2147_Ki                 | reg    |  1456 |            0 |                    0 |                    0 |                  0 |                  0.0034 |
| mace:CHEMBL214_Ki                  | reg    |  3317 |            0 |                    0 |                    0 |                  0 |                  0.0588 |
| mace:CHEMBL218_EC50                | reg    |  1031 |            0 |                    0 |                    0 |                  0 |                  0.0583 |
| mace:CHEMBL219_Ki                  | reg    |  1865 |            0 |                    0 |                    0 |                  0 |                  0.0456 |
| mace:CHEMBL228_Ki                  | reg    |  1704 |            0 |                    0 |                    0 |                  0 |                  0.0293 |
| mace:CHEMBL231_Ki                  | reg    |   973 |            0 |                    0 |                    0 |                  0 |                  0.0256 |
| mace:CHEMBL233_Ki                  | reg    |  3142 |            0 |                    0 |                    0 |                  0 |                  0.0414 |
| mace:CHEMBL234_Ki                  | reg    |  3657 |            0 |                    0 |                    0 |                  0 |                  0.0438 |
| mace:CHEMBL235_EC50                | reg    |  2349 |            0 |                    0 |                    0 |                  0 |                  0.0404 |
| mace:CHEMBL236_Ki                  | reg    |  2598 |            0 |                    0 |                    0 |                  0 |                  0.0577 |
| mace:CHEMBL237_EC50                | reg    |   955 |            0 |                    0 |                    0 |                  0 |                  0.0419 |
| mace:CHEMBL237_Ki                  | reg    |  2603 |            1 |                    0 |                    1 |                  0 |                  0.0441 |
| mace:CHEMBL238_Ki                  | reg    |  1052 |            0 |                    0 |                    0 |                  0 |                  0.0238 |
| mace:CHEMBL239_EC50                | reg    |  1721 |            0 |                    0 |                    0 |                  0 |                  0.0262 |
| mace:CHEMBL244_Ki                  | reg    |  3097 |            2 |                    0 |                    2 |                  2 |                  0.0372 |
| mace:CHEMBL262_Ki                  | reg    |   856 |            0 |                    0 |                    0 |                  0 |                  0.0117 |
| mace:CHEMBL264_Ki                  | reg    |  2862 |            0 |                    0 |                    0 |                  0 |                  0.0892 |
| mace:CHEMBL2835_Ki                 | reg    |   615 |            0 |                    0 |                    0 |                  0 |                  0.0244 |
| mace:CHEMBL287_Ki                  | reg    |  1328 |            0 |                    0 |                    0 |                  0 |                  0.1165 |
| mace:CHEMBL2971_Ki                 | reg    |   976 |            0 |                    0 |                    0 |                  0 |                  0.0256 |
| mace:CHEMBL3979_EC50               | reg    |  1125 |            0 |                    0 |                    0 |                  0 |                  0.0267 |
| mace:CHEMBL4005_Ki                 | reg    |   960 |            0 |                    0 |                    0 |                  0 |                  0.0104 |
| mace:CHEMBL4203_Ki                 | reg    |   731 |            0 |                    0 |                    0 |                  0 |                  0      |
| mace:CHEMBL4616_EC50               | reg    |   682 |            0 |                    0 |                    0 |                  0 |                  0.0147 |
| mace:CHEMBL4792_Ki                 | reg    |  1471 |            1 |                    0 |                    1 |                  0 |                  0.0034 |
| tdc:ames                           | clf    |  7278 |          618 |                    4 |                  611 |                 18 |                  0.008  |
| tdc:bbb_martins                    | clf    |  2030 |           84 |                   15 |                   38 |                  3 |                  0.0049 |
| tdc:bioavailability_ma             | clf    |   640 |            4 |                    0 |                    4 |                  0 |                  0      |
| tdc:cyp2c9_veith                   | clf    | 12092 |          268 |                   22 |                  291 |                  5 |                  0      |
| tdc:cyp2d6_veith                   | clf    | 13130 |          251 |                   33 |                  273 |                  9 |                  0      |
| tdc:cyp3a4_veith                   | clf    | 12328 |          215 |                    7 |                  231 |                  5 |                  0.002  |
| tdc:cyp2c9_substrate_carbonmangels | clf    |   669 |            4 |                    0 |                    1 |                  0 |                  0      |
| tdc:cyp2d6_substrate_carbonmangels | clf    |   667 |            4 |                    1 |                    1 |                  0 |                  0      |
| tdc:cyp3a4_substrate_carbonmangels | clf    |   670 |            4 |                    1 |                    1 |                  0 |                  0      |
| tdc:dili                           | clf    |   475 |            1 |                    0 |                    1 |                  0 |                  0      |
| tdc:herg                           | clf    |   655 |           87 |                   14 |                   83 |                 10 |                  0      |
| tdc:hia_hou                        | clf    |   578 |            1 |                    0 |                    1 |                  0 |                  0.0086 |
| tdc:pgp_broccatelli                | clf    |  1218 |           41 |                    0 |                   39 |                  0 |                  0.0082 |

**Datasets excluded from the method/invariance/shortcut analyses (>6000 molecules, memory):** tdc:solubility_aqsoldb (9982), tdc:ld50_zhu (7385), tdc:ames (7278), tdc:cyp2c9_veith (12092), tdc:cyp2d6_veith (13130), tdc:cyp3a4_veith (12328).


## Table S2. Aggregate-metric sensitivity to the standardization pipeline (RQ1)

Relative spread (%) of the primary metric across the six pipelines, per model x featurization.

| model   | feat   |   median |   mean |    max |
|:--------|:-------|---------:|-------:|-------:|
| lgbm    | desc   |    0.594 |  1.174 | 15.938 |
| lgbm    | ecfp   |    0.762 |  1.277 | 13.754 |
| linear  | desc   |    0.325 |  1.803 | 57.275 |
| linear  | ecfp   |    1.008 |  1.722 |  8.646 |
| rf      | desc   |    0.377 |  0.852 | 13.593 |
| rf      | ecfp   |    0.33  |  0.607 |  4.132 |

## Table S3. Operation ablation (RQ3)

Median/mean absolute change in the primary metric attributable to each single standardization operation.

| operation    |   median |   mean |
|:-------------|---------:|-------:|
| +tautomer    |   0.0059 | 0.0524 |
| +mdlarom     |   0      | 0.0006 |
| +ph74        |   0      | 0.0042 |
| +stripstereo |   0      | 0.0019 |
| -neutralize  |   0      | 0.0355 |

## Table S4. SORB accuracy effect, per model x featurization x task x variant (RQ4)

Improvement = positive means better than baseline (RMSE reduced / AUROC increased). holm_p_reject marks the per-cell Wilcoxon test significant after Holm-Bonferroni.

| model   | feat   | task   | variant        |   n_datasets |   median_improvement |   frac_datasets_improved |   wilcoxon_p | holm_p_reject   |
|:--------|:-------|:-------|:---------------|-------------:|---------------------:|-------------------------:|-------------:|:----------------|
| lgbm    | desc   | clf    | sorb_trainaug  |            9 |               0.0031 |                   0.5556 |       0.5703 | False           |
| lgbm    | desc   | clf    | sorb_testorbit |            9 |              -0.0018 |                   0.4444 |       0.7344 | False           |
| lgbm    | desc   | clf    | sorb_full      |            9 |              -0.0031 |                   0.4444 |       0.9102 | False           |
| lgbm    | desc   | reg    | sorb_trainaug  |           39 |              -0.0088 |                   0.359  |       0.0289 | False           |
| lgbm    | desc   | reg    | sorb_testorbit |           39 |               0.0027 |                   0.641  |       0.336  | False           |
| lgbm    | desc   | reg    | sorb_full      |           39 |               0.001  |                   0.5128 |       0.8631 | False           |
| lgbm    | ecfp   | clf    | sorb_trainaug  |            9 |               0.0093 |                   0.6667 |       0.3008 | False           |
| lgbm    | ecfp   | clf    | sorb_testorbit |            9 |              -0.0011 |                   0.4444 |       0.8203 | False           |
| lgbm    | ecfp   | clf    | sorb_full      |            9 |               0.0195 |                   0.7778 |       0.1641 | False           |
| lgbm    | ecfp   | reg    | sorb_trainaug  |           39 |              -0.0162 |                   0.2308 |       0      | True            |
| lgbm    | ecfp   | reg    | sorb_testorbit |           39 |              -0.0025 |                   0.3846 |       0.3949 | False           |
| lgbm    | ecfp   | reg    | sorb_full      |           39 |              -0.0069 |                   0.359  |       0.0372 | False           |
| linear  | desc   | clf    | sorb_trainaug  |            9 |              -0.0047 |                   0.3333 |       0.0977 | False           |
| linear  | desc   | clf    | sorb_testorbit |            9 |               0.0008 |                   0.5556 |       0.5703 | False           |
| linear  | desc   | clf    | sorb_full      |            9 |               0.0013 |                   0.5556 |       0.5703 | False           |
| linear  | desc   | reg    | sorb_trainaug  |           39 |              -0.0011 |                   0.4872 |       0.5387 | False           |
| linear  | desc   | reg    | sorb_testorbit |           39 |              -0.002  |                   0.3846 |       0.2586 | False           |
| linear  | desc   | reg    | sorb_full      |           39 |               0.0007 |                   0.5128 |       0.5296 | False           |
| linear  | ecfp   | clf    | sorb_trainaug  |            9 |              -0.0203 |                   0.1111 |       0.0117 | False           |
| linear  | ecfp   | clf    | sorb_testorbit |            9 |              -0.0031 |                   0.4444 |       0.9102 | False           |
| linear  | ecfp   | clf    | sorb_full      |            9 |              -0.0192 |                   0.1111 |       0.0977 | False           |
| linear  | ecfp   | reg    | sorb_trainaug  |           39 |              -0.3029 |                   0.1795 |       0.0043 | False           |
| linear  | ecfp   | reg    | sorb_testorbit |           39 |               0.0194 |                   0.5641 |       0.3222 | False           |
| linear  | ecfp   | reg    | sorb_full      |           39 |              -0.2135 |                   0.2564 |       0.0372 | False           |

## Table S5. Classification per-prediction probability SD (RQ5)

Baseline probability SD across the orbit and across the realistic pipelines (0-1 scale), by model x featurization.

| model   | feat   |   orbit_prob_sd |   pipe_prob_sd |
|:--------|:-------|----------------:|---------------:|
| lgbm    | desc   |          0.0311 |         0.0035 |
| lgbm    | ecfp   |          0.0329 |         0.0066 |
| linear  | desc   |          0.0179 |         0.0026 |
| linear  | ecfp   |          0.0642 |         0.0115 |

## Table S6. Semi-synthetic standardization-shortcut experiment (RQ5b)

Median drop when the shortcut is removed: baseline (exploits) vs SORB (resists).

| model   | feat   | task   |   n_datasets |   median_baseline_shortcut_drop |   median_baseline_cf_vs_honest |   median_sorb_cf_vs_honest |
|:--------|:-------|:-------|-------------:|--------------------------------:|-------------------------------:|---------------------------:|
| lgbm    | ecfp   | clf    |            9 |                          0.0578 |                         0.0158 |                    -0.0093 |
| lgbm    | ecfp   | reg    |           39 |                          0.1439 |                         0.1138 |                     0.0164 |
| linear  | ecfp   | clf    |            9 |                          0.0329 |                        -0.0017 |                     0.0203 |
| linear  | ecfp   | reg    |           39 |                          0.0964 |                         0.0007 |                     0.3029 |

## S7. Across-dataset omnibus and multiplicity

- Friedman over {baseline, train-aug, test-orbit, full} (regression RMSE, n=39 datasets): p = 2.3e-09; mean ranks {'baseline': 1.95, 'sorb_trainaug': 3.56, 'sorb_testorbit': 1.87, 'sorb_full': 2.62}; Nemenyi CD = 0.75.

- Multiplicity: of 24 per-cell SORB-vs-baseline tests, 6 are significant at 0.05 raw and 1 after Holm-Bonferroni.


## S8. Reproducibility
Exact versions in `code/requirements.txt`; one-command reproduction `python run_all.py --mode full`; seeds 0-2; revision experiments in `code/scripts/revision_experiments.py`. Hardware: 4 CPU cores, no GPU.
