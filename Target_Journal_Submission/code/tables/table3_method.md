**SORB variants vs baseline: paired accuracy improvement across datasets.**

| model   | feat   | task   | variant        | metric   |   n_datasets |   mean_improvement |   median_improvement |   frac_datasets_improved |   wilcoxon_p |   cohen_dz |   ci_lo |   ci_hi |
|:--------|:-------|:-------|:---------------|:---------|-------------:|-------------------:|---------------------:|-------------------------:|-------------:|-----------:|--------:|--------:|
| lgbm    | desc   | clf    | sorb_trainaug  | auroc    |            9 |              0.005 |                0.003 |                    0.556 |        0.570 |      0.285 |  -0.006 |   0.017 |
| lgbm    | desc   | clf    | sorb_testorbit | auroc    |            9 |             -0.001 |               -0.002 |                    0.444 |        0.734 |     -0.128 |  -0.006 |   0.004 |
| lgbm    | desc   | clf    | sorb_full      | auroc    |            9 |              0.004 |               -0.003 |                    0.444 |        0.910 |      0.154 |  -0.010 |   0.018 |
| lgbm    | desc   | reg    | sorb_trainaug  | rmse     |           40 |             -0.124 |               -0.010 |                    0.350 |        0.020 |     -0.210 |  -0.332 |   0.008 |
| lgbm    | desc   | reg    | sorb_testorbit | rmse     |           40 |              0.044 |                0.003 |                    0.650 |        0.295 |      0.247 |  -0.003 |   0.103 |
| lgbm    | desc   | reg    | sorb_full      | rmse     |           40 |             -0.147 |                0.002 |                    0.525 |        0.795 |     -0.131 |  -0.528 |   0.081 |
| lgbm    | ecfp   | clf    | sorb_trainaug  | auroc    |            9 |              0.013 |                0.009 |                    0.667 |        0.301 |      0.364 |  -0.010 |   0.034 |
| lgbm    | ecfp   | clf    | sorb_testorbit | auroc    |            9 |             -0.000 |               -0.001 |                    0.444 |        0.820 |     -0.021 |  -0.006 |   0.006 |
| lgbm    | ecfp   | clf    | sorb_full      | auroc    |            9 |              0.014 |                0.019 |                    0.778 |        0.164 |      0.379 |  -0.009 |   0.035 |
| lgbm    | ecfp   | reg    | sorb_trainaug  | rmse     |           40 |             -0.076 |               -0.015 |                    0.225 |        0.000 |     -0.392 |  -0.141 |  -0.026 |
| lgbm    | ecfp   | reg    | sorb_testorbit | rmse     |           40 |              0.049 |               -0.003 |                    0.375 |        0.327 |      0.171 |  -0.008 |   0.145 |
| lgbm    | ecfp   | reg    | sorb_full      | rmse     |           40 |              0.016 |               -0.006 |                    0.375 |        0.041 |      0.060 |  -0.045 |   0.107 |
| linear  | desc   | clf    | sorb_trainaug  | auroc    |            9 |             -0.004 |               -0.005 |                    0.333 |        0.098 |     -0.791 |  -0.006 |  -0.001 |
| linear  | desc   | clf    | sorb_testorbit | auroc    |            9 |              0.002 |                0.001 |                    0.556 |        0.570 |      0.248 |  -0.002 |   0.006 |
| linear  | desc   | clf    | sorb_full      | auroc    |            9 |             -0.001 |                0.001 |                    0.556 |        0.570 |     -0.151 |  -0.006 |   0.003 |
| linear  | desc   | reg    | sorb_trainaug  | rmse     |           40 |             -0.011 |               -0.000 |                    0.500 |        0.572 |     -0.070 |  -0.064 |   0.038 |
| linear  | desc   | reg    | sorb_testorbit | rmse     |           40 |             -0.003 |               -0.002 |                    0.400 |        0.361 |     -0.056 |  -0.019 |   0.017 |
| linear  | desc   | reg    | sorb_full      | rmse     |           40 |              0.003 |                0.001 |                    0.525 |        0.420 |      0.012 |  -0.062 |   0.078 |
| linear  | ecfp   | clf    | sorb_trainaug  | auroc    |            9 |             -0.023 |               -0.020 |                    0.111 |        0.012 |     -0.995 |  -0.036 |  -0.009 |
| linear  | ecfp   | clf    | sorb_testorbit | auroc    |            9 |             -0.002 |               -0.003 |                    0.444 |        0.910 |     -0.148 |  -0.009 |   0.005 |
| linear  | ecfp   | clf    | sorb_full      | auroc    |            9 |             -0.019 |               -0.019 |                    0.111 |        0.098 |     -0.604 |  -0.037 |   0.001 |
| linear  | ecfp   | reg    | sorb_trainaug  | rmse     |           40 |             -3.289 |               -0.292 |                    0.200 |        0.009 |     -0.347 |  -6.407 |  -0.827 |
| linear  | ecfp   | reg    | sorb_testorbit | rmse     |           40 |              0.315 |                0.020 |                    0.575 |        0.295 |      0.331 |   0.071 |   0.629 |
| linear  | ecfp   | reg    | sorb_full      | rmse     |           40 |             -2.516 |               -0.188 |                    0.275 |        0.066 |     -0.318 |  -5.080 |  -0.464 |
