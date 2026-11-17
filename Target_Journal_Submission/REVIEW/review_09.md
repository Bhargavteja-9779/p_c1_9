# Review 09: Associate Editor — writing, structure, figures/tables, references, JCIM scope & author-guideline compliance

**Recommendation:** major_revision

**Summary:** The core study is in-scope for JCIM, clearly motivated, and the honest-reporting framing (invariance violation hidden by aggregate metrics; SORB as a reliability not accuracy method) is a genuine strength. However, as an editorial package the submission is incomplete and would not clear screening: the Supporting Information is an unpopulated template (literal {{N_DATASETS}} tokens, no actual Tables S1-S4 or Figures S1-S4) even though the main text leans on it for the leakage audit, full results, ablation, and robustness; and the main-text table apparatus is broken (a referenced "Table 3" is absent, the only embedded table is never cited, Table 1 is missing). Several claims also over-reach or contradict the paper's own numbers (abstract "benefits tree-based learners," "smallest for GBT on fingerprints," an undisclosed 49/55-dataset restriction, a model family with no reported results, inconsistent seed counts). None of the fixes require new experiments except, optionally, reporting the already-named SVM model.

## Findings

### [FATAL] Supporting Information / editorial screening
- **Issue:** The Supporting Information (Target_Journal_Submission/supplementary/supplementary.md) is a 39-line unpopulated skeleton: it contains a literal unresolved placeholder ({{N_DATASETS}}) and only prose descriptions of what Tables S1-S4 and Figures S1-S4 'will' contain — no actual tables, numbers, or figures. Yet the main text depends on the SI for load-bearing content: the data-integrity/leakage audit (Table S1, Sections 5 & 6.1), full sensitivity results (Table S2, Section 6.2), the operation ablation (Table S3, Section 6.5), the orbit-size sensitivity (Table S4, Section 8), and Figures S1-S4 (Section 8 robustness, incl. the Friedman-Nemenyi analysis cited as supporting the method conclusions). A submission whose SI is an empty template would be returned at JCIM editorial screening without review.
- **Evidence:** supplementary.md line 6 reads 'Table S1 lists all {{N_DATASETS}} datasets' with the template token never filled; the file has zero embedded result tables or figures. Main text makes 6 Table S# references and 4 Figure S# references that resolve to nothing.
- **Fix:** Populate the SI from the existing machine-readable result files (tables/rq1_sensitivity.csv, rq3_ablation.csv, table1_datasets.*, rq4/rq5 CSVs, and the Friedman-Nemenyi/random-split outputs), fill {{N_DATASETS}}, and embed the actual Tables S1-S4 and Figures S1-S4 before submission.
- **Requires new experiment:** False

### [MAJOR] Main-text tables / numbering
- **Issue:** The main-body table apparatus is incoherent. Section 7.1 and 7.2 both reference 'Table 3' (the SORB results table) as a main table, but no Table 3 is embedded in the manuscript body — only table3_method.* exist as external files. The only table actually in the body (the invariance table captioned with {#tbl:inv}) is never referenced by number in the prose (no 'Table 1' or 'Table 2' appears anywhere), so under pandoc auto-numbering it becomes 'Table 1' while the cited 'Table 3' is a dangling reference. The datasets table (table1_datasets.*) is likewise neither embedded nor cited.
- **Evidence:** grep of manuscript_filled.md: table references are only 'Table S1/S2/S3/S4' and 'Table 3' (lines 294, 317); 'Table 1', 'Table 2', and any cite of {#tbl:inv} return NONE. Only one inline table exists (line 267).
- **Fix:** Embed and number all main tables consistently: Table 1 = datasets (cite in Section 5), Table 2 = invariance-by-representation (cite tbl:inv in Section 6.4), Table 3 = SORB accuracy results (embed a condensed version from table3_method and cite in Section 7.1/7.2).
- **Requires new experiment:** False

### [MAJOR] Transparency / undisclosed sample restriction
- **Issue:** All method, invariance, and shortcut results (RQ4, RQ5, RQ5b) were run on a 49-dataset subset (40 regression + 9 classification; molecules <=6000), but this restriction is disclosed nowhere in the manuscript — 'subset', '6000', and the 49-count never appear. Worse, Section 8 claims the invariance violation and SORB's effect are consistent '(iv) across the 55 datasets', which is inaccurate since SORB never ran on all 55 (13 clf datasets exist but only 9 were used; 42 reg but only 40).
- **Evidence:** metrics.json: reg SORB n=160 (=40 datasets x 4 model-feat), clf n=36 (=9 x 4); RQ5_invariance n_datasets=49; RQ5b n_datasets=49. Task spec confirms method/shortcut experiments ran on 'the 49 datasets <=6000 molecules.' Manuscript Section 8 line 336 says 'across the 55 datasets.'
- **Fix:** State the <=6000-molecule filter in Section 5, list which datasets were excluded and why, and correct Section 8 to say the SORB/invariance results cover the 49-dataset subset (aggregate RQ1 covers all 55).
- **Requires new experiment:** False

### [MAJOR] Model zoo / unreported experiment
- **Issue:** Section 5 lists 'four model families: gradient-boosted trees (LightGBM), random forests, an RBF-kernel SVM via Nystrom approximation, and a linear model', and Section 11 Limitations also invokes 'kernel approximations' as part of the evaluated zoo. But no SVM/Nystrom result appears anywhere, and RQ1 n_cells=330 = 55 datasets x 3 models x 2 featurizations confirms only lgbm/rf/linear were run. Meanwhile the abstract, introduction, and conclusion all say '3 model families'. The manuscript claims a model family it never reports.
- **Evidence:** rq1_sensitivity.csv unique models = {lgbm, rf, linear}; metrics.json RQ1 n_cells=330 (=55x3x2); 'svm'/'Nystrom' appear only in the Section 5 model list and Section 9's 'kernel approximations'. Abstract/intro/conclusion say '3 model families' (lines 19, 71, 381).
- **Fix:** Either report the SVM/Nystrom results (if run) or delete SVM from Section 5 and the 'kernel approximations' phrase from Section 9, reconciling the zoo to the 3 models actually reported.
- **Requires new experiment:** True

### [MAJOR] Abstract / overclaiming vs. data
- **Issue:** The abstract states training augmentation 'teaches invariance and benefits robust (tree-based) learners', and Section 7.2 calls its accuracy effect 'roughly accuracy-neutral' for the gradient-boosted-tree learner. The data show train-aug significantly *worsens* LightGBM regression accuracy, not benefits/neutral: ecfp regression median improvement -0.015 RMSE, Wilcoxon p=5.8e-6, Cohen dz=-0.39, 95% CI entirely on the worse side; desc regression p=0.02, dz=-0.21. Its only non-negative effects for lgbm are on classification and are non-significant.
- **Evidence:** statistics.json RQ4_method_table: lgbm/ecfp/reg/sorb_trainaug mean -0.076, median -0.015, wilcoxon_p 5.8e-6, ci_lo -0.141, ci_hi -0.026; lgbm/desc/reg median -0.010, p=0.020. metrics.json reg_sorb_trainaug pooled median -0.013, frac_improved 0.32, p=5e-8.
- **Fix:** Reframe train-aug's benefit strictly as invariance reduction, and state plainly that it carries a small but statistically significant accuracy cost even for tree learners on regression. Change the abstract wording from 'benefits robust (tree-based) learners' accordingly.
- **Requires new experiment:** False

### [MODERATE] Results text contradicts own table
- **Issue:** Section 6.4 states the invariance violation 'is smallest for gradient-boosted trees on fingerprints (0.12x of the label SD)', but the paper's own Table 2 shows two combinations are smaller: linear+descriptors (0.071) and lgbm+descriptors (0.084), both below lgbm+ecfp (0.119). The genuinely smallest is linear+descriptors.
- **Evidence:** table2_invariance.md / tbl:inv: lgbm-desc 0.084, lgbm-ecfp 0.119, linear-desc 0.071, linear-ecfp 0.425. metrics.json RQ5_invariance matches.
- **Fix:** Correct Section 6.4 to state the violation is smallest for descriptor featurizations (linear+desc 0.071) and largest for fingerprint+linear (0.425), i.e., driven by featurization more than model class.
- **Requires new experiment:** False

### [MODERATE] Figure cross-references
- **Issue:** Section 6.2 (the aggregate-metric result) cites 'Figure 2a', but Figure 2 by definition order is the invariance-by-representation figure (fig:inv), which is a single panel with no (a)/(b). The aggregate-vs-orbit contrast that actually has panels (a)/(b) is the silent-contrast figure (fig:silent), which is the third figure (Figure 3). The manuscript hard-codes 'Figure N' numbers in prose instead of using its own {#fig:...} labels, so references are not auto-synced and at least this one points to the wrong figure.
- **Evidence:** Figure definition order in manuscript_filled.md: fig:concept(1), fig:inv(2), fig:silent(3), fig:abl(4), fig:method(5). Section 6.2 line 241 '(Figure 2a; Table S2)' describes aggregate-metric mildness, which is fig:silent panel (a) = Figure 3a.
- **Fix:** Replace hard-coded figure numbers with pandoc-crossref references ([@fig:silent] etc.) throughout, or at minimum correct 'Figure 2a' to 'Figure 3a'.
- **Requires new experiment:** False

### [MODERATE] Methods internal consistency (seeds)
- **Issue:** Section 5 (Statistics) states 'Each configuration is run with five seeds', but Section 7.1 reports the SORB experiments with 'three seeds'. The blanket 'five seeds' claim is contradicted for the headline method results.
- **Evidence:** manuscript_filled.md line 211 'Each configuration is run with five seeds'; line 294 'Table 3 and Figure 5 summarize SORB across datasets (three seeds; ...)'. Research plan specified R=5 for the grid and R=10 for headline SORB, so the realized 3-seed SORB run matches neither.
- **Fix:** Replace the blanket statement with per-experiment seed counts (e.g., sensitivity/ablation grid = 5 seeds; SORB = 3 seeds) and note the deviation from the pre-registered R=10 for SORB.
- **Requires new experiment:** False

### [MODERATE] Data and Software Availability (JCIM guideline compliance)
- **Issue:** Open code and a reporting protocol are headline contributions (C4/C5), but the availability statement gives only a placeholder repository URL ('to be inserted on acceptance'), so reviewers cannot access or verify the implementation during review, and no package versions are pinned in the statement itself (only an RDKit version note lives in the bib). JCIM expects software underpinning the claims to be available to reviewers and reproducibility details to be stated.
- **Evidence:** manuscript_filled.md lines 387-394: 'released under the MIT License at [repository URL — to be inserted on acceptance]'; no sklearn/LightGBM/Dimorphite-DL version numbers in the statement.
- **Fix:** Provide an anonymized code repository accessible to reviewers and list exact versions (RDKit 2026.09, Dimorphite-DL, scikit-learn, LightGBM) in the Data and Software Availability statement.
- **Requires new experiment:** False

### [MINOR] Unsupported adjective / overclaiming
- **Issue:** Section 7.1 calls test-time orbit averaging 'a free reliability improvement', but it is not free: Section 9 states it multiplies inference cost by the orbit size (up to 6x).
- **Evidence:** Section 7.1 line ~300 'a free reliability improvement'; Section 9 line ~346 'Test-time averaging multiplies inference cost by the orbit size (at most 6)'.
- **Fix:** Replace 'free' with 'low-cost' or 'inference-only (no retraining)', consistent with the up-to-6x inference overhead reported in Section 9.
- **Requires new experiment:** False

### [MINOR] Imprecise reporting (unit of analysis)
- **Issue:** Section 6.3 says 'on 5% of datasets the identity of the best model changes', but the ranking analysis unit is a model leaderboard per dataset-featurization (n=110 = 55 datasets x 2 featurizations), not 55 datasets. Calling the 110 leaderboards 'datasets' is imprecise.
- **Evidence:** metrics.json RQ2_ranking n_datasets=110, frac_top_model_changes=0.0545; 55 datasets x 2 featurizations = 110.
- **Fix:** State the unit explicitly, e.g., 'on 5% of dataset-featurization leaderboards the top-model identity changes'.
- **Requires new experiment:** False
