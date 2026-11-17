# Claim–Evidence Audit (Phase 22)

Every quantitative claim in the manuscript is produced by `manuscript/fill_manuscript.py` substituting a
`{{TOKEN}}` from machine-readable results (`code/results/metrics.json`, `code/tables/*.csv`,
`code/results/leakage_audit.csv`); `manuscript/values.json` is the audit trail. No result number is
hand-typed. Mapping of major claims to evidence:

| Claim (manuscript) | Token / source | Evidence file |
|---|---|---|
| 55 benchmarks, 3 model families | N_DATASETS, N_MODELS | results/raw_results.csv |
| Aggregate metric spread median 0.51%, p90 2.42% | RQ1_MEDIAN_SPREAD, RQ1_P90_SPREAD | metrics.json RQ1_sensitivity; tables/rq1_sensitivity.csv |
| Rankings mostly stable (Kendall τ median 1.0; top-model changes on 5%) | RQ2_MEDIAN_TAU, RQ2_FRAC_TOPCHANGE | metrics.json RQ2_ranking; tables/rq2_ranking.csv |
| Tautomer canonicalization is the dominant operation | RQ3_TOP_OP | metrics.json RQ3_ablation; tables/rq3_ablation.csv |
| Per-prediction orbit SD: median 0.10×, up to 0.42× (linear-fp), 0.12× (tree) | RQ5_ORBIT_* | metrics.json RQ5_invariance; tables/rq5_baseline_invariance_by_representation.csv |
| Test-orbit averaging accuracy-neutral (−0.0001, p=0.80, 50%) | SORB_TESTORBIT_* | metrics.json RQ4_method_pooled; tables/rq4_method_accuracy.csv |
| Train-aug reduces orbit SD 36% on 91% of cases | SORB_TRAINAUG_* | metrics.json RQ4_method_pooled.trainaug_invariance |
| Baseline exploits shortcut (0.102 drop); SORB resists (−0.003, p=5e-10) | SHORTCUT_DROP, SHORTCUT_SORB_DROP | metrics.json RQ5b_shortcut; tables/rq5b_shortcut.csv |
| Canonical-identity leakage 82 (0.06%); collapse/dup stats | LEAK_ID_LEAK, LEAK_PCT, N_COLLAPSE, N_DUP_DISAGREE | results/leakage_audit.csv |
| Largest dataset 13,130 molecules | BIG_DATASET_N | results/raw_results.csv |

Qualitative claims are supported by cited literature (Phase-23 citation audit: all 22 references verified
via Crossref/arXiv; see references.bib). No claim asserts an accuracy improvement from SORB (the data do
not support one); SORB is framed as a reproducibility/invariance method throughout.
