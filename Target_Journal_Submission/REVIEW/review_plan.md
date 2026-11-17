# Adversarial review plan (Phases 24-33): 10 independent reviewer perspectives

Each round adopts a distinct, hostile perspective, finds NEW issues (no repeats), classifies them
FATAL/MAJOR/MODERATE/MINOR, and specifies concrete fixes. After each round, author-mode applies all
feasible fixes (rerunning experiments where needed) before the next round. Recorded in review_01..10.md
with the fixes applied and residual risks. Final adversarial pass + audit in final_audit.md.

1. review_01 — Novelty & prior-art: is the contribution genuinely new vs Kikuchi (string axis),
   tautomer-augmentation (2021), curation workflows, assay-heterogeneity work? Is the gap real?
2. review_02 — Methodology & math: is the invariance definition sound? Is SORB correct? Are the orbit
   and tautomer-scoring choices justified? Any circularity in the invariance measurement?
3. review_03 — Experimental design & fairness: splits, leakage, hyperparameters, baselines,
   model/featurizer choices, pipeline realism, seed count.
4. review_04 — Statistics: unit of analysis, Wilcoxon/Friedman assumptions, multiplicity, effect sizes,
   CIs, over-claiming significance, 3-seed adequacy.
5. review_05 — Dataset validity & generalization: benchmark pre-standardization caveat, property
   invariance assumption, classification vs regression, external validity.
6. review_06 — Robustness & effect magnitude: is the per-prediction effect large enough to matter? Does
   it hold under random splits, orbit-size changes? Is "comparable to assay noise" defensible?
7. review_07 — Computational cost: is SORB's inference/training overhead worth it? Scalability.
8. review_08 — Hostile competitor: can a simpler method (e.g., just canonicalize consistently) achieve
   the same? Is SORB necessary? Is the framing inflated?
9. review_09 — Writing, figures, tables, references, JCIM author-guideline compliance.
10. review_10 — Associate-editor + final adversarial: would this survive peer review at JCIM? One fatal
    reason; then 10 major, 20 moderate, all minor (only real ones).
