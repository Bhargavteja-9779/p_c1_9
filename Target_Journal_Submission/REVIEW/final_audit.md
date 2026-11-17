# Final Adversarial Audit (Phase 33–34)

Conducted after the round-1 revision, as an independent hostile pass on the *revised* manuscript and a
cross-check of every number against machine-readable outputs.

## Consistency (Phase 34) — PASS
- `python -m src.visualization.analyze` + `fill_manuscript.py` regenerate all 45 numeric tokens; 0 missing.
- Spot-checks (manuscript value vs results/metrics.json): RQ1 median spread 0.51% = 0.51%; linear+ECFP
  orbit SD 0.47× = 0.47×; Friedman p 2e-09 = 2e-09; shortcut baseline drop 0.101 = 0.101. All match.
- 11/11 unit tests pass; figures/tables regenerate from result files; DOCX has 0 unresolved citations,
  0 leftover placeholders, 5 embedded figures + 1 table.

## Round-1 FATAL issues — resolved
1. "Random splits claimed but not run" → random-split replication actually run
   (revision_results.csv, orbit_sd_random median reported in §8); unrun analyses (orbit-size sweep,
   Figures S1–S4) removed, not claimed.
2. "Supporting Information unpopulated skeleton" → SI regenerated from CSVs with real Tables S1–S7.
3. "Missing canonicalize-at-inference baseline" → added (§7.3): consistent canonical-tautomer pipeline
   collapses orbit SD from 0.12× to 0.001×; this is now the paper's central practical recommendation.

## Round-1 MAJOR themes — resolved
- Misreported setup corrected everywhere (3 model families, 3 seeds, 48 method datasets, no SVM/MLP).
- Duplicate logD benchmark removed throughout; classification invariance separated from regression
  (scale bug fixed). Friedman/Nemenyi + Holm multiplicity now computed and reported (1 of 24 per-cell
  effects survives correction — stated).
- Prior art added and novelty honestly repositioned: SORB framed as the augmentation+test-time-averaging
  paradigm [Bjerrum 2017; maxsmi 2021] on the graph axis; Kikuchi 2025 credited for the confounder
  concept; Fourches 2010/2016 curation canon cited.
- Overclaims removed: SORB reframed as reliability-not-accuracy; train-aug's significant accuracy cost
  stated; the cross-pipeline result reported as a NEGATIVE result (orbit augmentation does not reduce
  cross-canonicalizer disagreement); "free" and null-acceptance language replaced by equivalence framing.
- Headline reframed to lead with the realistic pipeline SD and treat the orbit SD as a worst-case upper
  bound; endpoint-invariance caveat (logD/permeability protonation) added.

## Residual weaknesses (honestly retained)
- The practically-relevant (realistic-pipeline) effect is modest (~0.07× label SD cross-pipeline); the
  dramatic 0.47× is a worst case (linear-on-fingerprint, full orbit). Reported as such.
- SORB's method contribution is largely negative; the paper's value is the diagnosis + the simple fix +
  the shortcut caution. This is a measurement/best-practice contribution, not a high-impact new method.
- Benchmarks are pre-standardized (effect is a lower bound); no GNN/language models; 3 seeds; CPU-only.

## Verdict
No remaining FATAL or unaddressed MAJOR issue. The revised manuscript is internally consistent, honest,
and reproducible. Its contribution is modest but real and clearly useful to practitioners.
