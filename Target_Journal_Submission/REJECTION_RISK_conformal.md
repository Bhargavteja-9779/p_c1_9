# Rejection-risk report — conformal paper

Honest pre-submission assessment of the most likely reviewer objections, our current defense, and the
residual risk. "Defense" points to evidence already in the manuscript/SI; "Residual" states what we
cannot fully rule out.

---

### R1. "Marginal-vs-conditional coverage is textbook; where is the novelty?"
**Risk: HIGH (this is the exact failure mode of the previous, rejected paper — novelty).**
- **Defense.** The *theoretical* caveat is known; what is new is (i) the first benchmark-scale
  *quantification* of the failure along the applicability-domain axis for molecular regression (41
  datasets × 3 models × 2 featurizations, scaffold splits), (ii) the demonstration that the community's
  default remedy (normalized conformal) does **not** fix it, and (iii) a concrete, guaranteed,
  no-retraining method (Mondrian-AD + shrinkage) with a clear recommendation. Novelty is the
  *empirical finding + method + recommendation*, not a new theorem.
- **Residual.** A reviewer may still consider the method a direct application of Mondrian conformal
  (true) and want a sharper methodological delta. Mitigation in text: the smoothed-AD shrinkage variant
  and the AD-stratum design are the methodological contributions; we are explicit that the guarantee is
  inherited, not newly proven.

### R2. "Mondrian conformal on a chosen covariate is not new (Svensson/Norinder already use it)."
**Risk: HIGH.**
- **Defense.** Prior cheminformatics Mondrian conformal conditions on **class** (imbalanced
  classification) or on predefined categories; we condition on a *continuous applicability-domain
  descriptor* discretized into strata, targeting the novel-chemotype coverage failure specifically, and
  we show this choice of taxonomy is what matters (scaffold-Mondrian does **not** help; AD-Mondrian
  does). The contribution is the *axis* (applicability domain) and the evidence that it is the right
  one, plus the shrinkage variant for sparse novel strata.
- **Residual.** Reviewer may still see it as incremental. Framing emphasizes the finding + the negative
  controls (normalized, scaffold-Mondrian) that make the positive result meaningful.

### R3. "You crippled the normalized baseline."
**Risk: MEDIUM (mitigated).**
- **Defense.** Normalized conformal uses the standard additive-regularized score
  |r|/(σ(x)+β) with β = median(σ) (hyperparameter-free), half-width q·(σ+β) — the textbook form
  [Papadopoulos 2002; Lei et al. 2018; Angelopoulos & Bates 2023]. It attains marginal coverage (as a
  valid conformal procedure must); our claim is only that it does not repair *conditional* coverage,
  which we show empirically. No width blow-up, no unfair handicap.
- **Residual.** A reviewer might want conformalized quantile regression (CQR) as an additional
  heteroscedastic baseline. We cite CQR; a limitation note states CQR with a strong conditional base
  model is future work. (Could add if reviewer insists.)

### R4. "AD distance (1 − top-k Tanimoto) is arbitrary; results may hinge on it."
**Risk: MEDIUM.**
- **Defense.** It is a standard, widely-used AD descriptor; the strata are quintiles (no tuned
  threshold). The finding (monotone coverage decay with distance) is descriptor-consistent, and
  Mondrian-AD's guarantee holds *per stratum regardless of how strata are defined*, so the method does
  not depend on the descriptor being "correct."
- **Residual.** We do not sweep alternative AD descriptors (embedding distance, density). Explicitly
  listed as a limitation and future work.

### R5. "Random-split results are weak / the effect is a scaffold-split artifact."
**Risk: LOW–MEDIUM.**
- **Defense.** We report both splits. The AD gradient is milder under random splits (expected: less
  novelty), yet marginal conformal still under-covers the far stratum and Mondrian-AD still restores it
  (SI S2.2). This strengthens, not weakens, the claim: the failure tracks genuine novelty.
- **Residual.** None major; the framing is honest about the milder random-split effect.

### R6. "Width cost makes the method impractical."
**Risk: LOW.**
- **Defense.** The width increase is modest and *concentrated in the novel strata where it is
  warranted* (Figure 2, Section 7.2). Reporting larger honest uncertainty for extrapolation is the
  desired behavior, not a defect; familiar-compound intervals are barely affected.
- **Residual.** Some applications want narrow intervals everywhere; we state the trade-off plainly and
  offer smoothed-AD as the tighter default.

### R7. "Only CPU-scale models (fingerprints/descriptors); no deep models."
**Risk: LOW–MEDIUM.**
- **Defense.** The claim is about the *calibration layer*, which is model-agnostic; we validate across
  three model families and two featurizations. Deep/GNN base models are an extension, not a refutation —
  the conditional-coverage failure is a property of marginal calibration, not of the base learner.
- **Residual.** Listed as a limitation; a reviewer may still want one deep example.

### R8. "Exchangeability may be violated (temporal/assay shift)."
**Risk: LOW.**
- **Defense.** We state the exchangeability assumption explicitly and that the guarantee degrades
  gracefully under its violation; scaffold splits are a conservative proxy for prospective novelty.
- **Residual.** A true time-split prospective study is future work (noted).

### R9. Reproducibility / statistics.
**Risk: LOW.**
- **Defense.** Every number is generated by a token-substitution pipeline from machine-readable result
  files (no hand-typed values). Paired Wilcoxon + bootstrap CIs + Friedman/Nemenyi omnibus; multiple
  seeds; both splits; multiple α in SI. Open code, seeds, environment.
- **Residual.** None material.

---

## Overall assessment
The dominant risk is **novelty perception (R1, R2)** — the same axis that sank the prior paper. The
defense rests on the combination that none of the prior works provide: a quantified benchmark-scale
*failure*, *negative controls* that rule out the obvious fixes (normalized, scaffold-Mondrian), a
*working* guaranteed method, and an actionable recommendation. The result is a strong positive finding
(not a null/best-practice note), which is what distinguishes it from the earlier submission. Framing
throughout foregrounds the finding and the negative controls so the contribution cannot be read as "just
apply Mondrian conformal."

**Recommended pre-submission actions:** (a) keep the title's "false sense of security" framing (it names
the finding, not a method tweak); (b) lead the abstract and intro with the *quantified failure*; (c)
ensure the negative controls (normalized retains marginal but not conditional; scaffold-Mondrian fails)
are prominent; (d) consider adding CQR as a baseline only if a reviewer requests it.
