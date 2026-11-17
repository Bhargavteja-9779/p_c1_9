"""Generate supplementary_filled.md with real tables from the machine-readable result files."""
import os, json
import pandas as pd
H = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(H, "..")
RES = os.path.join(ROOT, "code", "results"); TAB = os.path.join(ROOT, "code", "tables")
def csv(d, n):
    p = os.path.join(d, n); return pd.read_csv(p) if os.path.exists(p) else None
def md(df, **kw): return df.to_markdown(index=False, **kw)
m = json.load(open(os.path.join(RES, "metrics.json")))
DROP = {"molnet:lipophilicity"}
out = ["# Supporting Information\n",
       "## Chemically Equivalent, Computationally Different: Standardization-Invariance Violations in Molecular Property Prediction\n",
       "All tables are generated from machine-readable result files by `supplementary/build_supplementary.py`.\n"]

# S1 datasets + leakage audit
la = csv(RES, "leakage_audit.csv")
if la is not None:
    la2 = la[~la.dataset.isin(DROP)].copy()
    out.append("## Table S1. Datasets and data-integrity audit\n")
    out.append("Columns: n molecules; exact-duplicate groups by standardization-invariant canonical identity "
               "(dup_groups) and groups with disagreeing labels (dup_label_disagree); distinct raw SMILES "
               "collapsing to one molecule after standardization (std_collapse_pairs); canonical-identity overlap "
               "across the scaffold split (scaffold_id_leak); train/test near-duplicate rate at Tanimoto>=0.95.\n")
    out.append(md(la2[["dataset","task","n","dup_groups","dup_label_disagree","std_collapse_pairs","scaffold_id_leak","train_test_nn_ge95pct"]]))
    # excluded-from-method datasets
    big = la2[la2.n > 6000][["dataset","n"]]
    out.append(f"\n**Datasets excluded from the method/invariance/shortcut analyses (>6000 molecules, memory):** "
               + ", ".join(f"{r.dataset} ({r.n})" for r in big.itertuples()) + ".\n")

# S2 sensitivity (RQ1) summary by model x feat
r1 = csv(TAB, "rq1_sensitivity.csv")
if r1 is not None:
    r1 = r1[~r1.dataset.isin(DROP)]
    s = r1.groupby(["model","feat"]).rel_spread_pct.agg(["median","mean","max"]).round(3).reset_index()
    out.append("\n## Table S2. Aggregate-metric sensitivity to the standardization pipeline (RQ1)\n")
    out.append("Relative spread (%) of the primary metric across the six pipelines, per model x featurization.\n")
    out.append(md(s))

# S3 ablation (RQ3)
r3 = csv(TAB, "rq3_ablation.csv")
if r3 is not None:
    r3 = r3[~r3.dataset.isin(DROP)]
    s = r3.groupby("operation").abs_delta.agg(["median","mean"]).round(4).reset_index().sort_values("median", ascending=False)
    out.append("\n## Table S3. Operation ablation (RQ3)\n")
    out.append("Median/mean absolute change in the primary metric attributable to each single standardization operation.\n")
    out.append(md(s))

# S-method (RQ4) full table
r4 = csv(TAB, "rq4_method_accuracy.csv")
if r4 is not None:
    out.append("\n## Table S4. SORB accuracy effect, per model x featurization x task x variant (RQ4)\n")
    out.append("Improvement = positive means better than baseline (RMSE reduced / AUROC increased). "
               "holm_p_reject marks the per-cell Wilcoxon test significant after Holm-Bonferroni.\n")
    cols = [c for c in ["model","feat","task","variant","n_datasets","median_improvement","frac_datasets_improved","wilcoxon_p","holm_p_reject"] if c in r4.columns]
    out.append(md(r4[cols].round(4)))

# S-clf invariance
clf = csv(TAB, "rq5_baseline_invariance_classification.csv")
if clf is not None:
    out.append("\n## Table S5. Classification per-prediction probability SD (RQ5)\n")
    out.append("Baseline probability SD across the orbit and across the realistic pipelines (0-1 scale), by model x featurization.\n")
    out.append(md(clf.round(4)))

# S-shortcut
r5b = csv(TAB, "rq5b_shortcut.csv")
if r5b is not None:
    out.append("\n## Table S6. Semi-synthetic standardization-shortcut experiment (RQ5b)\n")
    out.append("Median drop when the shortcut is removed: baseline (exploits) vs SORB (resists).\n")
    out.append(md(r5b.round(4)))

# omnibus
om = m.get("RQ4_omnibus", {})
if om:
    out.append("\n## S7. Across-dataset omnibus and multiplicity\n")
    out.append(f"- Friedman over {{baseline, train-aug, test-orbit, full}} (regression RMSE, n={om['friedman']['n_datasets']} datasets): "
               f"p = {om['friedman']['p']:.2g}; mean ranks {({k:round(v,2) for k,v in om['friedman']['mean_ranks'].items()})}; "
               f"Nemenyi CD = {om['friedman']['nemenyi_cd']:.2f}.\n")
    out.append(f"- Multiplicity: of {om['multiplicity']['n_tests']} per-cell SORB-vs-baseline tests, "
               f"{om['multiplicity']['n_sig_raw_0.05']} are significant at 0.05 raw and "
               f"{om['multiplicity']['n_sig_holm_0.05']} after Holm-Bonferroni.\n")

out.append("\n## S8. Reproducibility\nExact versions in `code/requirements.txt`; one-command reproduction `python run_all.py --mode full`; "
           "seeds 0-2; revision experiments in `code/scripts/revision_experiments.py`. Hardware: 4 CPU cores, no GPU.\n")
open(os.path.join(H, "supplementary_filled.md"), "w").write("\n".join(out))
print("supplementary_filled.md written")
