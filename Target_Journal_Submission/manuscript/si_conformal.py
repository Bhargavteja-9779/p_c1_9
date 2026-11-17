"""Generate the Supporting Information for the conformal paper from machine-readable results.
Emits conformal_SI.md with data-derived tables (per-method summary for each split, per-dataset
worst-bin coverage, per-model and per-featurization breakdowns, and alpha-robustness if present).
No number is hand-typed. Run after analyze.py (+ optional alpha-robustness analysis).
"""
from __future__ import annotations
import os, json
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "code", "results")
TAB = os.path.join(HERE, "..", "code", "tables")
LAB = {"marginal": "Marginal (standard)", "normalized": "Normalized",
       "mondrian_scaffold": "Mondrian-scaffold", "smoothed_ad": "Smoothed-AD (ours)",
       "mondrian_ad": "Mondrian-AD (ours)"}
ORDER = ["marginal", "normalized", "mondrian_scaffold", "smoothed_ad", "mondrian_ad"]


def _md_table(df, cols, headers, fmt):
    out = "| " + " | ".join(headers) + " |\n|" + "|".join(["---"] * len(headers)) + "|\n"
    for _, r in df.iterrows():
        cells = [fmt.get(c, lambda x: str(x))(r[c]) for c in cols]
        out += "| " + " | ".join(cells) + " |\n"
    return out


def summary_table(split):
    p = os.path.join(TAB, f"conformal_summary_{split}.csv")
    if not os.path.exists(p):
        return "(summary unavailable)\n"
    s = pd.read_csv(p).set_index("method").reindex(ORDER).dropna(how="all").reset_index()
    s["name"] = s["method"].map(LAB)
    f3 = lambda x: f"{x:.3f}"
    f2 = lambda x: f"{x:.2f}"
    return _md_table(
        s, ["name", "marginal_cov", "near_cov", "far_cov", "worst_bin_cov", "cond_gap_ad", "rel_width"],
        ["Method", "Marginal cov.", "Near-AD cov.", "Far-AD cov.", "Worst-bin cov.", "AD gap", "Rel. width"],
        {"marginal_cov": f3, "near_cov": f3, "far_cov": f3, "worst_bin_cov": f3, "cond_gap_ad": f3, "rel_width": f2})


def main():
    m = json.load(open(os.path.join(RES, "conformal_metrics.json")))
    d = pd.read_csv(os.path.join(RES, "conformal_results.csv"))
    d = d[d.alpha == 0.1]
    doc = []
    doc.append("# Supporting Information\n")
    doc.append("## A False Sense of Security: Standard Conformal Prediction Under-Covers Novel "
               "Chemotypes, and an Applicability-Domain-Conditional Fix\n")
    doc.append("All tables are generated directly from the machine-readable result files; no value is "
               "hand-transcribed. Coverage is reported as a fraction at the 90% target (nominal 0.90). "
               "\"Worst-bin\" is the lowest coverage across the five applicability-domain strata; "
               "\"AD gap\" is the mean absolute deviation of per-stratum coverage from nominal; "
               "\"Rel. width\" is the mean interval width divided by the label standard deviation.\n")

    # S1 methods recap
    doc.append("## S1. Method details\n")
    doc.append(
        "**Applicability-domain descriptor.** For each molecule we compute the applicability-domain "
        "(AD) distance as one minus the mean Tanimoto similarity (ECFP4, 2048 bits) to its five nearest "
        "training neighbours; larger distance = more novel. Test and calibration molecules are assigned "
        "to five AD strata by quintiles of this distance estimated on the calibration set.\n\n"
        "**Marginal** split-conformal uses the absolute residual as the nonconformity score and a single "
        "calibration quantile. **Normalized** conformal uses the textbook additive-regularized score "
        "|r|/(sigma(x)+beta) with half-width q*(sigma(x)+beta), where sigma(x) is a gradient-boosted "
        "model of the absolute residual fit on OUT-OF-FOLD training residuals (two-fold cross-fit), "
        "disjoint from the calibration and test sets so the procedure stays valid, and beta is the median "
        "calibration difficulty. **Mondrian-scaffold** calibrates a "
        "separate quantile per Bemis-Murcko scaffold cluster. **Mondrian-AD** (ours) calibrates a "
        "separate quantile per AD stratum. **Smoothed-AD** (ours) shrinks each stratum quantile toward "
        "the pooled marginal quantile, protecting sparse strata. All methods share the same base models, "
        "features, splits and calibration folds; only the nonconformity calibration differs, so "
        "differences are attributable to the calibration scheme alone.\n")

    # S2 per-method summary, both splits
    for split in ("scaffold", "random"):
        if split in [s.split("_")[-1] for s in m if s.startswith("summary_")] or os.path.exists(os.path.join(TAB, f"conformal_summary_{split}.csv")):
            doc.append(f"## S2.{'1' if split=='scaffold' else '2'} Per-method summary — {split} split "
                       f"(mean over datasets)\n")
            doc.append(summary_table(split) + "\n")

    # S3 per-dataset worst-bin coverage (scaffold): marginal vs mondrian_ad
    doc.append("## S3. Per-dataset worst-bin coverage (scaffold split)\n")
    sc = d[d.split == "scaffold"]
    piv = sc.groupby(["dataset", "method"]).worst_bin_cov.mean().unstack()
    keep = [c for c in ["marginal", "normalized", "mondrian_ad", "smoothed_ad"] if c in piv.columns]
    piv = piv[keep].round(3).reset_index()
    labels = ["Dataset"] + [LAB[c] for c in keep]
    hdr = "| " + " | ".join(labels) + " |\n|" + "|".join(["---"] * len(labels)) + "|\n"
    for _, r in piv.iterrows():
        hdr += "| " + " | ".join([str(r["dataset"])] + [f"{r[c]:.3f}" for c in keep]) + " |\n"
    doc.append(hdr + "\n")

    # S4 robustness by model and featurization (scaffold): worst-bin cov marginal vs mondrian_ad
    doc.append("## S4. Robustness by model family and featurization (scaffold split)\n")
    doc.append("Worst-bin coverage (mean over datasets/seeds), marginal vs Mondrian-AD, for each base "
               "model and featurization. The under-coverage of marginal conformal and its restoration by "
               "Mondrian-AD hold across every configuration.\n\n")
    rows = []
    for (mdl, ft), g in sc.groupby(["model", "feat"]):
        gg = g.groupby(["dataset", "method"]).worst_bin_cov.mean().unstack()
        if {"marginal", "mondrian_ad"} <= set(gg.columns):
            rows.append({"Model": mdl, "Feat": ft,
                         "Marginal worst-bin": f"{gg['marginal'].mean():.3f}",
                         "Mondrian-AD worst-bin": f"{gg['mondrian_ad'].mean():.3f}",
                         "Gain (pts)": f"{100*(gg['mondrian_ad'].mean()-gg['marginal'].mean()):.1f}"})
    rdf = pd.DataFrame(rows)
    h = "| " + " | ".join(rdf.columns) + " |\n|" + "|".join(["---"] * len(rdf.columns)) + "|\n"
    for _, r in rdf.iterrows():
        h += "| " + " | ".join(str(r[c]) for c in rdf.columns) + " |\n"
    doc.append(h + "\n")

    # S5 alpha robustness (if the sweep file exists)
    arp = os.path.join(RES, "conformal_alpha_results.csv")
    if os.path.exists(arp):
        doc.append("## S5. Robustness across significance levels\n")
        doc.append("Worst-bin coverage against the nominal target for additional confidence levels "
                   "(scaffold split, all datasets and base models, seed 0). Marginal conformal "
                   "under-covers the worst AD stratum at every level; Mondrian-AD restores it.\n\n")
        a = pd.read_csv(arp)
        rows = []
        for al, g in a.groupby("alpha"):
            gg = g.groupby(["dataset", "method"]).worst_bin_cov.mean().unstack()
            if {"marginal", "mondrian_ad"} <= set(gg.columns):
                rows.append({"Target (1-alpha)": f"{1-al:.2f}",
                             "Marginal worst-bin": f"{gg['marginal'].mean():.3f}",
                             "Mondrian-AD worst-bin": f"{gg['mondrian_ad'].mean():.3f}"})
        adf = pd.DataFrame(rows)
        h = "| " + " | ".join(adf.columns) + " |\n|" + "|".join(["---"] * len(adf.columns)) + "|\n"
        for _, r in adf.iterrows():
            h += "| " + " | ".join(str(r[c]) for c in adf.columns) + " |\n"
        doc.append(h + "\n")

    # S5b strata-count sensitivity
    srp = os.path.join(RES, "conformal_strata_results.csv")
    if os.path.exists(srp):
        doc.append("## S5b. Sensitivity to the number of applicability-domain strata\n")
        doc.append("Mondrian-AD and smoothed-AD calibrated with m ∈ {3,5,10} strata, evaluated on a fixed "
                   "5-quintile reference grid (lgbm + ECFP, scaffold split, seed 0, all datasets). Marginal "
                   "conformal is shown for reference. The restoration of worst-stratum coverage is robust to "
                   "m; more strata give finer conditioning at the cost of noisier sparse strata.\n\n")
        s = pd.read_csv(srp)
        rows = []
        mref = s[s.method == "marginal"]
        if len(mref):
            rows.append({"Method": "Marginal (ref)", "m": "—",
                         "Worst-bin cov.": f"{mref.worst_bin_cov.mean():.3f}",
                         "Far-AD cov.": f"{mref.far_cov.mean():.3f}"})
        for mth in ["mondrian_ad", "smoothed_ad"]:
            for M in [3, 5, 10]:
                g = s[(s.method == mth) & (s.m == M)]
                if len(g):
                    rows.append({"Method": LAB[mth], "m": str(M),
                                 "Worst-bin cov.": f"{g.worst_bin_cov.mean():.3f}",
                                 "Far-AD cov.": f"{g.far_cov.mean():.3f}"})
        sdf = pd.DataFrame(rows)
        h = "| " + " | ".join(sdf.columns) + " |\n|" + "|".join(["---"] * len(sdf.columns)) + "|\n"
        for _, r in sdf.iterrows():
            h += "| " + " | ".join(str(r[c]) for c in sdf.columns) + " |\n"
        doc.append(h + "\n")

    # S6 omnibus
    fr = m.get("friedman_worstcov_scaffold", {})
    if fr:
        doc.append("## S6. Omnibus test across methods (scaffold split)\n")
        doc.append(f"Friedman test over the five conformal strategies on closeness of worst-bin coverage "
                   f"to nominal (n = {fr['n_datasets']} datasets): p = {fr['p']:.2g}. Mean ranks "
                   f"(lower = closer to nominal): " +
                   ", ".join(f"{LAB[k]} {v:.2f}" for k, v in sorted(fr['mean_ranks'].items(), key=lambda x: x[1])) +
                   f". Nemenyi critical difference (alpha=0.05): {fr['nemenyi_cd']:.2f}.\n")

    open(os.path.join(HERE, "conformal_SI.md"), "w").write("\n".join(doc))
    print(f"SI written: {os.path.join(HERE, 'conformal_SI.md')}")


if __name__ == "__main__":
    main()
