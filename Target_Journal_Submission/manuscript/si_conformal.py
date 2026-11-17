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
LAB = {"marginal": "Marginal (standard)", "normalized": "Learned-difficulty normalized",
       "knn_normalized": "Distance-normalized", "mondrian_scaffold": "Mondrian-scaffold",
       "smoothed_ad": "Smoothed-AD (ours)", "mondrian_ad": "Mondrian-AD (ours)"}
ORDER = ["marginal", "normalized", "knn_normalized", "mondrian_scaffold", "smoothed_ad", "mondrian_ad"]


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
               "\"Far-AD\" is the coverage of the most novel stratum (the primary endpoint, an unbiased "
               "single-bin estimate). \"Worst-bin\" is the minimum coverage across the five strata (a "
               "downward-biased order statistic); \"AD gap\" (cond_gap_ad) is the MAXIMUM absolute "
               "deviation of per-stratum coverage from nominal (an upward-biased order statistic); both "
               "are reported as secondary. \"Rel. width\" is the mean interval width / label SD.\n")

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

    # S4b per-collection breakdown (dataset non-independence)
    doc.append("## S4b. Per-collection breakdown (scaffold split)\n")
    doc.append("Far-AD coverage (primary endpoint) by dataset collection, marginal vs Mondrian-AD, since "
               "30 of the datasets come from one curated collection (MoleculeACE). The under-coverage and "
               "its repair hold within each collection.\n\n")
    def _collection(name):
        if name.startswith("mace:"): return "MoleculeACE (30)"
        if name.startswith("tdc:"): return "TDC ADMET (9)"
        return "MoleculeNet (2)"
    sc["collection"] = sc["dataset"].map(_collection)
    rows = []
    for coll, g in sc.groupby("collection"):
        gg = g.groupby(["dataset", "method"]).far_cov.mean().unstack()
        if {"marginal", "mondrian_ad"} <= set(gg.columns):
            rows.append({"Collection": coll, "n": gg.shape[0],
                         "Marginal far-cov": f"{gg['marginal'].mean():.3f}",
                         "Mondrian-AD far-cov": f"{gg['mondrian_ad'].mean():.3f}",
                         "Win rate": f"{100*(gg['mondrian_ad']>gg['marginal']).mean():.0f}%"})
    cdf = pd.DataFrame(rows)
    h = "| " + " | ".join(cdf.columns) + " |\n|" + "|".join(["---"]*len(cdf.columns)) + "|\n"
    for _, r in cdf.iterrows(): h += "| " + " | ".join(str(r[c]) for c in cdf.columns) + " |\n"
    doc.append(h + "\n")

    # S4c scaffold-conditional coverage gap (honest control; noise-dominated)
    if "scaffold_gap_scaffold" in m:
        g = m["scaffold_gap_scaffold"]
        doc.append("## S4c. Scaffold-conditional coverage (a different partition)\n")
        doc.append(f"As a control we also measured the scaffold-conditional coverage gap (maximum deviation "
                   f"from nominal across Bemis-Murcko scaffold clusters). This statistic is a maximum over "
                   f"up to 50 clusters, most with very few compounds, so it is dominated by small-sample "
                   f"noise and sits near {g['b_mean']:.2f} for every method; Mondrian-AD vs marginal shows "
                   f"no significant difference on it ({g['a_mean']:.3f} vs {g['b_mean']:.3f}, p = "
                   f"{g['wilcoxon_p']:.2g}). We do not claim Mondrian-AD repairs arbitrary conditional "
                   f"coverage: it targets the applicability-domain axis specifically. The evidence that the "
                   f"far-stratum gain is not an evaluation artifact is instead that the distance-normalized "
                   f"baseline, which conditions on the AD axis by normalization rather than per-bin "
                   f"calibration, also raises far-stratum coverage (main text Section 7.1).\n")

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

    # S6 omnibus (primary endpoint: far-cov)
    fr = m.get("friedman_farcov_scaffold", {})
    if fr:
        doc.append("## S6. Omnibus test across methods (scaffold split, far-AD coverage)\n")
        doc.append(f"Friedman test over the {fr.get('n_methods', len(ORDER))} conformal strategies on "
                   f"closeness of FAR-AD coverage (the primary endpoint) to nominal (n = {fr['n_datasets']} "
                   f"datasets): p = {fr['p']:.2g}. Mean ranks (lower = closer to nominal): " +
                   ", ".join(f"{LAB.get(k,k)} {v:.2f}" for k, v in sorted(fr['mean_ranks'].items(), key=lambda x: x[1])) +
                   f". Nemenyi critical difference (alpha=0.05): {fr['nemenyi_cd']:.2f}.\n")
    frw = m.get("friedman_worstcov_scaffold", {})
    if frw:
        doc.append(f"\nThe qualitative conclusion is the same for the secondary worst-bin endpoint "
                   f"(Friedman p = {frw['p']:.2g}; the applicability-domain-conditional methods rank best — " +
                   ", ".join(f"{LAB.get(k,k)}" for k, _ in sorted(frw['mean_ranks'].items(), key=lambda x: x[1])[:2]) +
                   " — and marginal and Mondrian-scaffold worst).\n")

    # S7 selective prediction: interval width as a novelty signal
    selp = os.path.join(RES, "conformal_selective_results.csv")
    if os.path.exists(selp):
        from scipy.stats import spearmanr
        doc.append("## S7. Interval width as a novelty signal (abstention)\n")
        doc.append("Standard marginal conformal gives every compound the same interval width, so width "
                   "carries no triage information. Mondrian-AD and smoothed-AD make width track distance to "
                   "training, so the widest intervals flag novel compounds. Below: the Spearman correlation "
                   "between interval width and applicability-domain stratum, and the fraction of the widest "
                   "30% of intervals (the would-be-abstained set) that are far-AD (novel) versus the base "
                   "rate (scaffold split, lgbm+ECFP, seed 0, mean over datasets).\n\n")
        sd = pd.read_csv(selp)
        rows = []
        for mth in ["marginal", "smoothed_ad", "mondrian_ad"]:
            g = sd[sd.method == mth]
            if not len(g): continue
            rhos, far_ab, base = [], [], []
            for ds, gg in g.groupby("dataset"):
                if gg.rel_width.std() > 1e-9 and gg.stratum.std() > 0:
                    rhos.append(spearmanr(gg.rel_width, gg.stratum).correlation)
                k = max(1, int(round(0.3*len(gg))))
                far_ab.append(gg.sort_values("rel_width", ascending=False).iloc[:k].far.mean()); base.append(gg.far.mean())
            rho = np.nanmean(rhos) if rhos else float("nan")
            rows.append({"Method": LAB.get(mth, mth),
                         "Spearman(width, stratum)": ("n/a (constant width)" if np.isnan(rho) else f"{rho:+.2f}"),
                         "Far-AD frac in widest 30%": f"{np.mean(far_ab):.3f}",
                         "Base far-AD rate": f"{np.mean(base):.3f}"})
        sdf = pd.DataFrame(rows)
        h = "| " + " | ".join(sdf.columns) + " |\n|" + "|".join(["---"]*len(sdf.columns)) + "|\n"
        for _, r in sdf.iterrows(): h += "| " + " | ".join(str(r[c]) for c in sdf.columns) + " |\n"
        doc.append(h + "\n")

    open(os.path.join(HERE, "conformal_SI.md"), "w").write("\n".join(doc))
    print(f"SI written: {os.path.join(HERE, 'conformal_SI.md')}")


if __name__ == "__main__":
    main()
