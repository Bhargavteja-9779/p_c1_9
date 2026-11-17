"""Analyze the AD-conditional conformal study -> metrics.json, tables, per-dataset summaries.

Primary endpoint: FAR-stratum coverage (far_cov) — a single pre-specified bin (the most novel
applicability-domain stratum), which is unbiased (unlike worst-bin = min over 5 bins, a downward-biased
order statistic) and is the scientifically motivated 'novel chemotype' quantity. Worst-bin coverage and
the AD gap (max stratum-wise deviation) are reported as secondary, with their selection bias flagged.
We also aggregate SCAFFOLD-conditional coverage (an axis NOT used to calibrate Mondrian-AD) to test
whether the fix generalizes beyond its own calibration partition.
"""
from __future__ import annotations
import os, sys, json
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ))
import numpy as np, pandas as pd
from stat_tests.compare import paired_wilcoxon, bootstrap_ci, friedman, nemenyi_cd, holm_bonferroni
ROOT = os.path.join(os.path.dirname(__file__), "..", "..")
RES = os.path.join(ROOT, "results"); TAB = os.path.join(ROOT, "tables"); os.makedirs(TAB, exist_ok=True)
METHODS = ["marginal", "normalized", "knn_normalized", "mondrian_scaffold", "smoothed_ad", "mondrian_ad"]
TARGET = 0.9


def load():
    return pd.read_csv(os.path.join(RES, "conformal_results.csv"))


def _paired_gain(agg, col, split, name, a="mondrian_ad", b="marginal", higher_better=True):
    """Paired comparison of `col` between method a and method b over datasets."""
    piv = agg[col].unstack(0)
    if not ({a, b} <= set(piv.columns)):
        return None
    pair = piv[[a, b]].dropna()
    if len(pair) < 3:
        return None
    diff = (pair[a] - pair[b]).values
    w = paired_wilcoxon(pair[a].values, pair[b].values)
    ci = bootstrap_ci(diff, np.mean, n_boot=5000)
    return {"metric": col, "a": a, "b": b, "n_datasets": len(pair),
            "a_mean": float(pair[a].mean()), "b_mean": float(pair[b].mean()),
            "a_median": float(pair[a].median()), "b_median": float(pair[b].median()),
            "mean_gain": float(np.mean(diff)), "median_gain": float(np.median(diff)),
            "frac_improved": float((diff > 0).mean() if higher_better else (diff < 0).mean()),
            "wilcoxon_p": w["p"], "ci_lo": ci[1], "ci_hi": ci[2],
            "mean_ci": bootstrap_ci(diff, np.mean, n_boot=5000)[1:]}


def main():
    d = load(); m = {}
    d = d[d.alpha == 0.1]
    for split in sorted(d.split.unique()):
        ds = d[d.split == split]
        cols = ["marginal_cov", "far_cov", "near_cov", "worst_bin_cov", "cond_gap_ad",
                "cond_gap_scaffold", "rel_width", "far_width"]
        agg = ds.groupby(["method", "dataset"])[cols].mean()
        summ = agg.groupby("method").mean().round(4)
        summ.to_csv(os.path.join(TAB, f"conformal_summary_{split}.csv"))
        m[f"summary_{split}"] = summ.reset_index().to_dict("records")
        # also store medians over datasets (for honest median reporting)
        med = agg.groupby("method").median().round(4)
        m[f"summary_median_{split}"] = med.reset_index().to_dict("records")
        # PRIMARY endpoint: far-stratum coverage, Mondrian-AD vs marginal
        fg = _paired_gain(agg, "far_cov", split, "far")
        if fg: m[f"far_cov_gain_{split}"] = fg
        # also far_cov for smoothed-AD and knn_normalized vs marginal (to credit AD-aware baselines)
        for a in ("smoothed_ad", "knn_normalized"):
            g = _paired_gain(agg, "far_cov", split, "far", a=a)
            if g: m[f"far_cov_gain_{a}_{split}"] = g
        # SECONDARY: worst-bin coverage (biased min-statistic), Mondrian-AD vs marginal
        wg = _paired_gain(agg, "worst_bin_cov", split, "worst")
        if wg: m[f"worst_cov_gain_{split}"] = wg
        # INDEPENDENT-PARTITION check: scaffold-conditional coverage gap (axis NOT used for AD calibration)
        # lower cond_gap_scaffold is better -> compare marginal vs mondrian_ad on the gap
        sg = _paired_gain(agg, "cond_gap_scaffold", split, "scafgap", higher_better=False)
        if sg: m[f"scaffold_gap_{split}"] = sg
        # Friedman over methods on far_cov (closeness to target; higher=better = closer)
        full = agg["far_cov"].unstack(0).dropna()
        common = [c for c in METHODS if c in full.columns]
        if len(full) >= 3 and len(common) >= 3:
            score = -np.abs(full[common].values - TARGET)
            fr = friedman(score)
            m[f"friedman_farcov_{split}"] = {"p": fr["p"], "mean_ranks": dict(zip(common, fr["mean_ranks"])),
                                             "nemenyi_cd": nemenyi_cd(len(common), len(full)), "n_datasets": len(full),
                                             "n_methods": len(common)}
        # also Friedman on worst_bin (secondary)
        fullw = agg["worst_bin_cov"].unstack(0).dropna()
        if len(fullw) >= 3 and len(common) >= 3:
            fr = friedman(-np.abs(fullw[common].values - TARGET))
            m[f"friedman_worstcov_{split}"] = {"p": fr["p"], "mean_ranks": dict(zip(common, fr["mean_ranks"])),
                                               "nemenyi_cd": nemenyi_cd(len(common), len(fullw)), "n_datasets": len(fullw),
                                               "n_methods": len(common)}
    json.dump(m, open(os.path.join(RES, "conformal_metrics.json"), "w"), indent=1, default=float)
    print(json.dumps({k: v for k, v in m.items()
                      if any(t in k for t in ("far_cov_gain_scaffold", "worst_cov_gain_scaffold",
                                              "scaffold_gap_scaffold", "friedman_farcov_scaffold"))},
                     indent=1, default=float)[:1800])


if __name__ == "__main__":
    main()
