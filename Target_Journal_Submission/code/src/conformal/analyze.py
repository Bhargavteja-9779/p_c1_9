"""Analyze the AD-conditional conformal study -> metrics.json, tables, per-dataset summaries."""
from __future__ import annotations
import os, sys, json
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ))
import numpy as np, pandas as pd
from stat_tests.compare import paired_wilcoxon, bootstrap_ci, friedman, nemenyi_cd, holm_bonferroni
ROOT = os.path.join(os.path.dirname(__file__), "..", "..")
RES = os.path.join(ROOT, "results"); TAB = os.path.join(ROOT, "tables"); os.makedirs(TAB, exist_ok=True)
METHODS = ["marginal", "normalized", "mondrian_scaffold", "smoothed_ad", "mondrian_ad"]
TARGET = {0.1: 0.9}

def load():
    d = pd.read_csv(os.path.join(RES, "conformal_results.csv"))
    return d

def main():
    d = load(); m = {}
    d = d[d.alpha == 0.1]
    for split in sorted(d.split.unique()):
        ds = d[d.split == split]
        # per-method: average over models/feats/seeds to one value per dataset, then mean over datasets
        agg = ds.groupby(["method", "dataset"])[["marginal_cov", "far_cov", "near_cov", "worst_bin_cov", "cond_gap_ad", "rel_width", "far_width"]].mean()
        summ = agg.groupby("method").mean().round(4)
        summ.to_csv(os.path.join(TAB, f"conformal_summary_{split}.csv"))
        m[f"summary_{split}"] = summ.reset_index().to_dict("records")
        # headline: worst-bin coverage marginal vs mondrian_ad, paired over datasets
        piv = agg["worst_bin_cov"].unstack(0)  # index dataset, cols method
        if {"marginal", "mondrian_ad"} <= set(piv.columns):
            pair = piv[["marginal", "mondrian_ad"]].dropna()
            w = paired_wilcoxon(pair["mondrian_ad"].values, pair["marginal"].values)
            ci = bootstrap_ci((pair["mondrian_ad"] - pair["marginal"]).values, np.mean, n_boot=5000)
            m[f"worst_cov_gain_{split}"] = {"n_datasets": len(pair),
                "marginal_median": float(pair["marginal"].median()), "mondrian_median": float(pair["mondrian_ad"].median()),
                "median_gain": float((pair["mondrian_ad"] - pair["marginal"]).median()),
                "frac_improved": float((pair["mondrian_ad"] > pair["marginal"]).mean()),
                "wilcoxon_p": w["p"], "ci_lo": ci[1], "ci_hi": ci[2]}
        # Friedman over methods on worst_bin_cov (closeness to target -> use -|cov-target| so higher=better)
        full = agg["worst_bin_cov"].unstack(0).dropna()
        common = [c for c in METHODS if c in full.columns]
        if len(full) >= 3 and len(common) >= 3:
            score = -np.abs(full[common].values - 0.9)  # higher = closer to nominal
            fr = friedman(score)
            m[f"friedman_worstcov_{split}"] = {"p": fr["p"], "mean_ranks": dict(zip(common, fr["mean_ranks"])),
                                               "nemenyi_cd": nemenyi_cd(len(common), len(full)), "n_datasets": len(full)}
    json.dump(m, open(os.path.join(RES, "conformal_metrics.json"), "w"), indent=1, default=float)
    print(json.dumps({k: v for k, v in m.items() if "worst_cov_gain" in k or "friedman" in k}, indent=1, default=float)[:1500])

if __name__ == "__main__":
    main()
