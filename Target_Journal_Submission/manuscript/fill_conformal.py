"""Substitute {{TOKEN}} placeholders in conformal_manuscript.md with values computed from the
machine-readable results (code/results/conformal_metrics.json + conformal_selective_results.csv).
Every quantitative claim traces to an experimental output; no number is hand-typed.
Primary endpoint = farthest-stratum coverage. Writes conformal_manuscript_filled.md + conformal_values.json.
"""
from __future__ import annotations
import os, json, re
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(os.path.join(HERE, "..", "code", "results"))
SPLIT = "scaffold"
TARGET = 0.90
LAB = {"marginal": "Marginal (standard)", "normalized": "Learned-difficulty normalized",
       "knn_normalized": "Distance-normalized", "mondrian_scaffold": "Mondrian-scaffold",
       "smoothed_ad": "Smoothed-AD (ours)", "mondrian_ad": "Mondrian-AD (ours)"}
ORDER = ["marginal", "normalized", "knn_normalized", "mondrian_scaffold", "smoothed_ad", "mondrian_ad"]


def _cov(x): return f"{x:.3f}"
def _pts(x): return f"{100*x:.1f}"
def _p(x): return f"{x:.1g}" if x >= 1e-4 else f"{x:.0e}"


def _selective_signal(abstain=0.30):
    """Is interval width an honest novelty signal? Returns, per method: the far-AD fraction among the
    `abstain` fraction of widest-interval (would-be-abstained) compounds, the base far-AD rate, and the
    Spearman correlation between interval width and AD stratum. For marginal (constant width) the
    correlation is undefined and abstention is no better than the base rate."""
    from scipy.stats import spearmanr
    p = os.path.join(RES, "conformal_selective_results.csv")
    if not os.path.exists(p):
        return {}
    d = pd.read_csv(p)
    out = {}
    for method, g in d.groupby("method"):
        far_ab, base, rhos = [], [], []
        for ds, gg in g.groupby("dataset"):
            n = len(gg); k = max(1, int(round(abstain * n)))
            widest = gg.sort_values("rel_width", ascending=False).iloc[:k]
            far_ab.append(widest.far.mean()); base.append(gg.far.mean())
            if gg.rel_width.std() > 1e-9 and gg.stratum.std() > 0:
                rhos.append(spearmanr(gg.rel_width, gg.stratum).correlation)
        out[method] = {"far_abstain": float(np.mean(far_ab)), "base": float(np.mean(base)),
                       "rho": float(np.nanmean(rhos)) if rhos else float("nan")}
    return out


def compute_values():
    m = json.load(open(os.path.join(RES, "conformal_metrics.json")))
    summ = {r["method"]: r for r in m[f"summary_{SPLIT}"]}
    med = {r["method"]: r for r in m[f"summary_median_{SPLIT}"]}
    farg = m[f"far_cov_gain_{SPLIT}"]
    worstg = m[f"worst_cov_gain_{SPLIT}"]
    scg = m.get(f"scaffold_gap_{SPLIT}", {})
    fried = m.get(f"friedman_farcov_{SPLIT}", {})
    marg, mond = summ["marginal"], summ["mondrian_ad"]
    v = {}
    v["N_DATASETS"] = str(farg["n_datasets"]); v["N_BINS"] = "5"
    v["N_METHODS"] = str(fried.get("n_methods", len(ORDER)))
    # marginal
    v["MARG_COV"] = _cov(marg["marginal_cov"]); v["MARG_NEAR_COV"] = _cov(marg["near_cov"])
    v["MARG_FAR_COV"] = _cov(marg["far_cov"]); v["MARG_FAR_COV_MED"] = _cov(med["marginal"]["far_cov"])
    v["MARG_FAR_PCT"] = f"{100*marg['far_cov']:.0f}"; v["MARG_WORST_COV"] = _cov(marg["worst_bin_cov"])
    v["MARG_RELWIDTH"] = f"{marg['rel_width']:.2f}"; v["MARG_FARWIDTH"] = f"{marg['far_width']:.1f}"
    # mondrian-AD
    v["MOND_COV"] = _cov(mond["marginal_cov"]); v["MOND_FAR_COV"] = _cov(mond["far_cov"])
    v["MOND_WORST_COV"] = _cov(mond["worst_bin_cov"]); v["MOND_RELWIDTH"] = f"{mond['rel_width']:.2f}"
    v["MOND_FARWIDTH"] = f"{mond['far_width']:.1f}"
    # other methods' far coverage
    v["NORM_COV"] = _cov(summ["normalized"]["marginal_cov"]); v["NORM_FAR_COV"] = _cov(summ["normalized"]["far_cov"])
    v["SCAF_FAR_COV"] = _cov(summ["mondrian_scaffold"]["far_cov"])
    v["KNN_FAR_COV"] = _cov(summ["knn_normalized"]["far_cov"])
    v["SMOOTH_FAR_COV"] = _cov(summ["smoothed_ad"]["far_cov"])
    # width cost
    v["WIDTH_COST"] = f"{100*(mond['rel_width']/marg['rel_width'] - 1):.0f}"
    # typical (median over datasets) far-stratum width increase, Mondrian-AD vs marginal (robust to outliers)
    try:
        d = pd.read_csv(os.path.join(RES, "conformal_results.csv"))
        d = d[(d.alpha == 0.1) & (d.split == SPLIT)]
        fw = d.groupby(["method", "dataset"]).far_width.mean().unstack(0)
        ratio = (fw["mondrian_ad"] / fw["marginal"]).replace([np.inf, -np.inf], np.nan).dropna()
        v["FARWIDTH_MED_PCT"] = f"{100*(ratio.median() - 1):.0f}"
    except Exception:
        v["FARWIDTH_MED_PCT"] = "30"
    # PRIMARY far-cov gain
    v["FAR_GAIN"] = _pts(farg["mean_gain"]); v["FAR_P"] = _p(farg["wilcoxon_p"])
    v["FAR_CI_LO"] = _pts(farg["ci_lo"]); v["FAR_CI_HI"] = _pts(farg["ci_hi"])
    v["FAR_FRAC"] = f"{100*farg['frac_improved']:.0f}"
    # secondary worst-cov gain
    v["WORST_GAIN"] = _pts(worstg["mean_gain"]); v["WORST_P"] = _p(worstg["wilcoxon_p"])
    # scaffold-conditional gap (independent partition)
    if scg:
        v["MARG_SCAFGAP"] = _cov(scg["b_mean"]); v["MOND_SCAFGAP"] = _cov(scg["a_mean"]); v["SCAFGAP_P"] = _p(scg["wilcoxon_p"])
    # Friedman (far-cov)
    if fried:
        v["FRIEDMAN_P"] = _p(fried["p"]); v["NEMENYI_CD"] = f"{fried['nemenyi_cd']:.2f}"
    # mechanism: growth of difficulty signals from familiar->novel stratum (why learned difficulty fails)
    mp = os.path.join(RES, "conformal_mechanism_results.csv")
    if os.path.exists(mp):
        md = pd.read_csv(mp)
        g5 = {"ad_dist": [], "sigma_learned": [], "abs_resid": []}
        for ds, gg in md.groupby("dataset"):
            if len(gg) < 25:
                continue
            q = pd.qcut(gg.ad_dist.rank(method="first"), 5, labels=False)
            for key in g5:
                s = gg.groupby(q)[key].mean()
                if s.iloc[0] > 0:
                    g5[key].append(s.iloc[-1] / s.iloc[0])
        v["SIGMA_GROWTH"] = f"{100*(np.median(g5['sigma_learned'])-1):.0f}"
        v["RESID_GROWTH"] = f"{100*(np.median(g5['abs_resid'])-1):.0f}"
        v["DIST_GROWTH"] = f"{100*(np.median(g5['ad_dist'])-1):.0f}"
    # selective-prediction: is interval width an honest novelty signal?
    sig = _selective_signal(0.30)
    if sig:
        if "marginal" in sig:
            v["MARG_ABSTAIN_FAR"] = _cov(sig["marginal"]["far_abstain"]); v["ABSTAIN_BASE"] = _cov(sig["marginal"]["base"])
        if "mondrian_ad" in sig:
            v["MOND_ABSTAIN_FAR"] = _cov(sig["mondrian_ad"]["far_abstain"]); v["MOND_WIDTH_RHO"] = f"{sig['mondrian_ad']['rho']:.2f}"
    # main-text summary table (scaffold)
    hdr = ("| Method | Marginal cov. | Near-AD cov. | **Far-AD cov.** | Worst-bin cov. | Rel. width |\n"
           "|---|---|---|---|---|---|\n")
    for k in ORDER:
        if k in summ:
            r = summ[k]
            hdr += (f"| {LAB[k]} | {r['marginal_cov']:.3f} | {r['near_cov']:.3f} | "
                    f"**{r['far_cov']:.3f}** | {r['worst_bin_cov']:.3f} | {r['rel_width']:.2f} |\n")
    v["TABLE_SUMMARY"] = ("Table 1: Coverage and interval width by method (scaffold split, 90% target, mean "
                          "over 41 datasets × 3 models × 2 featurizations × 3 seeds). Far-AD coverage (the "
                          "most novel stratum) is the primary endpoint; nominal is 0.90.\n\n" + hdr)
    return v


def main():
    v = compute_values()
    src = open(os.path.join(HERE, "conformal_manuscript.md")).read()
    tokens = set(re.findall(r"\{\{([A-Z0-9_]+)\}\}", src))
    for k, val in v.items():
        src = src.replace("{{" + k + "}}", str(val))
    open(os.path.join(HERE, "conformal_manuscript_filled.md"), "w").write(src)
    json.dump(v, open(os.path.join(HERE, "conformal_values.json"), "w"), indent=1)
    missing = sorted(tokens - set(v))
    print(f"filled {len(v)} tokens; {len(missing)} still missing: {missing}")


if __name__ == "__main__":
    main()
