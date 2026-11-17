"""Substitute {{TOKEN}} placeholders in conformal_manuscript.md with values computed from the
machine-readable conformal results (code/results/conformal_metrics.json + tables/conformal_summary_*.csv).
Every quantitative claim must trace to an experimental output; no number is hand-typed.
Run after run_conformal.py (full) + src/conformal/analyze.py.
Writes conformal_manuscript_filled.md and conformal_values.json (audit trail).
"""
from __future__ import annotations
import os, json, re

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
RES = os.path.join(ROOT, "code", "results")
SPLIT = "scaffold"   # primary split for the headline numbers; random is reported as robustness in SI
TARGET = 0.90


def _cov(x):
    """coverage as 3-decimal fraction, e.g. 0.866"""
    return f"{x:.3f}"


def _pts(x):
    """points (percentage points), 1 decimal, e.g. 8.8"""
    return f"{100*x:.1f}"


def compute_values():
    m = json.load(open(os.path.join(RES, "conformal_metrics.json")))
    summ = {r["method"]: r for r in m[f"summary_{SPLIT}"]}
    gain = m[f"worst_cov_gain_{SPLIT}"]
    fried = m.get(f"friedman_worstcov_{SPLIT}", {})
    marg, mond = summ["marginal"], summ["mondrian_ad"]
    norm = summ.get("normalized", {})
    scaf = summ.get("mondrian_scaffold", {})
    smooth = summ.get("smoothed_ad", {})
    v = {}
    v["N_DATASETS"] = str(gain["n_datasets"])
    v["N_BINS"] = "5"
    # marginal conformal
    v["MARG_COV"] = _cov(marg["marginal_cov"])
    v["MARG_NEAR_COV"] = _cov(marg["near_cov"])
    v["MARG_FAR_COV"] = _cov(marg["far_cov"])
    v["MARG_FAR_PCT"] = f"{100*marg['far_cov']:.0f}"
    v["MARG_WORST_COV"] = _cov(marg["worst_bin_cov"])
    v["MARG_GAP"] = _pts(TARGET - marg["worst_bin_cov"])          # max deviation from nominal, points
    v["MARG_RELWIDTH"] = f"{marg['rel_width']:.2f}"
    # mondrian-AD
    v["MOND_COV"] = _cov(mond["marginal_cov"])
    v["MOND_FAR_COV"] = _cov(mond["far_cov"])
    v["MOND_WORST_COV"] = _cov(mond["worst_bin_cov"])
    v["FAR_GAIN"] = _pts(mond["far_cov"] - marg["far_cov"])
    v["MOND_GAP"] = _pts(TARGET - mond["worst_bin_cov"])
    v["MOND_RELWIDTH"] = f"{mond['rel_width']:.2f}"
    # normalized / scaffold / smoothed
    if norm:
        v["NORM_COV"] = _cov(norm["marginal_cov"])
        v["NORM_WORST_COV"] = _cov(norm["worst_bin_cov"])
    if scaf:
        v["SCAF_WORST_COV"] = _cov(scaf["worst_bin_cov"])
    if smooth:
        v["SMOOTH_WORST_COV"] = _cov(smooth["worst_bin_cov"])
    # width cost of mondrian-AD over marginal (%)
    v["WIDTH_COST"] = f"{100*(mond['rel_width']/marg['rel_width'] - 1):.0f}"
    # headline gain (paired Wilcoxon, worst-bin coverage)
    v["WORST_GAIN"] = _pts(gain["median_gain"])
    p = gain["wilcoxon_p"]
    v["WORST_P"] = f"{p:.1g}" if p >= 1e-4 else f"{p:.0e}"
    v["WORST_CI_LO"] = _pts(gain["ci_lo"])
    v["WORST_CI_HI"] = _pts(gain["ci_hi"])
    v["WORST_FRAC"] = f"{100*gain['frac_improved']:.0f}"
    # Friedman omnibus over methods
    if fried:
        fp = fried["p"]
        v["FRIEDMAN_P"] = f"{fp:.1g}" if fp >= 1e-4 else f"{fp:.0e}"
        v["NEMENYI_CD"] = f"{fried['nemenyi_cd']:.2f}"
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
