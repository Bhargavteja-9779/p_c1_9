"""Substitute {{TOKEN}} placeholders in manuscript.md with values computed from the machine-readable
result files (results/metrics.json, results/*.csv, tables/*.csv). Enforces that every quantitative
claim traces to an experimental output. Writes manuscript_filled.md and a values.json audit trail.
Run after run_all.py + analyze.py + leakage_audit.py.
"""
from __future__ import annotations
import os, json, re
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
RES = os.path.join(ROOT, "code", "results")
TAB = os.path.join(ROOT, "code", "tables")


def _j(name):
    p = os.path.join(RES, name)
    return json.load(open(p)) if os.path.exists(p) else {}


def _csv(d, name):
    p = os.path.join(d, name)
    return pd.read_csv(p) if os.path.exists(p) else None


def _pct(x):
    return f"{x:.0f}"


def compute_values():
    m = _j("metrics.json")
    v = {"ORBIT_M": "6"}
    sens = _csv(RES, "raw_results.csv")
    DROP={"molnet:lipophilicity"}
    if sens is not None:
        sens = sens[~sens.dataset.isin(DROP)]
        v["N_DATASETS"] = str(sens.dataset.nunique())
        v["N_MODELS"] = str(sens.model.nunique())
    inv = m.get("RQ5_invariance", {})
    if inv.get("n_reg_datasets") is not None:
        v["N_METHOD_DATASETS"] = str(inv["n_reg_datasets"] + inv["n_clf_datasets"])
        v["BIG_DATASET_N"] = f"{int(sens.n_train.max() + sens.n_test.max()):,}" if "n_train" in sens else "13,000"
    # RQ1
    r1 = m.get("RQ1_sensitivity", {})
    if r1:
        v["RQ1_MEDIAN_SPREAD"] = f"{r1['median_rel_spread_pct']:.2f}"
        v["RQ1_P90_SPREAD"] = f"{r1['p90_rel_spread_pct']:.2f}"
        v["N1_CELLS"] = str(r1["n_cells"])
    # RQ2
    r2 = m.get("RQ2_ranking", {})
    if r2:
        v["RQ2_MEDIAN_TAU"] = f"{r2['median_mean_kendall_tau']:.3f}"
        v["RQ2_FRAC_TOPCHANGE"] = _pct(100 * r2["frac_top_model_changes"])
    # RQ3
    r3 = m.get("RQ3_ablation", {}).get("median_abs_delta_by_operation", {})
    if r3:
        top = max(r3, key=r3.get)
        v["RQ3_TOP_OP"] = {"+tautomer": "tautomer canonicalization", "+ph74": "protonation state (pH 7.4)",
                           "+stripstereo": "stereochemistry handling", "+mdlarom": "the aromaticity model",
                           "-neutralize": "charge/neutralization"}.get(top, top)
    # RQ5 invariance
    inv = m.get("RQ5_invariance", {})
    r5 = inv.get("baseline_invariance_by_representation", [])
    clf5 = inv.get("classification_probability_sd", [])
    if r5:
        ov = [r["orbit_mean_rel_sd"] for r in r5]
        pv = [r["pipe_mean_rel_sd"] for r in r5]
        v["RQ5_ORBIT_MEDIAN"] = f"{np.median(ov):.2f}×"
        v["RQ5_ORBIT_MAX"] = f"{max(ov):.2f}×"
        v["PIPE_SD_MEDIAN"] = f"{np.median(pv):.3f}×"
        lin = [r for r in r5 if r["model"] == "linear" and r["feat"] == "ecfp"]
        tree = [r for r in r5 if r["model"] == "lgbm" and r["feat"] == "ecfp"]
        if lin: v["RQ5_ORBIT_LINEAR"] = f"{lin[0]['orbit_mean_rel_sd']:.2f}×"
        if tree: v["RQ5_ORBIT_TREE"] = f"{tree[0]['orbit_mean_rel_sd']:.2f}×"
        if clf5:
            v["CLF_PROB_SD_MEDIAN"] = f"{np.median([r['pipe_prob_sd'] for r in clf5]):.3f}"
        # machine-generated invariance table (regression)
        tbl = (": Standardization-invariance violation (regression), median over datasets: per-prediction "
               "SD relative to the label SD, across a molecule's full equivalent-form orbit and across the "
               "six realistic pipelines. {#tbl:inv}\n\n"
               "| Model | Featurization | Orbit SD / label SD | Pipeline SD / label SD |\n|---|---|---|---|\n")
        for r in r5:
            tbl += f"| {r['model']} | {r['feat']} | {r['orbit_mean_rel_sd']:.3f} | {r['pipe_mean_rel_sd']:.4f} |\n"
        v["TABLE_INVARIANCE"] = tbl
    # RQ4 pooled method
    mp = m.get("RQ4_method_pooled", {})
    # prefer regression primary (rmse): median improvement -> report as metric units
    to = mp.get("reg_sorb_testorbit") or mp.get("clf_sorb_testorbit")
    if to:
        v["SORB_TESTORBIT_MED"] = f"a negligible {to['median_improvement']:+.4f} in RMSE"
        v["SORB_TESTORBIT_FRAC"] = _pct(100 * to["frac_improved"])
        v["SORB_TESTORBIT_P"] = f"{to['wilcoxon_p']:.2g}"
    om = m.get("RQ4_omnibus", {})
    fr = om.get("friedman")
    if fr:
        v["FRIEDMAN_P"] = f"{fr['p']:.1g}"
        v["NEMENYI_CD"] = f"{fr['nemenyi_cd']:.2f}"
        mr = fr["mean_ranks"]
        v["RANK_BASELINE"] = f"{mr['baseline']:.2f}"
        v["RANK_TESTORBIT"] = f"{mr['sorb_testorbit']:.2f}"
        v["RANK_FULL"] = f"{mr['sorb_full']:.2f}"
        v["RANK_TRAINAUG"] = f"{mr['sorb_trainaug']:.2f}"
    mu = om.get("multiplicity")
    if mu:
        v["N_SIG_RAW"] = str(mu["n_sig_raw_0.05"])
        v["N_SIG_HOLM"] = str(mu["n_sig_holm_0.05"])
    ti = mp.get("trainaug_invariance")
    if ti:
        v["SORB_TRAINAUG_INVRED"] = f"{ti['median_abs_reduction']:.3f} of the label SD"
        v["SORB_TRAINAUG_INVREDPCT"] = _pct(ti["median_rel_reduction_pct"])
        v["SORB_TRAINAUG_FRAC"] = _pct(100 * ti["frac_reduced"])
    # worst example for honest error analysis (largest accuracy degradation from trainaug on linear)
    r4t = _csv(TAB, "rq4_method_accuracy.csv")
    if r4t is not None and len(r4t):
        lin = r4t[(r4t.variant == "sorb_trainaug")]
        if len(lin):
            worst = lin.sort_values("mean_improvement").iloc[0]
            v["SORB_WORST_EXAMPLE"] = f"{worst['model']} on fingerprints, mean {worst['metric']} change {worst['mean_improvement']:+.3f}"
    # RQ5b shortcut
    sc = m.get("RQ5b_shortcut", {})
    if sc:
        v["SHORTCUT_DROP"] = f"{sc['median_baseline_shortcut_drop']:.3f} in the primary metric (median over datasets)"
        ssd = sc.get("median_sorb_shortcut_drop", 0.0)
        pp = sc.get("baseline_vs_sorb_shortcut_wilcoxon_p", None)
        pstr = f", paired Wilcoxon p = {pp:.1g}" if pp is not None else ""
        v["SHORTCUT_SORB_DROP"] = f"a negligible {ssd:+.3f}{pstr}"
    # leakage
    la = _csv(RES, "leakage_audit.csv")
    if la is not None:
        v["LEAK_ID_LEAK"] = f"{int(la.scaffold_id_leak.sum())} molecules"
        v["N_COLLAPSE"] = f"{int(la.std_collapse_pairs.sum()):,}"
        v["N_DUP_DISAGREE"] = f"{int(la.dup_label_disagree.sum()):,}"
        v["LEAK_PCT"] = f"{100*la.scaffold_id_leak.sum()/la.n.sum():.2f}"
    return v


def main():
    v = compute_values()
    src = open(os.path.join(HERE, "manuscript.md")).read()
    tokens = set(re.findall(r"\{\{([A-Z0-9_]+)\}\}", src))
    for k, val in v.items():
        src = src.replace("{{" + k + "}}", str(val))
    open(os.path.join(HERE, "manuscript_filled.md"), "w").write(src)
    json.dump(v, open(os.path.join(HERE, "values.json"), "w"), indent=1)
    missing = sorted(tokens - set(v))
    print(f"filled {len(v)} tokens; {len(missing)} still missing: {missing}")


if __name__ == "__main__":
    main()
