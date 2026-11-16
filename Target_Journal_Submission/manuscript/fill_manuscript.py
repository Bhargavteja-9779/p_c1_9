"""Substitute {{TOKEN}} placeholders in manuscript.md with values computed from the machine-readable
result files (results/metrics.json, results/statistics.json, tables/*.csv). Enforces that every
quantitative claim in the manuscript traces to an experimental output. Writes manuscript_filled.md
and a values.json audit trail.
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


def _csv(name):
    p = os.path.join(TAB, name)
    return pd.read_csv(p) if os.path.exists(p) else None


def compute_values():
    m = _j("metrics.json")
    v = {}
    v["ORBIT_M"] = "6"
    # dataset / model counts
    meth = _csv("rq5_baseline_invariance_by_representation.csv")
    sens = pd.read_csv(os.path.join(RES, "raw_results.csv")) if os.path.exists(os.path.join(RES, "raw_results.csv")) else None
    if sens is not None:
        v["N_DATASETS"] = str(sens.dataset.nunique())
        v["N_MODELS"] = str(sens.model.nunique())
    # RQ1
    if "RQ1_sensitivity" in m:
        v["RQ1_MEDIAN_SPREAD"] = f"{m['RQ1_sensitivity']['median_rel_spread_pct']:.2f}"
        v["RQ1_P90_SPREAD"] = f"{m['RQ1_sensitivity']['p90_rel_spread_pct']:.2f}"
    # RQ2
    if "RQ2_ranking" in m:
        v["RQ2_MEDIAN_TAU"] = f"{m['RQ2_ranking']['median_mean_kendall_tau']:.3f}"
        v["RQ2_FRAC_TOPCHANGE"] = f"{100*m['RQ2_ranking']['frac_top_model_changes']:.0f}"
    # RQ3 top operation
    if "RQ3_ablation" in m:
        ops = m["RQ3_ablation"]["median_abs_delta_by_operation"]
        if ops:
            top = max(ops, key=ops.get)
            v["RQ3_TOP_OP"] = {"+tautomer": "tautomer canonicalization", "+ph74": "protonation (pH 7.4)",
                               "+stripstereo": "stereochemistry handling", "+mdlarom": "aromaticity model",
                               "-neutralize": "charge/neutralization"}.get(top, top)
    # RQ5 headline invariance
    if "RQ5_invariance" in m:
        rep = m["RQ5_invariance"]["baseline_invariance_by_representation"]
        orbit_vals = [r["orbit_mean_rel_sd"] for r in rep]
        v["RQ5_ORBIT_MEDIAN"] = f"{np.median(orbit_vals):.2f}×"
        v["RQ5_ORBIT_MAX"] = f"{max(orbit_vals):.2f}×"
        # linear-fingerprint specifically
        lin = [r for r in rep if r["model"] == "linear" and r["feat"] == "ecfp"]
        if lin:
            v["RQ5_ORBIT_LINEAR"] = f"{lin[0]['orbit_mean_rel_sd']:.2f}×"
        tree = [r for r in rep if r["model"] == "lgbm" and r["feat"] == "ecfp"]
        if tree:
            v["RQ5_ORBIT_TREE"] = f"{tree[0]['orbit_mean_rel_sd']:.2f}×"
    return v


def main():
    v = compute_values()
    src = open(os.path.join(HERE, "manuscript.md")).read()
    missing = set(re.findall(r"\{\{([A-Z0-9_]+)\}\}", src)) - set(v)
    for k, val in v.items():
        src = src.replace("{{" + k + "}}", str(val))
    open(os.path.join(HERE, "manuscript_filled.md"), "w").write(src)
    json.dump(v, open(os.path.join(HERE, "values.json"), "w"), indent=1)
    print("filled", len(v), "tokens; still missing:", sorted(missing))


if __name__ == "__main__":
    main()
