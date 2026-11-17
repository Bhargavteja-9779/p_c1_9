"""Compute all RQ summaries and statistics from the result CSVs. Writes:
  results/metrics.json, results/statistics.json, tables/*.csv (auto-generated, machine-readable).
Every number cited in the manuscript originates here.
"""
from __future__ import annotations
import os, sys, json
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ))
import numpy as np
import pandas as pd
from scipy import stats as sps
from stat_tests.compare import paired_wilcoxon, cohen_dz, bootstrap_ci, paired_diff_ci, friedman, nemenyi_cd, holm_bonferroni

ROOT = os.path.join(os.path.dirname(__file__), "..", "..")
RES = os.path.join(ROOT, "results")
TAB = os.path.join(ROOT, "tables")
os.makedirs(TAB, exist_ok=True)

PRIMARY = {"reg": "rmse", "clf": "auroc"}            # primary accuracy metric per task
SECONDARY = {"reg": "spearman", "clf": "auprc"}
LOWER_BETTER = {"rmse": True, "mae": True, "ece": True, "auroc": False, "auprc": False,
                "spearman": False, "mcc": False}


# molnet:lipophilicity is the same AstraZeneca logD7.4 benchmark as tdc:lipophilicity_astrazeneca
# (confirmed identical by the leakage audit); drop the duplicate from all analyses.
DROP_DATASETS = {"molnet:lipophilicity"}


def _load(name):
    p = os.path.join(RES, name)
    if not os.path.exists(p):
        return None
    df = pd.read_csv(p)
    if "dataset" in df.columns:
        df = df[~df.dataset.isin(DROP_DATASETS)].copy()
    return df


def rq1_sensitivity(sens):
    """Spread of the primary metric across the 6 pipelines, per dataset x model x feat (mean over seeds)."""
    rows = []
    for (ds, task, model, feat), g in sens.groupby(["dataset", "task", "model", "feat"]):
        metric = PRIMARY[task]
        pm = g.groupby("pipeline")[metric].mean()
        if pm.notna().sum() < 2:
            continue
        spread = float(pm.max() - pm.min())
        ref = pm.get("P_rdkit", pm.mean())
        rows.append({"dataset": ds, "task": task, "model": model, "feat": feat, "metric": metric,
                     "min": float(pm.min()), "max": float(pm.max()), "spread": spread,
                     "rel_spread_pct": float(100 * spread / (abs(ref) + 1e-9))})
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(TAB, "rq1_sensitivity.csv"), index=False)
    summary = {"n_cells": len(df),
               "median_rel_spread_pct": float(df.rel_spread_pct.median()),
               "p90_rel_spread_pct": float(df.rel_spread_pct.quantile(0.9)),
               "by_model_feat": df.groupby(["model", "feat"]).rel_spread_pct.median().round(3).to_dict().__str__()}
    return df, summary


def rq2_ranking(sens):
    """Model-ranking stability across pipelines: Kendall-tau between pipeline rankings, top-model change."""
    rows = []
    for (ds, task, feat), g in sens.groupby(["dataset", "task", "feat"]):
        metric = PRIMARY[task]; sign = 1 if LOWER_BETTER[metric] else -1
        piv = g.groupby(["pipeline", "model"])[metric].mean().unstack()
        if piv.shape[0] < 2 or piv.shape[1] < 3 or piv.isna().any().any():
            continue
        pipes = list(piv.index)
        rankings = {p: sps.rankdata(sign * piv.loc[p].values) for p in pipes}
        taus = []
        for i in range(len(pipes)):
            for j in range(i + 1, len(pipes)):
                t = sps.kendalltau(rankings[pipes[i]], rankings[pipes[j]]).statistic
                if not np.isnan(t):
                    taus.append(t)
        best = {p: piv.loc[p].idxmin() if LOWER_BETTER[metric] else piv.loc[p].idxmax() for p in pipes}
        top_change = len(set(best.values())) > 1
        rows.append({"dataset": ds, "task": task, "feat": feat, "mean_kendall_tau": float(np.mean(taus)),
                     "min_kendall_tau": float(np.min(taus)), "top_model_changes": bool(top_change),
                     "n_distinct_top_models": len(set(best.values()))})
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(TAB, "rq2_ranking.csv"), index=False)
    summary = {"n_datasets": len(df), "median_mean_kendall_tau": float(df.mean_kendall_tau.median()),
               "frac_top_model_changes": float(df.top_model_changes.mean()),
               "frac_min_tau_lt_1": float((df.min_kendall_tau < 1.0).mean())}
    return df, summary


def rq3_ablation(abl):
    """Effect of each single standardization operation vs base (mean primary-metric change)."""
    rows = []
    for (ds, task, model), g in abl.groupby(["dataset", "task", "model"]):
        metric = PRIMARY[task]
        base = g[g.ablation_pipe == "base"][metric].mean()
        for op in g.ablation_pipe.unique():
            if op == "base":
                continue
            v = g[g.ablation_pipe == op][metric].mean()
            rows.append({"dataset": ds, "task": task, "model": model, "operation": op,
                         "delta_vs_base": float(v - base),
                         "abs_delta": float(abs(v - base))})
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(TAB, "rq3_ablation.csv"), index=False)
    summary = df.groupby("operation").abs_delta.median().round(4).to_dict()
    return df, {"median_abs_delta_by_operation": summary}


def rq4_method(meth):
    """SORB variants vs baseline: paired deltas of accuracy across datasets (per model x feat)."""
    out = {}; tables = []
    for (model, feat, task), g in meth.groupby(["model", "feat", "task"]):
        metric = PRIMARY[task]; lower = LOWER_BETTER[metric]
        wide = g.groupby(["dataset", "seed", "method"])[metric].mean().unstack()
        wide = wide.dropna(subset=["baseline"])
        # average over seeds per dataset
        ds_mean = wide.groupby(level=0).mean()
        for variant in ["sorb_trainaug", "sorb_testorbit", "sorb_full"]:
            if variant not in ds_mean.columns:
                continue
            pair = ds_mean[["baseline", variant]].dropna()
            if len(pair) < 3:
                continue
            delta = (pair[variant] - pair["baseline"]).values  # positive = variant higher metric
            improve = -delta if lower else delta                # positive = improvement
            w = paired_wilcoxon(pair[variant].values, pair["baseline"].values)
            ci = bootstrap_ci(improve, np.mean, n_boot=5000)
            tables.append({"model": model, "feat": feat, "task": task, "variant": variant, "metric": metric,
                           "n_datasets": len(pair), "mean_improvement": float(improve.mean()),
                           "median_improvement": float(np.median(improve)),
                           "frac_datasets_improved": float(np.mean(improve > 0)),
                           "wilcoxon_p": w["p"], "cohen_dz": cohen_dz(improve, np.zeros_like(improve)),
                           "ci_lo": ci[1], "ci_hi": ci[2]})
    df = pd.DataFrame(tables)
    df.to_csv(os.path.join(TAB, "rq4_method_accuracy.csv"), index=False)
    return df, {"n_rows": len(df)}


def rq5_invariance(meth):
    """Invariance violation (baseline) and SORB train-aug reduction, by model x feat (representation dep.)."""
    base = meth[meth.method == "baseline"]
    ta = meth[meth.method == "sorb_trainaug"]
    rows = []
    for (model, feat), gb in base.groupby(["model", "feat"]):
        gt = ta[(ta.model == model) & (ta.feat == feat)]
        # per dataset mean over seeds
        for col in ["orbit_mean_rel_sd", "mean_rel_sd"]:
            b = gb.groupby("dataset")[col].mean().dropna()
            t = gt.groupby("dataset")[col].mean().dropna()
            common = b.index.intersection(t.index)
            if len(common) < 3:
                continue
            red = (b[common] - t[common]).values
            w = paired_wilcoxon(b[common].values, t[common].values)
            rows.append({"model": model, "feat": feat, "measure": col,
                         "baseline_median": float(b.median()), "trainaug_median": float(t.median()),
                         "median_reduction": float(np.median(red)),
                         "frac_reduced": float(np.mean(red > 0)), "wilcoxon_p": w["p"], "n_datasets": len(common)})
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(TAB, "rq5_invariance.csv"), index=False)
    # headline baseline invariance by representation. Regression only: orbit/pipe SD are divided by the
    # label SD, which is meaningful only for regression. Classification probability SD is reported
    # separately (raw, on the 0-1 probability scale) to avoid the scale conflation flagged in review.
    reg = base[base.task == "reg"]
    inv_head = reg.groupby(["model", "feat"]).agg(
        orbit_mean_rel_sd=("orbit_mean_rel_sd", "median"),
        pipe_mean_rel_sd=("mean_rel_sd", "median")).round(4).reset_index()
    inv_head.to_csv(os.path.join(TAB, "rq5_baseline_invariance_by_representation.csv"), index=False)
    clf = base[base.task == "clf"]
    clf_head = clf.groupby(["model", "feat"]).agg(
        orbit_prob_sd=("orbit_mean_rel_sd", "median"),   # label_sd was 1.0 for clf -> this is raw prob SD
        pipe_prob_sd=("mean_rel_sd", "median")).round(4).reset_index()
    clf_head.to_csv(os.path.join(TAB, "rq5_baseline_invariance_classification.csv"), index=False)
    return df, {"baseline_invariance_by_representation": inv_head.to_dict("records"),
                "classification_probability_sd": clf_head.to_dict("records"),
                "n_reg_datasets": int(reg.dataset.nunique()), "n_clf_datasets": int(clf.dataset.nunique())}



def rq5b_shortcut(sc):
    """Semi-synthetic shortcut: baseline exploits standardization; SORB resists."""
    rows = []
    for (model, feat, task), g in sc.groupby(["model", "feat", "task"]):
        gd = g.groupby("dataset").mean(numeric_only=True)
        rows.append({"model": model, "feat": feat, "task": task, "n_datasets": len(gd),
                     "median_baseline_shortcut_drop": float(gd.baseline_shortcut_drop.median()),
                     "median_baseline_cf_vs_honest": float(gd.baseline_cf_vs_honest.median()),
                     "median_sorb_cf_vs_honest": float(gd.sorb_cf_vs_honest.median())})
    df = pd.DataFrame(rows); df.to_csv(os.path.join(TAB, "rq5b_shortcut.csv"), index=False)
    pooled = sc.groupby("dataset").mean(numeric_only=True)
    out = {"median_baseline_shortcut_drop": float(pooled.baseline_shortcut_drop.median()),
           "median_baseline_cf_vs_honest": float(pooled.baseline_cf_vs_honest.median()),
           "n_datasets": int(len(pooled))}
    if "sorb_shortcut_drop" in pooled:
        out["median_sorb_shortcut_drop"] = float(pooled.sorb_shortcut_drop.median())
        # how much more the baseline relies on the shortcut than SORB (paired)
        pw = paired_wilcoxon(sc.groupby("dataset").baseline_shortcut_drop.mean().values,
                             sc.groupby("dataset").sorb_shortcut_drop.mean().values)
        out["baseline_vs_sorb_shortcut_wilcoxon_p"] = float(pw["p"])
    return df, out


def method_pooled(meth):
    """Pooled (over model x feat x dataset) accuracy effect per SORB variant, and train-aug invariance
    reduction. Positive improvement = better; lower_better handled per task."""
    out = {}
    for task in meth.task.unique():
        g = meth[meth.task == task]; metric = PRIMARY[task]; lower = LOWER_BETTER[metric]
        base = g[g.method == "baseline"].groupby(["dataset", "model", "feat"])[metric].mean()
        for v in ["sorb_trainaug", "sorb_testorbit", "sorb_full"]:
            vv = g[g.method == v].groupby(["dataset", "model", "feat"])[metric].mean()
            d = (base - vv).dropna() if lower else (vv - base).dropna()  # +improvement
            if len(d):
                w = paired_wilcoxon(vv.loc[d.index].values, base.loc[d.index].values)
                out[f"{task}_{v}"] = {"median_improvement": float(d.median()),
                                      "frac_improved": float((d > 0).mean()),
                                      "wilcoxon_p": float(w["p"]), "n": int(len(d))}
    # train-aug invariance reduction (orbit SD), pooled
    b = meth[meth.method == "baseline"].groupby(["dataset", "model", "feat"]).orbit_mean_rel_sd.mean()
    t = meth[meth.method == "sorb_trainaug"].groupby(["dataset", "model", "feat"]).orbit_mean_rel_sd.mean()
    j = pd.concat([b.rename("b"), t.rename("t")], axis=1).dropna()
    if len(j):
        red = (j.b - j.t)
        out["trainaug_invariance"] = {"baseline_median": float(j.b.median()), "trainaug_median": float(j.t.median()),
                                      "median_abs_reduction": float(red.median()),
                                      "median_rel_reduction_pct": float(100 * (red / j.b.replace(0, np.nan)).median()),
                                      "frac_reduced": float((red > 0).mean())}
    return out


def method_omnibus(meth):
    """Friedman omnibus + Nemenyi CD across regression datasets over the 4 methods (per model x feat,
    then pooled by averaging to one value per dataset), plus Holm-Bonferroni on the per-cell SORB-vs-
    baseline Wilcoxon p-values. Addresses the review that these were claimed but not computed."""
    out = {}
    reg = meth[meth.task == "reg"]
    methods = ["baseline", "sorb_trainaug", "sorb_testorbit", "sorb_full"]
    # one RMSE per (dataset) averaged over model/feat/seed, per method
    piv = (reg.groupby(["dataset", "method"]).rmse.mean().unstack())
    piv = piv.dropna(subset=methods)[methods]
    if len(piv) >= 3:
        # RMSE: lower is better -> pass negated so friedman() ranks 1=best
        fr = friedman(-piv.values)
        out["friedman"] = {"p": fr["p"], "stat": fr["stat"], "mean_ranks": dict(zip(methods, fr["mean_ranks"])),
                           "nemenyi_cd": nemenyi_cd(len(methods), len(piv)), "n_datasets": len(piv)}
    # Holm-Bonferroni across the per-cell SORB-vs-baseline tests
    r4 = _load("../tables/rq4_method_accuracy.csv") if False else None
    import os as _os
    p = _os.path.join(TAB, "rq4_method_accuracy.csv")
    if _os.path.exists(p):
        t = pd.read_csv(p)
        pv = {f"{r.model}|{r.feat}|{r.task}|{r.variant}": float(r.wilcoxon_p) for r in t.itertuples()}
        hb = holm_bonferroni(pv)
        n_sig_raw = int((t.wilcoxon_p < 0.05).sum())
        n_sig_adj = int(sum(1 for v in hb.values() if v[1]))
        out["multiplicity"] = {"n_tests": len(pv), "n_sig_raw_0.05": n_sig_raw, "n_sig_holm_0.05": n_sig_adj}
        t["holm_p_reject"] = t.apply(lambda r: hb.get(f"{r.model}|{r.feat}|{r.task}|{r.variant}", (1, False))[1], axis=1)
        t.to_csv(p, index=False)
    return out


def main():
    metrics, statistics = {}, {}
    sens = _load("raw_results.csv")
    abl = _load("ablation_results.csv")
    meth = _load("method_results.csv")
    if sens is not None:
        _, metrics["RQ1_sensitivity"] = rq1_sensitivity(sens)
        _, metrics["RQ2_ranking"] = rq2_ranking(sens)
    if abl is not None:
        _, metrics["RQ3_ablation"] = rq3_ablation(abl)
    if meth is not None:
        dfm, _ = rq4_method(meth); statistics["RQ4_method_table"] = dfm.to_dict("records")
        dfi, metrics["RQ5_invariance"] = rq5_invariance(meth)
        statistics["RQ5_invariance_table"] = dfi.to_dict("records")
        metrics["RQ4_method_pooled"] = method_pooled(meth)
        metrics["RQ4_omnibus"] = method_omnibus(meth)
    sc = _load("shortcut_results.csv")
    if sc is not None:
        _, metrics["RQ5b_shortcut"] = rq5b_shortcut(sc)
    json.dump(metrics, open(os.path.join(RES, "metrics.json"), "w"), indent=1, default=float)
    json.dump(statistics, open(os.path.join(RES, "statistics.json"), "w"), indent=1, default=float)
    print(json.dumps(metrics, indent=1, default=float)[:3000])


if __name__ == "__main__":
    main()
