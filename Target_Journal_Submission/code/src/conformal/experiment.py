"""Per-cell conformal experiment: evaluate all methods on (dataset, model, feat, seed, split)."""
from __future__ import annotations
import numpy as np
from data import datasets as D
from data.splits import three_way, murcko_scaffold
from preprocessing.cache import features, feature_spec_to_pipe
from models.zoo import build_model
from conformal.methods import (ad_distance, marginal, normalized, mondrian, smoothed_mondrian,
                               coverage, group_coverage, conditional_gap, conformal_quantile)
import lightgbm as lgb

N_EVAL_BINS = 5  # evaluation AD bins (quintiles of calibration distance)


def _ad_bins(dist_cal, dist_te, n=N_EVAL_BINS):
    edges = np.quantile(dist_cal, np.linspace(0, 1, n + 1)); edges[0], edges[-1] = -np.inf, np.inf
    gcal = np.clip(np.digitize(dist_cal, edges[1:-1]), 0, n - 1)
    gte = np.clip(np.digitize(dist_te, edges[1:-1]), 0, n - 1)
    return gcal, gte


def run_cell(spec, seed, model_name="lgbm", feat="ecfp", split_kind="scaffold", alpha=0.1):
    df = D.load(spec); raw = df.smiles.tolist(); y = df.y.to_numpy().astype(float)
    tr, cal, te = three_way(raw, split_kind, 0.7, 0.15, seed)
    X = features(spec.name, raw, feature_spec_to_pipe("P_rdkit"), feat)
    m = build_model(model_name, spec.task, seed); m.fit(X[tr], y[tr])
    res_cal = y[cal] - m.predict(X[cal]); pred_te = m.predict(X[te]); yte = y[te]
    tr_smi = [raw[i] for i in tr]
    d_cal = ad_distance([raw[i] for i in cal], tr_smi); d_te = ad_distance([raw[i] for i in te], tr_smi)
    gcal, gte = _ad_bins(d_cal, d_te)
    target = 1 - alpha
    # normalized sigma(x): LightGBM predicting |residual| from calibration features (heteroscedastic CP)
    sig_model = lgb.LGBMRegressor(n_estimators=200, num_leaves=15, min_child_samples=20, random_state=seed, n_jobs=1, verbosity=-1)
    sig_model.fit(X[cal], np.abs(res_cal))
    sig_cal = np.maximum(sig_model.predict(X[cal]), 1e-6); sig_te = np.maximum(sig_model.predict(X[te]), 1e-6)
    # scaffold groups among cal/test (for scaffold-Mondrian + scaffold-conditional eval)
    sc_cal = np.array([hash(murcko_scaffold(raw[i])) % 50 for i in cal])
    sc_te = np.array([hash(murcko_scaffold(raw[i])) % 50 for i in te])

    halves = {
        "marginal": marginal(res_cal, len(te), alpha),
        "normalized": normalized(res_cal, sig_cal, sig_te, alpha),
        "mondrian_ad": mondrian(res_cal, gcal, gte, alpha),
        "smoothed_ad": smoothed_mondrian(res_cal, d_cal, d_te, alpha, n_bins=N_EVAL_BINS, shrink=0.3),
        "mondrian_scaffold": mondrian(res_cal, sc_cal, sc_te, alpha),
    }
    rows = []
    scale = float(np.std(y[tr])) + 1e-9
    for method, h in halves.items():
        gc = group_coverage(yte, pred_te, h, gte)
        far = gc.get(N_EVAL_BINS - 1, float("nan"))   # most-novel bin
        near = gc.get(0, float("nan"))
        worst = min(gc.values()) if gc else float("nan")
        row = {"dataset": spec.name, "task": spec.task, "model": model_name, "feat": feat,
               "seed": seed, "split": split_kind, "alpha": alpha, "method": method,
               "marginal_cov": coverage(yte, pred_te, h),
               "far_cov": far, "near_cov": near, "worst_bin_cov": worst,
               "cond_gap_ad": conditional_gap(yte, pred_te, h, gte, target),
               "cond_gap_scaffold": conditional_gap(yte, pred_te, h, sc_te, target),
               "mean_width": float(np.mean(2 * h)), "rel_width": float(np.mean(2 * h) / scale),
               "far_width": float(np.mean(2 * h[gte == N_EVAL_BINS - 1])) if (gte == N_EVAL_BINS - 1).any() else float("nan"),
               "n_test": len(te), "n_cal": len(cal)}
        rows.append(row)
    return rows
