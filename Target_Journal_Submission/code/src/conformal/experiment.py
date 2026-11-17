"""Per-cell conformal experiment: evaluate all methods on (dataset, model, feat, seed, split)."""
from __future__ import annotations
import numpy as np
from data import datasets as D
from data.splits import three_way, murcko_scaffold
from preprocessing.cache import features, feature_spec_to_pipe
from models.zoo import build_model
from conformal.methods import (ad_distance, marginal, normalized, distance_normalized, mondrian,
                               smoothed_mondrian, coverage, group_coverage, conditional_gap, conformal_quantile)
import lightgbm as lgb
import hashlib


def _scaffold_bucket(smiles, nbuckets=50):
    """Deterministic scaffold -> bucket id (hashlib, not Python's salted hash())."""
    sc = murcko_scaffold(smiles) or ""
    return int(hashlib.md5(sc.encode("utf-8")).hexdigest(), 16) % nbuckets

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
    # normalized sigma(x): a heteroscedastic difficulty model for normalized conformal. CRITICAL for
    # validity: sigma must be independent of the calibration/test labels, else the normalized scores are
    # not exchangeable and marginal coverage breaks. We therefore train sigma on OUT-OF-FOLD residuals of
    # the TRAINING set (2-fold cross-fit), disjoint from cal/test, giving a genuine (non-degenerate)
    # difficulty estimate while preserving the marginal guarantee. sigma_cal/sigma_te are out-of-sample.
    tr_arr = np.asarray(tr); Xtr = X[tr_arr]; ytr = y[tr_arr]
    oof_res = np.empty(len(tr_arr))
    _perm = np.random.default_rng(1000 + seed).permutation(len(tr_arr))
    _folds = np.array_split(_perm, 2)
    for _k in range(2):
        _val = _folds[_k]; _trn = np.concatenate([_folds[_j] for _j in range(2) if _j != _k])
        if len(_trn) < 2 or len(_val) == 0:
            oof_res[_val] = np.abs(res_cal).mean() if len(res_cal) else 0.0
            continue
        _bk = build_model(model_name, spec.task, seed); _bk.fit(Xtr[_trn], ytr[_trn])
        oof_res[_val] = np.abs(ytr[_val] - _bk.predict(Xtr[_val]))
    sig_model = lgb.LGBMRegressor(n_estimators=200, num_leaves=15, min_child_samples=20, random_state=seed, n_jobs=1, verbosity=-1)
    sig_model.fit(Xtr, oof_res)
    sig_cal = np.maximum(sig_model.predict(X[cal]), 1e-6); sig_te = np.maximum(sig_model.predict(X[te]), 1e-6)
    beta_norm = float(np.median(np.maximum(oof_res, 0.0)))  # normalized regularizer from OOF train (cal-independent)
    if not np.isfinite(beta_norm) or beta_norm <= 0:
        beta_norm = None
    # scaffold groups among cal/test (deterministic hash; for scaffold-Mondrian + scaffold-conditional eval)
    sc_cal = np.array([_scaffold_bucket(raw[i]) for i in cal])
    sc_te = np.array([_scaffold_bucket(raw[i]) for i in te])

    halves = {
        "marginal": marginal(res_cal, len(te), alpha),
        "normalized": normalized(res_cal, sig_cal, sig_te, alpha, beta=beta_norm),
        "knn_normalized": distance_normalized(res_cal, d_cal, d_te, alpha),
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
