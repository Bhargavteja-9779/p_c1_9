"""Conformal prediction methods and applicability-domain (AD) tooling for molecular regression.

Split (inductive) conformal gives finite-sample MARGINAL coverage under exchangeability. We compare:
  - marginal        : one quantile of |residual| over the whole calibration set (the standard method)
  - normalized      : heteroscedastic score |residual| / sigma(x) with a learned difficulty sigma(x)
  - mondrian(groups): a separate quantile per group; with group = AD-distance bin this yields
                      finite-sample group-conditional coverage (our method), with group = scaffold cluster
                      it yields scaffold-conditional coverage.
AD descriptors measure distance of a molecule to the training set (1 - mean top-k Tanimoto similarity).
"""
from __future__ import annotations
import numpy as np
from rdkit import DataStructs
from preprocessing.featurize import ecfp_bitset


# ---- applicability-domain descriptors --------------------------------------------------------
def ad_distance(query_smiles, train_smiles, k=5, n_bits=2048, radius=2):
    """1 - mean top-k Tanimoto similarity of each query molecule to the training set (in [0,1])."""
    tr = [ecfp_bitset(s, radius, n_bits) for s in train_smiles]
    tr = [f for f in tr if f is not None]
    out = np.ones(len(query_smiles))
    for i, s in enumerate(query_smiles):
        f = ecfp_bitset(s, radius, n_bits)
        if f is None or not tr:
            continue
        sims = sorted(DataStructs.BulkTanimotoSimilarity(f, tr), reverse=True)[:k]
        out[i] = 1.0 - float(np.mean(sims))
    return out


def conformal_quantile(scores, alpha):
    """The (1-alpha) conformal quantile with finite-sample correction. Returns +inf if too few points."""
    n = len(scores)
    if n == 0:
        return np.inf
    level = np.ceil((n + 1) * (1 - alpha)) / n
    if level > 1:
        return np.inf
    return float(np.quantile(scores, level, method="higher"))


# ---- conformal methods: return per-test half-width q_i (interval = pred ± q_i) ----------------
def marginal(res_cal, pred_te_dummy_len, alpha):
    q = conformal_quantile(np.abs(res_cal), alpha)
    return np.full(pred_te_dummy_len, q)


def normalized(res_cal, sigma_cal, sigma_te, alpha):
    """Heteroscedastic: score = |res|/sigma; half-width_i = q * sigma_i."""
    sc = np.abs(res_cal) / np.maximum(sigma_cal, 1e-6)
    q = conformal_quantile(sc, alpha)
    return q * np.maximum(sigma_te, 1e-6)


def mondrian(res_cal, group_cal, group_te, alpha, fallback_q=None):
    """Separate conformal quantile per group. Groups with too few calibration points fall back to the
    pooled quantile (fallback_q), preserving validity in expectation over the pooled set."""
    if fallback_q is None:
        fallback_q = conformal_quantile(np.abs(res_cal), alpha)
    q_by_group = {}
    for g in np.unique(group_cal):
        q_by_group[g] = conformal_quantile(np.abs(res_cal[group_cal == g]), alpha)
    out = np.array([q_by_group.get(g, fallback_q) for g in group_te], float)
    out[~np.isfinite(out)] = fallback_q
    return out


def smoothed_mondrian(res_cal, dist_cal, dist_te, alpha, n_bins=5, shrink=0.3):
    """AD-conditional conformal with shrinkage: per-distance-bin quantile shrunk toward the pooled
    quantile by `shrink`, mitigating the sparse-bin instability of hard Mondrian. Monotone-safe."""
    pooled = conformal_quantile(np.abs(res_cal), alpha)
    edges = np.quantile(dist_cal, np.linspace(0, 1, n_bins + 1))
    edges[0], edges[-1] = -np.inf, np.inf
    gcal = np.clip(np.digitize(dist_cal, edges[1:-1]), 0, n_bins - 1)
    gte = np.clip(np.digitize(dist_te, edges[1:-1]), 0, n_bins - 1)
    qg = {}
    for g in range(n_bins):
        m = gcal == g
        qg[g] = conformal_quantile(np.abs(res_cal[m]), alpha) if m.sum() >= 20 else pooled
        if not np.isfinite(qg[g]):
            qg[g] = pooled
    out = np.array([(1 - shrink) * qg[g] + shrink * pooled for g in gte], float)
    return out


# ---- metrics ---------------------------------------------------------------------------------
def coverage(y, pred, half):
    return float(np.mean((y >= pred - half) & (y <= pred + half)))


def group_coverage(y, pred, half, groups):
    cov = ((y >= pred - half) & (y <= pred + half)).astype(float)
    return {int(g): float(cov[groups == g].mean()) for g in np.unique(groups) if (groups == g).any()}


def conditional_gap(y, pred, half, groups, target):
    gc = group_coverage(y, pred, half, groups)
    return float(max(abs(c - target) for c in gc.values())) if gc else float("nan")
