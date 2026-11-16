"""Regression and classification metrics, plus calibration/UQ diagnostics."""
from __future__ import annotations
import numpy as np
from scipy import stats
from sklearn import metrics as skm

# ---- regression ----
def regression_metrics(y_true, y_pred) -> dict:
    y_true = np.asarray(y_true, float); y_pred = np.asarray(y_pred, float)
    out = {
        "rmse": float(np.sqrt(np.mean((y_true - y_pred) ** 2))),
        "mae": float(np.mean(np.abs(y_true - y_pred))),
        "medae": float(np.median(np.abs(y_true - y_pred))),
        "r2": float(skm.r2_score(y_true, y_pred)) if len(y_true) > 1 else float("nan"),
        "pearson": float(stats.pearsonr(y_true, y_pred)[0]) if len(y_true) > 1 else float("nan"),
        "spearman": float(stats.spearmanr(y_true, y_pred)[0]) if len(y_true) > 1 else float("nan"),
        "kendall": float(stats.kendalltau(y_true, y_pred)[0]) if len(y_true) > 1 else float("nan"),
    }
    return out

# ---- classification (binary, needs scores/probabilities) ----
def classification_metrics(y_true, y_score, threshold: float = 0.5) -> dict:
    y_true = np.asarray(y_true, int); y_score = np.asarray(y_score, float)
    y_hat = (y_score >= threshold).astype(int)
    out = {
        "auroc": float(skm.roc_auc_score(y_true, y_score)) if len(np.unique(y_true)) > 1 else float("nan"),
        "auprc": float(skm.average_precision_score(y_true, y_score)) if len(np.unique(y_true)) > 1 else float("nan"),
        "balanced_acc": float(skm.balanced_accuracy_score(y_true, y_hat)),
        "mcc": float(skm.matthews_corrcoef(y_true, y_hat)),
        "f1": float(skm.f1_score(y_true, y_hat, zero_division=0)),
        "brier": float(np.mean((y_score - y_true) ** 2)),
        "prevalence": float(np.mean(y_true)),
    }
    return out

# ---- calibration ----
def expected_calibration_error(y_true, y_prob, n_bins: int = 15) -> float:
    """ECE with equal-width probability bins."""
    y_true = np.asarray(y_true, float); y_prob = np.asarray(y_prob, float)
    bins = np.linspace(0, 1, n_bins + 1)
    idx = np.clip(np.digitize(y_prob, bins) - 1, 0, n_bins - 1)
    ece = 0.0; n = len(y_true)
    for b in range(n_bins):
        m = idx == b
        if m.any():
            ece += (m.sum() / n) * abs(y_prob[m].mean() - y_true[m].mean())
    return float(ece)

# ---- interval / UQ diagnostics ----
def interval_coverage(y_true, lower, upper) -> float:
    y_true = np.asarray(y_true, float)
    return float(np.mean((y_true >= np.asarray(lower, float)) & (y_true <= np.asarray(upper, float))))

def mean_interval_width(lower, upper) -> float:
    return float(np.mean(np.asarray(upper, float) - np.asarray(lower, float)))

def error_ranking_auc(abs_err, uncertainty) -> float:
    """AUC-like: probability that a higher-uncertainty prediction has higher abs error
    (Spearman between uncertainty and |error|; robust to scale). 1 = perfect ranking."""
    abs_err = np.asarray(abs_err, float); uncertainty = np.asarray(uncertainty, float)
    if len(abs_err) < 3 or np.allclose(uncertainty, uncertainty[0]):
        return float("nan")
    return float(stats.spearmanr(uncertainty, abs_err)[0])
