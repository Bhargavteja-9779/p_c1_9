"""Statistical comparison utilities for multi-dataset, multi-seed ML experiments.

Implements the ingredients of a practically-significant method-comparison protocol:
- paired Wilcoxon signed-rank test (per-dataset, across seeds/folds)
- bootstrap confidence intervals for a statistic and for paired differences
- effect sizes (Cohen's d_z for paired data; rank-biserial for Wilcoxon)
- Friedman omnibus test + Nemenyi critical difference across many datasets
- Holm-Bonferroni multiple-comparison correction
These are standard, well-established tests; the module selects and reports them, it does not invent significance.
"""
from __future__ import annotations
import numpy as np
from scipy import stats

def paired_wilcoxon(a, b):
    """Two-sided Wilcoxon signed-rank on paired samples a,b. Returns dict with stat, p, rank-biserial effect size, n."""
    a = np.asarray(a, float); b = np.asarray(b, float)
    d = a - b
    nz = d[d != 0]
    if len(nz) == 0:
        return {"stat": float("nan"), "p": 1.0, "rank_biserial": 0.0, "n": int(len(a)), "n_nonzero": 0}
    res = stats.wilcoxon(a, b, zero_method="wilcox", alternative="two-sided", method="auto")
    # rank-biserial correlation effect size
    r = res.statistic
    n = len(nz); total = n * (n + 1) / 2
    rb = 1.0 - 2.0 * r / total
    return {"stat": float(res.statistic), "p": float(res.pvalue), "rank_biserial": float(rb),
            "n": int(len(a)), "n_nonzero": int(n)}

def cohen_dz(a, b) -> float:
    """Cohen's d_z for paired samples."""
    d = np.asarray(a, float) - np.asarray(b, float)
    sd = d.std(ddof=1)
    return float(d.mean() / sd) if sd > 0 else 0.0

def bootstrap_ci(x, func=np.mean, n_boot: int = 10_000, alpha: float = 0.05, seed: int = 0):
    """Percentile bootstrap CI for func(x)."""
    x = np.asarray(x, float); rng = np.random.default_rng(seed)
    n = len(x)
    if n == 0:
        return (float("nan"), float("nan"), float("nan"))
    stats_ = np.array([func(x[rng.integers(0, n, n)]) for _ in range(n_boot)])
    lo, hi = np.percentile(stats_, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    return (float(func(x)), float(lo), float(hi))

def paired_diff_ci(a, b, n_boot: int = 10_000, alpha: float = 0.05, seed: int = 0):
    """Bootstrap CI of the mean paired difference a-b."""
    return bootstrap_ci(np.asarray(a, float) - np.asarray(b, float), np.mean, n_boot, alpha, seed)

def holm_bonferroni(pvals: dict, alpha: float = 0.05) -> dict:
    """Holm-Bonferroni correction. Returns {key: (p, reject)} sorted by ascending p."""
    items = sorted(pvals.items(), key=lambda kv: kv[1])
    m = len(items); out = {}
    for i, (k, p) in enumerate(items):
        thresh = alpha / (m - i)
        out[k] = (float(p), bool(p <= thresh))
        if p > thresh:
            for k2, p2 in items[i + 1:]:
                out[k2] = (float(p2), False)
            break
    return out

def friedman(score_matrix):
    """Friedman omnibus test. Rows = datasets, cols = methods. Returns stat, p, mean ranks (1=best for larger-is-better).
    Pass scores where LARGER IS BETTER."""
    M = np.asarray(score_matrix, float)
    stat, p = stats.friedmanchisquare(*[M[:, j] for j in range(M.shape[1])])
    ranks = np.array([stats.rankdata(-row) for row in M])  # 1 = best
    return {"stat": float(stat), "p": float(p), "mean_ranks": ranks.mean(0).tolist(), "n_datasets": M.shape[0], "n_methods": M.shape[1]}

def nemenyi_cd(n_methods: int, n_datasets: int, alpha: float = 0.05) -> float:
    """Nemenyi critical difference for average ranks."""
    # studentized range q_alpha / sqrt(2) tabulated for alpha=0.05, 0.10
    q05 = {2:1.960,3:2.343,4:2.569,5:2.728,6:2.850,7:2.949,8:3.031,9:3.102,10:3.164,
           11:3.219,12:3.268,13:3.313,14:3.354,15:3.391}
    q10 = {2:1.645,3:2.052,4:2.291,5:2.459,6:2.589,7:2.693,8:2.780,9:2.855,10:2.920}
    q = (q05 if alpha == 0.05 else q10).get(n_methods)
    if q is None:
        return float("nan")
    return float(q * np.sqrt(n_methods * (n_methods + 1) / (6.0 * n_datasets)))
