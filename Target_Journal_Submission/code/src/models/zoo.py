"""CPU model zoo for molecular property prediction. Each entry builds a fresh sklearn-compatible
estimator given a seed and task type ('reg' or 'clf'). Featurization is handled separately
(src/preprocessing) so models consume dense feature matrices.
"""
from __future__ import annotations
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from sklearn.svm import LinearSVR, LinearSVC
from sklearn.kernel_approximation import Nystroem
from sklearn.linear_model import Ridge, LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
import lightgbm as lgb

N_JOBS = 1  # per-model single-thread; parallelism is across the experiment grid


def build_model(name: str, task: str, seed: int):
    """Return a fresh estimator. name in {lgbm, rf, svm, linear}."""
    if name == "lgbm":
        common = dict(n_estimators=400, learning_rate=0.05, num_leaves=63, subsample=0.8,
                      subsample_freq=1, colsample_bytree=0.8, min_child_samples=20,
                      random_state=seed, n_jobs=N_JOBS, verbosity=-1)
        return lgb.LGBMRegressor(**common) if task == "reg" else lgb.LGBMClassifier(**common)
    if name == "rf":
        # max_features='sqrt' is essential for high-dim sparse fingerprints: the sklearn default
        # (1.0 for regression) evaluates all 2048 features per split and is ~100x slower.
        common = dict(n_estimators=200, max_depth=None, min_samples_leaf=2, max_features="sqrt",
                      random_state=seed, n_jobs=N_JOBS)
        return RandomForestRegressor(**common) if task == "reg" else RandomForestClassifier(**common)
    if name == "svm":
        # RBF-kernel SVM via Nystroem approximation -> scales linearly in n (consistent across all
        # dataset sizes, unlike exact SVC/SVR which is O(n^2) and blows up on 10k-molecule sets).
        nys = Nystroem(kernel="rbf", gamma=None, n_components=512, random_state=seed)
        if task == "reg":
            return make_pipeline(StandardScaler(with_mean=False), nys, LinearSVR(C=1.0, random_state=seed, max_iter=5000))
        # probability via decision_function+sigmoid in predict_scores
        return make_pipeline(StandardScaler(with_mean=False), nys, LinearSVC(C=1.0, random_state=seed, max_iter=5000))
    if name == "linear":
        if task == "reg":
            return make_pipeline(StandardScaler(with_mean=False), Ridge(alpha=1.0, random_state=seed))
        return make_pipeline(StandardScaler(with_mean=False),
                             LogisticRegression(C=1.0, max_iter=2000, random_state=seed))
    raise ValueError(name)


MODEL_NAMES = ["lgbm", "rf", "svm", "linear"]


def predict_scores(model, X, task: str):
    """Unified prediction: regression -> values; classification -> P(y=1)."""
    if task == "reg":
        return model.predict(X)
    if hasattr(model, "predict_proba"):
        return model.predict_proba(X)[:, 1]
    # fall back to decision_function squashed
    import numpy as np
    d = model.decision_function(X)
    return 1.0 / (1.0 + np.exp(-d))
