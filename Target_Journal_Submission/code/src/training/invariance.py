"""Standardization-invariance measurement and the SORB method evaluation.

Primary invariance measure: for a fitted model, present each test molecule under each of the
realistic named pipelines (P_raw..P_aggressive) and measure the spread of the model's prediction
using the cached per-pipeline feature matrices. Because the modeled property is invariant to these
chemically-equivalent standardization choices, any spread is spurious model sensitivity.

Method (SORB) variants:
  baseline       - trained canonicalize-once; single-presentation inference        (the status quo)
  sorb_trainaug  - trained on orbit augmentation; single-presentation inference     (LEARNED invariance)
  sorb_testorbit - baseline model; inference averaged over the orbit                (ENFORCED invariance)
  sorb_full      - trainaug model; inference averaged over the orbit
"""
from __future__ import annotations
import numpy as np
from data import datasets as D
from data.splits import scaffold_split, random_split
from preprocessing.cache import features, feature_spec_to_pipe, dataset_orbits
from preprocessing.standardize import NAMED_PIPELINES, OrbitConfig
from models.zoo import build_model, predict_scores
from evaluation.metrics import regression_metrics, classification_metrics, expected_calibration_error

INV_PIPELINES = list(NAMED_PIPELINES.keys())
ACC_CAP = 600  # cap test molecules for orbit-averaged accuracy to bound compute on large sets


def _split(raw, kind, seed):
    return scaffold_split(raw, 0.8, seed) if kind == "scaffold" else random_split(len(raw), 0.8, seed)


def _acc(task, y, s):
    if task == "reg":
        m = regression_metrics(y, s); return {"rmse": m["rmse"], "spearman": m["spearman"], "mae": m["mae"]}
    m = classification_metrics(y, s)
    return {"auroc": m["auroc"], "auprc": m["auprc"], "mcc": m["mcc"], "ece": expected_calibration_error(y, s)}


def _inv_stats(pred_matrix, scale):
    sd = pred_matrix.std(axis=1)
    rng = pred_matrix.max(axis=1) - pred_matrix.min(axis=1)
    scale = scale if scale and scale > 0 else 1.0
    return {"mean_rel_sd": float(sd.mean() / scale), "median_rel_sd": float(np.median(sd) / scale),
            "p90_rel_sd": float(np.percentile(sd, 90) / scale), "mean_rel_range": float(rng.mean() / scale),
            "frac_gt5pct": float(np.mean(sd > 0.05 * scale)), "frac_changed": float(np.mean(rng > 1e-9))}


def _pipeline_matrix(model, spec, raw, te, feat, task):
    return np.column_stack([
        predict_scores(model, features(spec.name, raw, feature_spec_to_pipe(p), feat)[te], task)
        for p in INV_PIPELINES])


def _orbit_invariance(model, orbit_subset, feat, task, scale):
    """Per-molecule prediction SD across its (cached) standardization orbit -> the pure invariance
    test (all forms are the same molecule). orbit_subset: list of SMILES-lists. Fast: cached SMILES."""
    from preprocessing.featurize import ecfp_matrix, rdkit_descriptors
    sds = []
    for orb in orbit_subset:
        if not orb or len(orb) < 2:
            sds.append(0.0); continue
        Xo = ecfp_matrix(orb, 2, 2048) if feat == "ecfp" else rdkit_descriptors(orb)[0]
        sds.append(float(np.std(predict_scores(model, Xo, task))))
    sds = np.array(sds); scale = scale if scale and scale > 0 else 1.0
    return {"orbit_mean_rel_sd": float(sds.mean() / scale),
            "orbit_p90_rel_sd": float(np.percentile(sds, 90) / scale),
            "orbit_frac_gt5pct": float(np.mean(sds > 0.05 * scale))}


def run_method_cell(spec, seed, model_name="lgbm", feat="ecfp", split_kind="scaffold",
                    train_aug=4, test_orbit=6, orbit_cfg: OrbitConfig | None = None):
    from models.sorb import SORB
    cfg = orbit_cfg or OrbitConfig()
    df = D.load(spec); raw = df.smiles.tolist(); y = df.y.to_numpy().astype(float)
    tr, te = _split(raw, split_kind, seed)
    scale = float(np.std(y[tr])) if spec.task == "reg" else 1.0
    orbits = dataset_orbits(spec.name, raw, cfg)
    tr_smi = [raw[i] for i in tr]; tr_orb = [orbits[i] for i in tr]
    Xcanon = features(spec.name, raw, feature_spec_to_pipe("P_rdkit"), feat)
    rng = np.random.default_rng(seed)
    cap = te if len(te) <= ACC_CAP else te[rng.choice(len(te), ACC_CAP, replace=False)]
    cap_smi = [raw[i] for i in cap]; cap_orb = [orbits[i] for i in cap]
    inv_cap = te if len(te) <= 300 else te[rng.choice(len(te), 300, replace=False)]
    inv_orb = [orbits[i] for i in inv_cap]
    rows = []

    def base_row(method, acc_scores, inv_model):
        row = {"dataset": spec.name, "task": spec.task, "model": model_name, "feat": feat,
               "method": method, "seed": seed, "split": split_kind, "label_sd": scale,
               "n_train": len(tr), "n_test": len(te)}
        row.update(acc_scores)
        if inv_model is not None:
            row.update(_inv_stats(_pipeline_matrix(inv_model, spec, raw, te, feat, spec.task), scale))
            row.update(_orbit_invariance(inv_model, inv_orb, feat, spec.task, scale))
        rows.append(row)

    base = build_model(model_name, spec.task, seed); base.fit(Xcanon[tr], y[tr])
    base_row("baseline", _acc(spec.task, y[te], predict_scores(base, Xcanon[te], spec.task)), base)

    est_ta = SORB(model_name, spec.task, feat, seed=seed, orbit_cfg=cfg, train_aug=train_aug, test_orbit=1)
    est_ta.fit(tr_smi, y[tr], orbits=tr_orb)
    base_row("sorb_trainaug",
             _acc(spec.task, y[te], predict_scores(est_ta.model, Xcanon[te], spec.task)), est_ta.model)

    for method, model_for in [("sorb_testorbit", base), ("sorb_full", est_ta.model)]:
        est = SORB(model_name, spec.task, feat, seed=seed, orbit_cfg=cfg, train_aug=1, test_orbit=test_orbit)
        est.model = model_for
        sc = est.predict(cap_smi, orbits=cap_orb)
        base_row(method, _acc(spec.task, y[cap], sc), None)
    return rows
