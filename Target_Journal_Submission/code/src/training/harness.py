"""Core experiment harness: train/evaluate one (dataset, pipeline, seed) over all models,
reusing a single pipeline-invariant split so that only featurization differs across pipelines.
Also measures SORB and standardization-invariance of predictions.
"""
from __future__ import annotations
import numpy as np
from data import datasets as D
from data.splits import scaffold_split, random_split
from preprocessing.cache import features, feature_spec_to_pipe
from preprocessing.standardize import NAMED_PIPELINES, standardization_orbit, OrbitConfig
from preprocessing.featurize import ecfp_matrix, rdkit_descriptors
from models.zoo import build_model, predict_scores, MODEL_NAMES
from evaluation.metrics import regression_metrics, classification_metrics, expected_calibration_error

ECFP_RADIUS, ECFP_BITS = 2, 2048


def _split_indices(raw_smiles, kind, seed, frac_train=0.8):
    """Pipeline-invariant split: computed once on the raw SMILES and reused for every pipeline."""
    if kind == "scaffold":
        return scaffold_split(raw_smiles, frac_train, seed)
    return random_split(len(raw_smiles), frac_train, seed)


def _eval(task, y_true, scores):
    if task == "reg":
        m = regression_metrics(y_true, scores)
        return {"rmse": m["rmse"], "mae": m["mae"], "spearman": m["spearman"], "r2": m["r2"]}
    m = classification_metrics(y_true, scores)
    return {"auroc": m["auroc"], "auprc": m["auprc"], "mcc": m["mcc"], "bacc": m["balanced_acc"],
            "ece": expected_calibration_error(y_true, scores)}


def run_cell(spec, pipeline_name, seed, models=MODEL_NAMES, feat="ecfp", split_kind="scaffold"):
    """Train+eval every model for one (dataset, pipeline, seed). Returns list of result rows."""
    df = D.load(spec)
    raw = df.smiles.tolist()
    y = df.y.to_numpy().astype(float)
    tr, te = _split_indices(raw, split_kind, seed)
    X = features(spec.name, raw, feature_spec_to_pipe(pipeline_name), feat)
    rows = []
    for mname in models:
        model = build_model(mname, spec.task, seed)
        model.fit(X[tr], y[tr])
        sc = predict_scores(model, X[te], spec.task)
        row = {"dataset": spec.name, "task": spec.task, "model": mname, "feat": feat,
               "pipeline": pipeline_name, "seed": seed, "split": split_kind,
               "n_train": len(tr), "n_test": len(te)}
        row.update(_eval(spec.task, y[te], sc))
        rows.append(row)
    return rows


def run_sorb_cell(spec, seed, models=("lgbm", "rf"), feat="ecfp", split_kind="scaffold",
                  train_aug=4, test_orbit=8, orbit_cfg: OrbitConfig | None = None):
    """SORB vs canonicalize-once baseline for one (dataset, seed). Returns rows with method in
    {baseline, sorb_trainaug, sorb_testorbit, sorb_full} to support the method ablation (RQ4)."""
    from models.sorb import SORB
    from preprocessing.cache import dataset_orbits
    df = D.load(spec)
    raw = df.smiles.tolist(); y = df.y.to_numpy().astype(float)
    tr, te = _split_indices(raw, split_kind, seed)
    cfg = orbit_cfg or OrbitConfig()
    orbits = dataset_orbits(spec.name, raw, cfg)  # cached; Dimorphite runs once per dataset
    tr_smi = [raw[i] for i in tr]; te_smi = [raw[i] for i in te]
    tr_orb = [orbits[i] for i in tr]; te_orb = [orbits[i] for i in te]
    rows = []
    for mname in models:
        variants = {
            "baseline": dict(train_aug=1, test_orbit=1),
            "sorb_trainaug": dict(train_aug=train_aug, test_orbit=1),
            "sorb_testorbit": dict(train_aug=1, test_orbit=test_orbit),
            "sorb_full": dict(train_aug=train_aug, test_orbit=test_orbit),
        }
        for vname, kw in variants.items():
            est = SORB(mname, spec.task, feat, seed=seed, orbit_cfg=cfg, **kw)
            est.fit(tr_smi, y[tr], orbits=tr_orb)
            sc = est.predict(te_smi, orbits=te_orb)
            row = {"dataset": spec.name, "task": spec.task, "model": mname, "method": vname,
                   "seed": seed, "split": split_kind, "train_aug": kw["train_aug"],
                   "test_orbit": kw["test_orbit"]}
            row.update(_eval(spec.task, y[te], sc))
            rows.append(row)
    return rows


def measure_prediction_invariance(spec, seed, model_name="lgbm", feat="ecfp",
                                  split_kind="scaffold", n_test=400, orbit_cfg=None):
    """RQ5(a): for a model trained under canonicalize-once, measure the spread of its predictions
    across each test molecule's standardization orbit. True property is orbit-invariant, so any
    spread is spurious standardization sensitivity. Returns per-molecule orbit-prediction SD stats.
    """
    from preprocessing.cache import dataset_orbits
    cfg = orbit_cfg or OrbitConfig()
    df = D.load(spec); raw = df.smiles.tolist(); y = df.y.to_numpy().astype(float)
    tr, te = _split_indices(raw, split_kind, seed)
    Xtr = features(spec.name, raw, feature_spec_to_pipe("P_rdkit"), feat)
    model = build_model(model_name, spec.task, seed)
    model.fit(Xtr[tr], y[tr])
    orbits = dataset_orbits(spec.name, raw, cfg)
    rng = np.random.default_rng(seed)
    sel = te if len(te) <= n_test else te[rng.choice(len(te), n_test, replace=False)]
    sds = []
    for i in sel:
        orb = orbits[i]
        if len(orb) < 2:
            sds.append(0.0); continue
        Xo = ecfp_matrix(orb, ECFP_RADIUS, ECFP_BITS) if feat == "ecfp" else rdkit_descriptors(orb)[0]
        pr = predict_scores(model, Xo, spec.task)
        sds.append(float(np.std(pr)))
    sds = np.array(sds)
    # express regression SD relative to the training-label SD for comparability across datasets
    scale = float(np.std(y[tr])) if spec.task == "reg" else 1.0
    return {"dataset": spec.name, "task": spec.task, "model": model_name, "seed": seed,
            "mean_orbit_pred_sd": float(sds.mean()), "median_orbit_pred_sd": float(np.median(sds)),
            "p90_orbit_pred_sd": float(np.percentile(sds, 90)), "label_sd": scale,
            "mean_rel_sd": float(sds.mean() / scale) if scale > 0 else float("nan"),
            "frac_orbit_gt1pct": float(np.mean(sds > 0.01 * scale)) if scale > 0 else float("nan")}
