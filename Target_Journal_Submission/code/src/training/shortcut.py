"""Semi-synthetic standardization-shortcut experiment (controlled ground truth).

We make the *standardization state* artificially predictive of the label: a molecule with >=2 orbit
members is presented as member A if its label is high and member B if low. A standard model can then
exploit which standardization was used as a spurious feature. We measure:
  - honest_ref        : baseline trained & tested on canonical form (no shortcut) -> reference acc
  - shortcut_biased   : model trained on shortcut data, tested WITH the shortcut (inflated)
  - shortcut_cf       : same model tested on COUNTERFACTUAL (canonical) form (shortcut removed)
  - sorb_cf           : SORB train-aug model (orbit, shared label) tested counterfactually
A model that exploited the shortcut shows a large drop shortcut_biased -> shortcut_cf; an invariant
model (SORB) does not. The ground truth is that the property is standardization-invariant, so the
honest reference is the correct target.
"""
from __future__ import annotations
import numpy as np
from data import datasets as D
from data.splits import scaffold_split
from preprocessing.cache import features, feature_spec_to_pipe, dataset_orbits
from preprocessing.featurize import ecfp_matrix, rdkit_descriptors
from preprocessing.standardize import OrbitConfig
from models.zoo import build_model, predict_scores
from evaluation.metrics import regression_metrics, classification_metrics


def _feat(smiles, feat):
    return ecfp_matrix(smiles, 2, 2048) if feat == "ecfp" else rdkit_descriptors(smiles)[0]


def _primary(task, y, s):
    return regression_metrics(y, s)["rmse"] if task == "reg" else classification_metrics(y, s)["auroc"]


def _present(orbits, idx, hi_mask, memberA=0, memberB=1):
    """Choose orbit member per molecule: A if hi, else B (clamped to available members)."""
    out = []
    for k, i in enumerate(idx):
        orb = orbits[i]
        if len(orb) >= 2:
            out.append(orb[memberA] if hi_mask[k] else orb[min(memberB, len(orb) - 1)])
        else:
            out.append(orb[0] if orb else "")
    return out


def run_shortcut_cell(spec, seed, model_name="lgbm", feat="ecfp", train_aug=4,
                      orbit_cfg: OrbitConfig | None = None):
    cfg = orbit_cfg or OrbitConfig()
    df = D.load(spec); raw = df.smiles.tolist(); y = df.y.to_numpy().astype(float)
    tr, te = scaffold_split(raw, 0.8, seed)
    orbits = dataset_orbits(spec.name, raw, cfg)
    thr = np.median(y[tr])
    hi_tr = y[tr] > thr if spec.task == "reg" else y[tr] > 0.5
    hi_te = y[te] > thr if spec.task == "reg" else y[te] > 0.5

    # canonical (member 0) featurization
    Xc = features(spec.name, raw, feature_spec_to_pipe("P_rdkit"), feat)
    # honest reference
    m0 = build_model(model_name, spec.task, seed); m0.fit(Xc[tr], y[tr])
    honest_ref = _primary(spec.task, y[te], predict_scores(m0, Xc[te], spec.task))

    # shortcut-biased training set
    tr_short = _present(orbits, tr, hi_tr)
    Xtr_short = _feat(tr_short, feat)
    ms = build_model(model_name, spec.task, seed); ms.fit(Xtr_short, y[tr])
    # tested WITH shortcut (members chosen by true test label) and COUNTERFACTUAL (canonical)
    te_short = _present(orbits, te, hi_te)
    biased = _primary(spec.task, y[te], predict_scores(ms, _feat(te_short, feat), spec.task))
    cf = _primary(spec.task, y[te], predict_scores(ms, Xc[te], spec.task))

    # SORB train-aug on the SHORTCUT training molecules (augment over orbit, shared label) -> cannot
    # use standardization state as signal; test counterfactually (canonical)
    from models.sorb import SORB
    est = SORB(model_name, spec.task, feat, seed=seed, orbit_cfg=cfg, train_aug=train_aug, test_orbit=1)
    tr_orb = [orbits[i] for i in tr]
    est.fit([raw[i] for i in tr], y[tr], orbits=tr_orb)
    sorb_cf = _primary(spec.task, y[te], predict_scores(est.model, Xc[te], spec.task))

    lower_better = spec.task == "reg"
    def drop(a, b):  # performance drop from a to b (positive = got worse)
        return (b - a) if lower_better else (a - b)
    return {"dataset": spec.name, "task": spec.task, "model": model_name, "feat": feat, "seed": seed,
            "honest_ref": honest_ref, "shortcut_biased": biased, "shortcut_cf": cf, "sorb_cf": sorb_cf,
            "baseline_shortcut_drop": float(drop(biased, cf)),
            "baseline_cf_vs_honest": float(drop(honest_ref, cf)),
            "sorb_cf_vs_honest": float(drop(honest_ref, sorb_cf))}
