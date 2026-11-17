"""Revision experiments addressing round-1 review (run on the 49 datasets <=6000 molecules).

E1 cross-training-pipeline disagreement (the key justification for SORB over 'just canonicalize'):
   train two models on two different but equally-valid pipelines (P_rdkit vs P_taut); on a common
   canonical test set, measure median |pred_A - pred_B| / label_sd. Baseline models (pipeline-specific
   training) vs SORB-trainaug models (orbit augmentation). If SORB shrinks the disagreement, it buys
   cross-lab reproducibility that canonicalize-at-inference cannot.
E2 canonicalize-at-inference control: baseline model; present each test molecule's orbit members but
   re-canonicalize each to P_rdkit before featurizing -> orbit SD should collapse to ~0 (the trivial
   fix for input-presentation variance; included so we can state it explicitly).
E3 random-split replication of the per-prediction orbit invariance (baseline), to support robustness.
Outputs results/revision_results.csv.
"""
from __future__ import annotations
import os, sys, json, argparse
from multiprocessing import Pool
import numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from data import datasets as D
from data.splits import scaffold_split, random_split
from preprocessing.cache import features, feature_spec_to_pipe, dataset_orbits
from preprocessing.featurize import ecfp_matrix, rdkit_descriptors
from preprocessing.standardize import NAMED_PIPELINES, OrbitConfig
from models.zoo import build_model, predict_scores
from models.sorb import SORB

ROOT = os.path.join(os.path.dirname(__file__), "..")
RES = os.path.join(ROOT, "results")
CFG = OrbitConfig(max_tautomers=3, max_members=6)
CANON = feature_spec_to_pipe("P_rdkit")


def _feat(smiles, feat):
    return ecfp_matrix(smiles, 2, 2048) if feat == "ecfp" else rdkit_descriptors(smiles)[0]


def _label_scale(y, task):
    return float(np.std(y)) if task == "reg" else float(np.sqrt(np.mean(y) * (1 - np.mean(y))) + 1e-9)


def run(spec, seed, model_name, feat):
    df = D.load(spec); raw = df.smiles.tolist(); y = df.y.to_numpy().astype(float)
    out = {}
    for split_kind in ("scaffold", "random"):
        tr, te = (scaffold_split(raw, 0.8, seed) if split_kind == "scaffold" else random_split(len(raw), 0.8, seed))
        scale = _label_scale(y[tr], spec.task)
        orbits = dataset_orbits(spec.name, raw, CFG)
        # common canonical test features
        Xc = features(spec.name, raw, CANON, feat)
        # --- E1 cross-pipeline disagreement (scaffold only to save compute) ---
        if split_kind == "scaffold":
            Xa = features(spec.name, raw, feature_spec_to_pipe("P_rdkit"), feat)
            Xb = features(spec.name, raw, feature_spec_to_pipe("P_taut"), feat)
            ma = build_model(model_name, spec.task, seed); ma.fit(Xa[tr], y[tr])
            mb = build_model(model_name, spec.task, seed); mb.fit(Xb[tr], y[tr])
            pa = predict_scores(ma, Xc[te], spec.task); pb = predict_scores(mb, Xc[te], spec.task)
            out["xpipe_baseline_med"] = float(np.median(np.abs(pa - pb)) / scale)
            # SORB-trainaug trained "from" each pipeline presentation but augmenting over the orbit
            tr_orb = [orbits[i] for i in tr]
            sa = SORB(model_name, spec.task, feat, seed=seed, orbit_cfg=CFG, train_aug=4, test_orbit=1)
            sa.fit([raw[i] for i in tr], y[tr], orbits=tr_orb)
            # a second SORB with a different augmentation seed (models that two labs would independently fit)
            sb = SORB(model_name, spec.task, feat, seed=seed + 101, orbit_cfg=CFG, train_aug=4, test_orbit=1)
            sb.fit([raw[i] for i in tr], y[tr], orbits=tr_orb)
            psa = predict_scores(sa.model, Xc[te], spec.task); psb = predict_scores(sb.model, Xc[te], spec.task)
            out["xpipe_sorb_med"] = float(np.median(np.abs(psa - psb)) / scale)
            # --- E2 the trivial fix: a canonical-TAUTOMER pipeline (P_taut) applied consistently at
            # train and inference. Baseline P_rdkit does NOT canonicalize tautomers (the dominant axis),
            # so its orbit SD persists; P_taut projects each presentation to one tautomer, collapsing it.
            ptaut = feature_spec_to_pipe("P_taut")
            base = build_model(model_name, spec.task, seed); base.fit(Xc[tr], y[tr])          # P_rdkit-trained
            baset = build_model(model_name, spec.task, seed); baset.fit(Xb[tr], y[tr])        # P_taut-trained
            rng = np.random.default_rng(seed); sel = te if len(te) <= 200 else te[rng.choice(len(te), 200, replace=False)]
            sds_raw, sds_canon = [], []
            for i in sel:
                orb = orbits[i]
                if len(orb) < 2:
                    sds_raw.append(0.0); sds_canon.append(0.0); continue
                sds_raw.append(float(np.std(predict_scores(base, _feat(orb, feat), spec.task))))
                canon_orb = [ptaut.apply(s) or s for s in orb]   # project each presentation to canonical tautomer
                sds_canon.append(float(np.std(predict_scores(baset, _feat(canon_orb, feat), spec.task))))
            out["orbit_sd_raw"] = float(np.mean(sds_raw) / scale)
            out["orbit_sd_canon_taut"] = float(np.mean(sds_canon) / scale)
        else:
            # --- E3 random-split baseline orbit SD (properly scaled) ---
            base = build_model(model_name, spec.task, seed); base.fit(Xc[tr], y[tr])
            rng = np.random.default_rng(seed); sel = te if len(te) <= 200 else te[rng.choice(len(te), 200, replace=False)]
            sds = []
            for i in sel:
                orb = orbits[i]
                if len(orb) < 2: sds.append(0.0); continue
                sds.append(float(np.std(predict_scores(base, _feat(orb, feat), spec.task))))
            out["orbit_sd_random"] = float(np.mean(sds) / scale)
    out.update({"dataset": spec.name, "task": spec.task, "model": model_name, "feat": feat, "seed": seed, "label_scale": scale})
    return out


import hashlib, json as _json
CKPT = os.path.join(RES, "_revcells"); os.makedirs(CKPT, exist_ok=True)
def _job(a):
    tag = f"{a[0].name}|{a[1]}|{a[2]}|{a[3]}"
    cp = os.path.join(CKPT, hashlib.md5(tag.encode()).hexdigest()[:16] + ".json")
    if os.path.exists(cp):
        try: return _json.load(open(cp))
        except Exception: pass
    try:
        r = run(*a)
    except Exception as e:
        r = {"dataset": a[0].name, "error": f"{type(e).__name__}:{e}"}
    _json.dump(r, open(cp, "w")); return r


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--nproc", type=int, default=2); args = ap.parse_args()
    specs = [s for s in D.regression_registry() + [D.DatasetSpec(f"tdc:{t}", "tdc", t, "clf") for t in D.TDC_CLF]
             if len(D.load(s)) <= 6000 and s.name != "molnet:lipophilicity"]  # drop duplicate logD benchmark
    jobs = [(s, 0, m, f) for s in specs for m in ("lgbm", "linear") for f in ("ecfp", "desc")]  # seed 0; dataset is the unit
    print(f"revision experiments: {len(jobs)} jobs on {len(specs)} datasets", flush=True)
    rows = []
    with Pool(args.nproc) as pool:
        for i, r in enumerate(pool.imap_unordered(_job, jobs), 1):
            rows.append(r)
            if i % 20 == 0: print(f"  {i}/{len(jobs)}", flush=True)
    pd.DataFrame(rows).to_csv(os.path.join(RES, "revision_results.csv"), index=False)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
