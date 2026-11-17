#!/usr/bin/env python3
"""Strata-count sensitivity for AD-conditional conformal (SI robustness).
For each regression dataset (lgbm + ECFP, scaffold split, seed 0): calibrate Mondrian-AD and smoothed-AD
with m in {3,5,10} strata, but evaluate conditional coverage on a FIXED 5-quintile reference grid so the
worst-stratum coverage is comparable across m. Also records marginal conformal as the reference.
Writes results/conformal_strata_results.csv. Checkpointed."""
from __future__ import annotations
import os, sys, json, hashlib
from multiprocessing import Pool
import numpy as np, pandas as pd
ROOT = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(ROOT, "src"))
from utils.seeds import set_global_seed
from data import datasets as D
from data.splits import three_way
from preprocessing.cache import features, feature_spec_to_pipe
from models.zoo import build_model
from conformal.methods import ad_distance, marginal, mondrian, smoothed_mondrian, group_coverage, conformal_quantile
RES = os.path.join(ROOT, "results"); CK = os.path.join(RES, "_stratacells"); os.makedirs(CK, exist_ok=True)
DROP = {"molnet:lipophilicity"}
M_LIST = [3, 5, 10]
EVAL_BINS = 5
ALPHA = 0.1


def _bins(dist_cal, dist_te, n):
    edges = np.quantile(dist_cal, np.linspace(0, 1, n + 1)); edges[0], edges[-1] = -np.inf, np.inf
    gcal = np.clip(np.digitize(dist_cal, edges[1:-1]), 0, n - 1)
    gte = np.clip(np.digitize(dist_te, edges[1:-1]), 0, n - 1)
    return gcal, gte


def run(spec):
    df = D.load(spec); raw = df.smiles.tolist(); y = df.y.to_numpy().astype(float)
    tr, cal, te = three_way(raw, "scaffold", 0.7, 0.15, 0)
    X = features(spec.name, raw, feature_spec_to_pipe("P_rdkit"), "ecfp")
    m = build_model("lgbm", spec.task, 0); m.fit(X[tr], y[tr])
    res_cal = y[cal] - m.predict(X[cal]); pred_te = m.predict(X[te]); yte = y[te]
    tr_smi = [raw[i] for i in tr]
    d_cal = ad_distance([raw[i] for i in cal], tr_smi); d_te = ad_distance([raw[i] for i in te], tr_smi)
    # fixed 5-bin reference evaluation grid
    _, gte_eval = _bins(d_cal, d_te, EVAL_BINS)

    def worst(h):
        gc = group_coverage(yte, pred_te, h, gte_eval)
        return (min(gc.values()) if gc else float("nan"),
                gc.get(EVAL_BINS - 1, float("nan")),
                float(np.mean((yte >= pred_te - h) & (yte <= pred_te + h))))
    rows = []
    hm = marginal(res_cal, len(te), ALPHA)
    w, f, mc = worst(hm)
    rows.append({"dataset": spec.name, "m": 0, "method": "marginal", "worst_bin_cov": w, "far_cov": f, "marginal_cov": mc})
    for M in M_LIST:
        gcal, gte = _bins(d_cal, d_te, M)
        hmo = mondrian(res_cal, gcal, gte, ALPHA)
        w, f, mc = worst(hmo)
        rows.append({"dataset": spec.name, "m": M, "method": "mondrian_ad", "worst_bin_cov": w, "far_cov": f, "marginal_cov": mc})
        hsm = smoothed_mondrian(res_cal, d_cal, d_te, ALPHA, n_bins=M, shrink=0.3)
        w, f, mc = worst(hsm)
        rows.append({"dataset": spec.name, "m": M, "method": "smoothed_ad", "worst_bin_cov": w, "far_cov": f, "marginal_cov": mc})
    return rows


def _job(spec):
    cp = os.path.join(CK, hashlib.md5(spec.name.encode()).hexdigest()[:16] + ".json")
    if os.path.exists(cp):
        try: return json.load(open(cp))
        except Exception: pass
    try:
        r = run(spec)
    except Exception as e:
        r = [{"dataset": spec.name, "error": f"{type(e).__name__}:{str(e)[:150]}"}]
    json.dump(r, open(cp, "w")); return r


def main():
    import argparse; ap = argparse.ArgumentParser(); ap.add_argument("--nproc", type=int, default=3); a = ap.parse_args()
    set_global_seed(0)
    specs = [s for s in D.regression_registry() if s.name not in DROP]
    print(f"strata sensitivity: {len(specs)} datasets x m in {M_LIST}", flush=True)
    out = []
    with Pool(a.nproc) as pool:
        for i, rows in enumerate(pool.imap_unordered(_job, specs), 1):
            out.extend([r for r in rows if "error" not in r])
            if i % 10 == 0 or i == len(specs): print(f"  {i}/{len(specs)}", flush=True)
    pd.DataFrame(out).to_csv(os.path.join(RES, "conformal_strata_results.csv"), index=False)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
