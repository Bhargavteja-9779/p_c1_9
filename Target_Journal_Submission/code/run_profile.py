#!/usr/bin/env python3
"""Per-stratum coverage profile across ALL five applicability-domain strata (for the headline
coverage-vs-AD figure). For each regression dataset (lgbm + ECFP, scaffold split, seed 0), compute the
coverage of marginal / normalized / mondrian_ad / smoothed_ad in each of the 5 AD strata (quintiles of
calibration distance), then aggregate across datasets. Writes results/conformal_profile_results.csv.
Reuses cached features; fast. Checkpointed."""
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
from conformal.methods import ad_distance, marginal, normalized, distance_normalized, mondrian, smoothed_mondrian
import lightgbm as lgb
RES = os.path.join(ROOT, "results"); CK = os.path.join(RES, "_profilecells"); os.makedirs(CK, exist_ok=True)
DROP = {"molnet:lipophilicity"}
NB = 5; ALPHA = 0.1


def _bins(dc, dt, n=NB):
    e = np.quantile(dc, np.linspace(0, 1, n + 1)); e[0], e[-1] = -np.inf, np.inf
    return (np.clip(np.digitize(dc, e[1:-1]), 0, n - 1), np.clip(np.digitize(dt, e[1:-1]), 0, n - 1))


def run(spec):
    df = D.load(spec); raw = df.smiles.tolist(); y = df.y.to_numpy().astype(float)
    tr, cal, te = three_way(raw, "scaffold", 0.7, 0.15, 0)
    X = features(spec.name, raw, feature_spec_to_pipe("P_rdkit"), "ecfp")
    m = build_model("lgbm", spec.task, 0); m.fit(X[tr], y[tr])
    res_cal = y[cal] - m.predict(X[cal]); pred_te = m.predict(X[te]); yte = y[te]
    tr_smi = [raw[i] for i in tr]
    d_cal = ad_distance([raw[i] for i in cal], tr_smi); d_te = ad_distance([raw[i] for i in te], tr_smi)
    gcal, gte = _bins(d_cal, d_te)
    # out-of-fold sigma for normalized (same recipe as experiment.py)
    tr_arr = np.asarray(tr); Xtr = X[tr_arr]; ytr = y[tr_arr]; oof = np.empty(len(tr_arr))
    _p = np.random.default_rng(1000).permutation(len(tr_arr)); fs = np.array_split(_p, 2)
    for k in range(2):
        v = fs[k]; t = np.concatenate([fs[j] for j in range(2) if j != k])
        if len(t) < 2 or len(v) == 0: oof[v] = np.abs(res_cal).mean() if len(res_cal) else 0.0; continue
        bk = build_model("lgbm", spec.task, 0); bk.fit(Xtr[t], ytr[t]); oof[v] = np.abs(ytr[v] - bk.predict(Xtr[v]))
    sm = lgb.LGBMRegressor(n_estimators=200, num_leaves=15, min_child_samples=20, random_state=0, n_jobs=1, verbosity=-1)
    sm.fit(Xtr, oof); sig_cal = np.maximum(sm.predict(X[cal]), 1e-6); sig_te = np.maximum(sm.predict(X[te]), 1e-6)
    halves = {"marginal": marginal(res_cal, len(te), ALPHA),
              "normalized": normalized(res_cal, sig_cal, sig_te, ALPHA),
              "knn_normalized": distance_normalized(res_cal, d_cal, d_te, ALPHA),
              "mondrian_ad": mondrian(res_cal, gcal, gte, ALPHA),
              "smoothed_ad": smoothed_mondrian(res_cal, d_cal, d_te, ALPHA, n_bins=NB, shrink=0.3)}
    scale = float(np.std(y[tr])) + 1e-9
    covered = {k: ((yte >= pred_te - h) & (yte <= pred_te + h)).astype(float) for k, h in halves.items()}
    prof = []
    for b in range(NB):
        mask = gte == b
        if not mask.any():
            continue
        for k in halves:
            prof.append({"dataset": spec.name, "stratum": b, "method": k,
                         "coverage": float(covered[k][mask].mean()), "n": int(mask.sum())})
    # per-compound selective-prediction records
    comp = []
    for k in ("marginal", "knn_normalized", "mondrian_ad", "smoothed_ad"):
        h = halves[k]
        for i in range(len(te)):
            comp.append({"dataset": spec.name, "method": k, "covered": int(covered[k][i]),
                         "rel_width": float(2 * h[i] / scale), "stratum": int(gte[i]),
                         "far": int(gte[i] == NB - 1)})
    # mechanism: per-compound AD distance, true abs residual, and the LEARNED difficulty sigma(x).
    # If the learned sigma does NOT track AD distance while the true error does, that explains why
    # learned-difficulty normalized conformal fails on novel chemotypes.
    abs_res_te = np.abs(yte - pred_te)
    mech = [{"dataset": spec.name, "ad_dist": float(d_te[i]), "abs_resid": float(abs_res_te[i]),
             "sigma_learned": float(sig_te[i])} for i in range(len(te))]
    return {"profile": prof, "compound": comp, "mech": mech}


def _job(spec):
    cp = os.path.join(CK, hashlib.md5(spec.name.encode()).hexdigest()[:16] + ".json")
    if os.path.exists(cp):
        try: return json.load(open(cp))
        except Exception: pass
    try: r = run(spec)
    except Exception as e: r = {"profile": [], "compound": [], "error": f"{type(e).__name__}:{str(e)[:120]}"}
    json.dump(r, open(cp, "w")); return r


def main():
    import argparse; ap = argparse.ArgumentParser(); ap.add_argument("--nproc", type=int, default=3); a = ap.parse_args()
    set_global_seed(0)
    specs = [s for s in D.regression_registry() if s.name not in DROP]
    print(f"profile: {len(specs)} datasets x {NB} strata", flush=True)
    prof, comp, mech = [], [], []
    with Pool(a.nproc) as pool:
        for i, r in enumerate(pool.imap_unordered(_job, specs), 1):
            if isinstance(r, dict):
                prof.extend(r.get("profile", [])); comp.extend(r.get("compound", [])); mech.extend(r.get("mech", []))
            if i % 10 == 0 or i == len(specs): print(f"  {i}/{len(specs)}", flush=True)
    pd.DataFrame(prof).to_csv(os.path.join(RES, "conformal_profile_results.csv"), index=False)
    pd.DataFrame(comp).to_csv(os.path.join(RES, "conformal_selective_results.csv"), index=False)
    pd.DataFrame(mech).to_csv(os.path.join(RES, "conformal_mechanism_results.csv"), index=False)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
