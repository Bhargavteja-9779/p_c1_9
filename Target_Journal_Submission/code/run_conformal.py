#!/usr/bin/env python3
"""Full applicability-domain-conditional conformal study. Checkpointed & resumable.
Usage: python run_conformal.py --mode {quick,full} --nproc N [--splits scaffold,random] [--alphas 0.1]"""
from __future__ import annotations
import os, sys, json, argparse, hashlib, time
from multiprocessing import Pool
import numpy as np, pandas as pd
ROOT = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(ROOT, "src"))
from utils.seeds import set_global_seed
from data import datasets as D
from conformal.experiment import run_cell
RES = os.path.join(ROOT, "results"); CK = os.path.join(RES, "_confcells"); LOG = os.path.join(RES, "experiment_logs")
for d in (RES, CK, LOG): os.makedirs(d, exist_ok=True)
DROP = {"molnet:lipophilicity"}

def specs_reg():
    return [s for s in D.regression_registry() if s.name not in DROP]

def _job(a):
    spec, seed, model, feat, split, alpha = a
    tag = f"conf|{spec.name}|{seed}|{model}|{feat}|{split}|{alpha}"
    cp = os.path.join(CK, hashlib.md5(tag.encode()).hexdigest()[:18] + ".json")
    if os.path.exists(cp):
        try: return json.load(open(cp))
        except Exception: pass
    try:
        rows = run_cell(spec, seed, model_name=model, feat=feat, split_kind=split, alpha=alpha)
    except Exception as e:
        with open(os.path.join(LOG, "conf_failures.log"), "a") as f: f.write(f"{tag}\t{type(e).__name__}\t{str(e)[:200]}\n")
        rows = []
    json.dump(rows, open(cp, "w")); return rows

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--mode", default="quick"); ap.add_argument("--nproc", type=int, default=3)
    ap.add_argument("--splits", default="scaffold,random"); ap.add_argument("--alphas", default="0.1")
    ap.add_argument("--models", default=""); ap.add_argument("--feats", default=""); ap.add_argument("--seeds", default="")
    ap.add_argument("--out", default="conformal_results.csv")
    a = ap.parse_args(); set_global_seed(0)
    splits = a.splits.split(","); alphas = [float(x) for x in a.alphas.split(",")]
    sp = specs_reg()
    if a.mode == "quick":
        sp = [s for s in sp if s.name in ("mace:CHEMBL204_Ki", "molnet:esol", "tdc:caco2_wang", "tdc:ppbr_az")]
        models, feats, seeds = ["lgbm", "linear"], ["ecfp"], [0, 1]
    else:
        models, feats, seeds = ["lgbm", "rf", "linear"], ["ecfp", "desc"], [0, 1, 2]
    # optional overrides for robustness grids (e.g. alpha sweeps on a representative subset)
    if a.models: models = a.models.split(",")
    if a.feats: feats = a.feats.split(",")
    if a.seeds: seeds = [int(x) for x in a.seeds.split(",")]
    jobs = [(s, sd, m, f, sp_, al) for s in sp for sd in seeds for m in models for f in feats for sp_ in splits for al in alphas]
    print(f"conformal study: {len(jobs)} cells on {len(sp)} datasets", flush=True)
    out = []; t0 = time.time()
    with Pool(a.nproc) as pool:
        for i, rows in enumerate(pool.imap_unordered(_job, jobs), 1):
            out.extend(rows)
            if i % 50 == 0 or i == len(jobs): print(f"  {i}/{len(jobs)} ({time.time()-t0:.0f}s)", flush=True)
    pd.DataFrame(out).to_csv(os.path.join(RES, a.out), index=False)
    print("DONE", flush=True)

if __name__ == "__main__":
    main()
