#!/usr/bin/env python3
"""Reproduce the full study. Usage:
    python run_all.py --mode quick   # fast pipeline check (~minutes)
    python run_all.py --mode full    # complete experiments (~2 h wall on 4 cores)

Outputs under results/: raw_results.csv (sensitivity + ablation grid), sorb_results.csv,
invariance_results.csv, semisynthetic_results.csv, plus per-cell checkpoints in results/_cells/.
Re-running skips completed cells. Deterministic given seeds.
"""
from __future__ import annotations
import os, sys, json, argparse, time, hashlib
from multiprocessing import Pool
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(ROOT, "src"))

from utils.seeds import set_global_seed, detect_device
from data import datasets as D
from training import grid as G
from training.harness import run_cell
from training.invariance import run_method_cell
from training.shortcut import run_shortcut_cell
from preprocessing.cache import features, feature_spec_to_pipe, dataset_orbits
from preprocessing.standardize import OrbitConfig, PROTONATION_BACKEND

RESULTS = os.path.join(ROOT, "results")
CELLS = os.path.join(RESULTS, "_cells")
LOGS = os.path.join(RESULTS, "experiment_logs")
for d in (RESULTS, CELLS, LOGS):
    os.makedirs(d, exist_ok=True)

ORBIT_CFG = OrbitConfig(max_tautomers=3, max_members=6)


def _cell_path(tag: str) -> str:
    return os.path.join(CELLS, hashlib.md5(tag.encode()).hexdigest()[:20] + ".json")


def _log_fail(tag, e):
    with open(os.path.join(LOGS, "failures.log"), "a") as f:
        f.write(f"{tag}\t{type(e).__name__}\t{str(e)[:200]}\n")


def _cached(tag: str):
    p = _cell_path(tag)
    if os.path.exists(p):
        try:
            return json.load(open(p))
        except Exception:
            return None
    return None


def _save(tag: str, rows):
    json.dump(rows, open(_cell_path(tag), "w"))


# ---- precompute all feature + orbit caches (parallel across datasets, no write races) --------
def _precompute_one(args):
    spec, pipes = args
    try:
        df = D.load(spec); raw = df.smiles.tolist()
        for pname in pipes:
            for feat in ("ecfp", "desc"):
                features(spec.name, raw, feature_spec_to_pipe(pname), feat)
        dataset_orbits(spec.name, raw, ORBIT_CFG)
        return spec.name
    except Exception as e:
        return f"ERROR:{spec.name}:{type(e).__name__}:{str(e)[:120]}"


def precompute(cfg, nproc):
    pipes = sorted(set(cfg["sens_pipes"]) | set(cfg["abl_pipes"]) | {"P_rdkit"})
    jobs = [(s, pipes) for s in cfg["specs"]]
    t0 = time.time()
    with Pool(nproc) as pool:
        for i, name in enumerate(pool.imap_unordered(_precompute_one, jobs), 1):
            print(f"  precompute {i}/{len(jobs)}: {name} ({time.time()-t0:.0f}s)", flush=True)


# ---- grid workers ----------------------------------------------------------------------------
def _sens_job(args):
    spec, pname, seed, models, split = args
    tag = f"sens|{spec.name}|{pname}|{seed}|{','.join(models)}|{split}"
    c = _cached(tag)
    if c is not None:
        return c
    try:
        rows = run_cell(spec, pname, seed, models=models, feat="ecfp", split_kind=split)
        rows += [{**r, "feat": "desc"} for r in run_cell(spec, pname, seed, models=models, feat="desc", split_kind=split)]
    except Exception as e:
        _log_fail(tag, e); return []
    _save(tag, rows)
    return rows


def _abl_job(args):
    spec, pname, seed, models, split = args
    tag = f"abl|{spec.name}|{pname}|{seed}|{','.join(models)}|{split}"
    c = _cached(tag)
    if c is not None:
        return c
    try:
        rows = run_cell(spec, pname, seed, models=models, feat="ecfp", split_kind=split)
        for r in rows:
            r["ablation_pipe"] = pname
    except Exception as e:
        _log_fail(tag, e); return []
    _save(tag, rows)
    return rows


def _shortcut_job(args):
    spec, seed, model, feat, split = args
    tag = f"shortcut|{spec.name}|{seed}|{model}|{feat}"
    c = _cached(tag)
    if c is not None:
        return c
    try:
        rows = [run_shortcut_cell(spec, seed, model_name=model, feat=feat, orbit_cfg=ORBIT_CFG)]
    except Exception as e:
        _log_fail(tag, e); return []
    _save(tag, rows)
    return rows


def _method_job(args):
    spec, seed, model, feat, split = args
    tag = f"method|{spec.name}|{seed}|{model}|{feat}|{split}"
    c = _cached(tag)
    if c is not None:
        return c
    try:
        rows = run_method_cell(spec, seed, model_name=model, feat=feat, split_kind=split, orbit_cfg=ORBIT_CFG)
    except Exception as e:
        _log_fail(tag, e); return []
    _save(tag, rows)
    return rows


def _run_pool(fn, jobs, nproc, label):
    out = []
    t0 = time.time()
    with Pool(nproc) as pool:
        for i, rows in enumerate(pool.imap_unordered(fn, jobs), 1):
            out.extend(rows)
            if i % 20 == 0 or i == len(jobs):
                print(f"  {label} {i}/{len(jobs)} ({time.time()-t0:.0f}s)", flush=True)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["quick", "full"], default="quick")
    ap.add_argument("--nproc", type=int, default=max(1, min(4, os.cpu_count() or 1)))
    ap.add_argument("--skip-precompute", action="store_true")
    args = ap.parse_args()
    set_global_seed(0)
    cfg = G.grid(args.mode)
    meta = {"mode": args.mode, "nproc": args.nproc, "device": detect_device(),
            "protonation_backend": PROTONATION_BACKEND, "n_datasets": len(cfg["specs"]),
            "orbit_cfg": vars(ORBIT_CFG)}
    json.dump(meta, open(os.path.join(LOGS, f"run_meta_{args.mode}.json"), "w"), indent=1)
    print("RUN", meta, flush=True)

    if not args.skip_precompute:
        print("[precompute caches]", flush=True)
        precompute(cfg, args.nproc)

    specs, models, split = cfg["specs"], cfg["models"], cfg["split"]
    print("[E1/E2 sensitivity grid]", flush=True)
    jobs = [(s, p, seed, models, split) for s in specs for p in cfg["sens_pipes"] for seed in cfg["seeds"]]
    sens = _run_pool(_sens_job, jobs, args.nproc, "sens")
    pd.DataFrame(sens).to_csv(os.path.join(RESULTS, "raw_results.csv"), index=False)

    print("[E3 operation ablation grid]", flush=True)
    jobs = [(s, p, seed, models, split) for s in specs for p in cfg["abl_pipes"] for seed in cfg["seeds"]]
    abl = _run_pool(_abl_job, jobs, args.nproc, "abl")
    pd.DataFrame(abl).to_csv(os.path.join(RESULTS, "ablation_results.csv"), index=False)

    # method/shortcut: cap to datasets <= METHOD_MAX_N molecules (orbit augmentation on very large
    # sets is memory-heavy; the largest sets remain in the sensitivity/ablation analyses).
    METHOD_MAX_N = 6000
    def _nmol(sp):
        try:
            return len(D.load(sp))
        except Exception:
            return 10**9
    method_specs = [s for s in specs if _nmol(s) <= METHOD_MAX_N]
    print(f"[E4/E5 method: SORB + invariance] on {len(method_specs)}/{len(specs)} datasets (<= {METHOD_MAX_N} mol)", flush=True)
    jobs = [(s, seed, m, f, split) for s in method_specs for seed in cfg["method_seeds"]
            for m in cfg["method_models"] for f in cfg["method_feats"]]
    meth = _run_pool(_method_job, jobs, args.nproc, "method")
    pd.DataFrame(meth).to_csv(os.path.join(RESULTS, "method_results.csv"), index=False)

    print("[E5b semi-synthetic shortcut]", flush=True)
    jobs = [(s, seed, m, "ecfp", split) for s in method_specs for seed in cfg["shortcut_seeds"]
            for m in cfg["shortcut_models"]]
    sc = _run_pool(_shortcut_job, jobs, args.nproc, "shortcut")
    pd.DataFrame(sc).to_csv(os.path.join(RESULTS, "shortcut_results.csv"), index=False)

    print("DONE. Results in", RESULTS, flush=True)


if __name__ == "__main__":
    main()
