"""Phase 8: data leakage / integrity audit of every dataset used.

Checks, per dataset and split (scaffold, seed 0):
  - exact duplicate molecules (by canonical identity) and label disagreement among them
  - train/test near-duplicate rate (max Tanimoto of each test molecule to the train set >= 0.95)
  - standardization-collapse: distinct raw SMILES that map to the same canonical identity
  - split integrity: no canonical-identity leakage across train/test (by construction of our splitter)
Writes results/leakage_audit.csv and a short markdown summary.
"""
from __future__ import annotations
import os, sys, json
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
import pandas as pd
from data import datasets as D
from data.splits import scaffold_split
from preprocessing.standardize import canonical_identity
from preprocessing.featurize import ecfp_bitset, tanimoto
from rdkit import DataStructs

OUT = os.path.join(os.path.dirname(__file__), "..", "results")


def nn_tanimoto_rate(train_smiles, test_smiles, thr=0.95, cap=1500):
    tr_fps = [ecfp_bitset(s) for s in train_smiles]
    tr_fps = [f for f in tr_fps if f is not None]
    rng = np.random.default_rng(0)
    idx = range(len(test_smiles)) if len(test_smiles) <= cap else rng.choice(len(test_smiles), cap, replace=False)
    hits = 0; n = 0
    for i in idx:
        f = ecfp_bitset(test_smiles[i])
        if f is None or not tr_fps:
            continue
        sims = DataStructs.BulkTanimotoSimilarity(f, tr_fps)
        hits += int(max(sims) >= thr); n += 1
    return hits / n if n else float("nan")


def audit(spec):
    df = D.load(spec); raw = df.smiles.tolist(); y = df.y.to_numpy()
    ids = [canonical_identity(s) for s in raw]
    valid = [i for i, c in enumerate(ids) if c]
    id_arr = np.array([ids[i] for i in valid])
    # exact duplicates by canonical identity
    uniq, counts = np.unique(id_arr, return_counts=True)
    n_dup_groups = int((counts > 1).sum())
    n_dup_mols = int(counts[counts > 1].sum())
    # label disagreement among duplicates (regression: SD; clf: mixed labels)
    disagree = 0
    for u in uniq[counts > 1]:
        yy = y[[valid[k] for k in range(len(valid)) if id_arr[k] == u]]
        if spec.task == "clf":
            disagree += int(len(set(yy.tolist())) > 1)
        else:
            disagree += int(np.std(yy) > 0.1 * (np.std(y) + 1e-9))
    # standardization-collapse: distinct raw -> same identity
    n_raw_unique = len(set(raw)); n_id_unique = len(set(id_arr))
    collapse = n_raw_unique - n_id_unique
    # split integrity + train/test near-dup
    tr, te = scaffold_split(raw, 0.8, 0)
    tr_ids = {ids[i] for i in tr if ids[i]}; te_ids = {ids[i] for i in te if ids[i]}
    id_leak = len(tr_ids & te_ids)
    nn_rate = nn_tanimoto_rate([raw[i] for i in tr], [raw[i] for i in te])
    return {"dataset": spec.name, "task": spec.task, "n": len(raw), "n_valid": len(valid),
            "dup_groups": n_dup_groups, "dup_mols": n_dup_mols, "dup_label_disagree": disagree,
            "std_collapse_pairs": collapse, "scaffold_id_leak": id_leak,
            "train_test_nn_ge95pct": round(nn_rate, 4)}


def main():
    specs = D.regression_registry() + [D.DatasetSpec(f"tdc:{t}", "tdc", t, "clf") for t in D.TDC_CLF]
    rows = [audit(s) for s in specs]
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, "leakage_audit.csv"), index=False)
    with open(os.path.join(OUT, "leakage_audit_summary.md"), "w") as f:
        f.write("# Data leakage / integrity audit (scaffold split, seed 0)\n\n")
        f.write(df.to_markdown(index=False))
        f.write("\n\n**Key:** dup_* = exact duplicates by standardization-invariant canonical identity; "
                "std_collapse_pairs = distinct raw SMILES collapsing to one molecule after standardization "
                "(these would silently merge/leak if not handled); scaffold_id_leak = canonical-identity "
                "overlap across the train/test split (0 by construction of our identity-aware splitter); "
                "train_test_nn_ge95pct = fraction of test molecules with a >=0.95-Tanimoto neighbour in train.\n")
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
