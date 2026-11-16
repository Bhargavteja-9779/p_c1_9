"""Reproducible data splits for molecular datasets: random, scaffold (Bemis-Murcko),
and repeated resampling. All splits are deterministic given a seed.
"""
from __future__ import annotations
from collections import defaultdict
from typing import Iterator
import numpy as np
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

def murcko_scaffold(smiles: str, include_chirality: bool = False) -> str:
    """Bemis-Murcko scaffold SMILES ('' if unparseable)."""
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return ""
    try:
        return MurckoScaffold.MurckoScaffoldSmiles(mol=m, includeChirality=include_chirality)
    except Exception:
        return ""

def random_split(n: int, frac_train: float = 0.8, seed: int = 0) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    idx = rng.permutation(n)
    k = int(round(frac_train * n))
    return np.sort(idx[:k]), np.sort(idx[k:])

def scaffold_split(smiles: list[str], frac_train: float = 0.8, seed: int = 0) -> tuple[np.ndarray, np.ndarray]:
    """Group molecules by Murcko scaffold; assign whole scaffold groups to train/test.
    Largest groups first (deterministic) so test holds rarer scaffolds -> a harder, leakage-free split.
    `seed` shuffles the order of equal-sized groups for repeated scaffold splits."""
    groups: dict[str, list[int]] = defaultdict(list)
    for i, s in enumerate(smiles):
        groups[murcko_scaffold(s)].append(i)
    rng = np.random.default_rng(seed)
    gl = list(groups.values())
    # deterministic: sort by (size desc, rng key) so ties break reproducibly per seed
    order = sorted(range(len(gl)), key=lambda j: (-len(gl[j]), rng.random()))
    n = len(smiles); n_train_target = int(round(frac_train * n))
    train: list[int] = []; test: list[int] = []
    for j in order:
        (train if len(train) + len(gl[j]) <= n_train_target else test).extend(gl[j])
    return np.sort(np.array(train, dtype=int)), np.sort(np.array(test, dtype=int))

def repeated_splits(smiles: list[str], kind: str = "random", n_repeats: int = 10,
                    frac_train: float = 0.8, base_seed: int = 0) -> Iterator[tuple[np.ndarray, np.ndarray]]:
    """Yield (train_idx, test_idx) for n_repeats deterministic resplits."""
    for r in range(n_repeats):
        seed = base_seed + r
        if kind == "random":
            yield random_split(len(smiles), frac_train, seed)
        elif kind == "scaffold":
            yield scaffold_split(smiles, frac_train, seed)
        else:
            raise ValueError(kind)

def three_way(smiles: list[str], kind: str, frac_train=0.7, frac_cal=0.15, seed: int = 0):
    """Train / calibration / test split (calibration carved from the train portion)."""
    if kind == "scaffold":
        tr, te = scaffold_split(smiles, frac_train + frac_cal, seed)
        sub = [smiles[i] for i in tr]
        a, b = scaffold_split(sub, frac_train / (frac_train + frac_cal), seed + 10_000)
    else:
        tr, te = random_split(len(smiles), frac_train + frac_cal, seed)
        a, b = random_split(len(tr), frac_train / (frac_train + frac_cal), seed + 10_000)
    return np.sort(tr[a]), np.sort(tr[b]), np.sort(te)
