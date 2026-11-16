"""Disk cache of standardized SMILES and feature matrices per (dataset, pipeline).

Standardizing with Dimorphite-DL (pH 7.4) is the compute bottleneck, so each dataset is
standardized once per pipeline and the resulting SMILES + ECFP + descriptor matrices are
cached to .npz under data_cache/std_cache/. All experiments then reuse the cache.
"""
from __future__ import annotations
import os, hashlib, json
import numpy as np
from preprocessing.standardize import NAMED_PIPELINES, ABLATION_PIPELINES, Pipeline
from preprocessing.featurize import ecfp_matrix, rdkit_descriptors

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE_DIR = os.path.abspath(os.path.join(HERE, "..", "..", "data_cache", "std_cache"))
os.makedirs(CACHE_DIR, exist_ok=True)

ECFP_RADIUS, ECFP_BITS = 2, 2048


def _key(dataset: str, pipe: str, feat: str) -> str:
    h = hashlib.md5(f"{dataset}|{pipe}|{feat}|r{ECFP_RADIUS}b{ECFP_BITS}".encode()).hexdigest()[:16]
    return os.path.join(CACHE_DIR, f"{h}.npz")


def standardized_smiles(dataset: str, raw_smiles: list[str], pipe: Pipeline) -> list[str]:
    """Return pipeline-standardized SMILES (cached). Unparseable -> '' placeholder."""
    path = _key(dataset, pipe.name, "smiles")
    if os.path.exists(path):
        return list(np.load(path, allow_pickle=True)["smiles"])
    out = [pipe.apply(s) or "" for s in raw_smiles]
    np.savez_compressed(path, smiles=np.array(out, dtype=object))
    return out


def features(dataset: str, raw_smiles: list[str], pipe: Pipeline, feat: str = "ecfp") -> np.ndarray:
    """Return the (n, d) feature matrix for a dataset under a pipeline (cached). feat in {ecfp, desc}."""
    path = _key(dataset, pipe.name, feat)
    if os.path.exists(path):
        return np.load(path)["X"]
    smi = standardized_smiles(dataset, raw_smiles, pipe)
    # molecules that failed to standardize: fall back to raw so rows stay aligned
    smi = [s if s else raw for s, raw in zip(smi, raw_smiles)]
    if feat == "ecfp":
        X = ecfp_matrix(smi, ECFP_RADIUS, ECFP_BITS)
    elif feat == "desc":
        X, _ = rdkit_descriptors(smi)
    else:
        raise ValueError(feat)
    np.savez_compressed(path, X=X.astype(np.float32))
    return X


def feature_spec_to_pipe(pipe_name: str) -> Pipeline:
    if pipe_name in NAMED_PIPELINES:
        return NAMED_PIPELINES[pipe_name]
    if pipe_name in ABLATION_PIPELINES:
        return ABLATION_PIPELINES[pipe_name]
    raise KeyError(pipe_name)


# ---- orbit cache (Dimorphite is slow; compute each molecule's orbit once per dataset) --------
from preprocessing.standardize import standardization_orbit, OrbitConfig


def _orbit_key(dataset: str, cfg: OrbitConfig) -> str:
    tag = f"{dataset}|orbit|t{cfg.max_tautomers}m{cfg.max_members}ph{int(cfg.include_ph74)}" \
          f"n{int(cfg.include_neutral)}s{int(cfg.include_stripstereo)}a{int(cfg.include_mdlarom)}"
    h = hashlib.md5(tag.encode()).hexdigest()[:16]
    return os.path.join(CACHE_DIR, f"orbit_{h}.npz")


def dataset_orbits(dataset: str, raw_smiles: list[str], cfg: OrbitConfig | None = None) -> list[list[str]]:
    """List (aligned to raw_smiles) of each molecule's standardization-orbit SMILES. Cached to disk;
    building it is the only place Dimorphite-DL is invoked, and only once per dataset."""
    cfg = cfg or OrbitConfig()
    path = _orbit_key(dataset, cfg)
    if os.path.exists(path):
        arr = np.load(path, allow_pickle=True)["orbits"]
        return [list(o) for o in arr]
    orbits = [standardization_orbit(s, cfg) for s in raw_smiles]
    np.savez_compressed(path, orbits=np.array([np.array(o, dtype=object) for o in orbits], dtype=object))
    return orbits
