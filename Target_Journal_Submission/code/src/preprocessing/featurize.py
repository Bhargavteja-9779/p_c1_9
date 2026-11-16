"""Molecular featurization: ECFP (Morgan) fingerprints and RDKit physchem descriptors.

Deterministic and dependency-light (RDKit + NumPy). SMILES are standardized
(largest organic fragment, neutralized, canonicalized) before featurization so
that duplicate-detection and near-duplicate audits are well defined.
"""
from __future__ import annotations
from functools import lru_cache
from typing import Iterable, Sequence
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, DataStructs
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem import Descriptors

RDLogger.DisableLog("rdApp.*")
_LFC = rdMolStandardize.LargestFragmentChooser()
_UNCHARGER = rdMolStandardize.Uncharger()

def standardize_smiles(smiles: str) -> str | None:
    """Return a canonical SMILES for the largest, neutralized fragment, or None if unparseable."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        mol = _LFC.choose(mol)
        mol = _UNCHARGER.uncharge(mol)
        return Chem.MolToSmiles(mol, canonical=True)
    except Exception:
        return None

@lru_cache(maxsize=200_000)
def _mol(smiles: str):
    return Chem.MolFromSmiles(smiles)

def ecfp(smiles: str, radius: int = 2, n_bits: int = 2048) -> np.ndarray | None:
    """ECFP/Morgan bit vector as float32 array of length n_bits, or None if unparseable."""
    m = _mol(smiles)
    if m is None:
        return None
    gen = AllChem.GetMorganGenerator(radius=radius, fpSize=n_bits)
    fp = gen.GetFingerprint(m)
    arr = np.zeros((n_bits,), dtype=np.float32)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr

def ecfp_matrix(smiles_list: Sequence[str], radius: int = 2, n_bits: int = 2048) -> np.ndarray:
    """Stack ECFP vectors into an (n, n_bits) matrix. Unparseable SMILES -> zero row."""
    out = np.zeros((len(smiles_list), n_bits), dtype=np.float32)
    gen = AllChem.GetMorganGenerator(radius=radius, fpSize=n_bits)
    for i, s in enumerate(smiles_list):
        m = _mol(s)
        if m is None:
            continue
        fp = gen.GetFingerprint(m)
        DataStructs.ConvertToNumpyArray(fp, out[i])
    return out

# A compact, numerically stable set of 2D physicochemical descriptors.
_DESC_NAMES = [
    "MolWt", "MolLogP", "TPSA", "NumHAcceptors", "NumHDonors", "NumRotatableBonds",
    "NumAromaticRings", "RingCount", "FractionCSP3", "NumHeteroatoms", "HeavyAtomCount",
    "NHOHCount", "NOCount", "NumSaturatedRings", "NumAliphaticRings", "LabuteASA",
    "BalabanJ", "BertzCT", "Chi0", "Chi1", "HallKierAlpha", "Kappa1", "Kappa2", "Kappa3",
    "qed",
]
_DESC_FUNCS = {n: f for n, f in Descriptors.descList if n in set(_DESC_NAMES)}

def rdkit_descriptors(smiles_list: Sequence[str]) -> tuple[np.ndarray, list[str]]:
    """Return (n, d) descriptor matrix and the ordered descriptor names. NaN/inf -> 0."""
    names = [n for n in _DESC_NAMES if n in _DESC_FUNCS]
    out = np.zeros((len(smiles_list), len(names)), dtype=np.float32)
    for i, s in enumerate(smiles_list):
        m = _mol(s)
        if m is None:
            continue
        for j, n in enumerate(names):
            try:
                v = _DESC_FUNCS[n](m)
                out[i, j] = v if np.isfinite(v) else 0.0
            except Exception:
                out[i, j] = 0.0
    out = np.nan_to_num(out, nan=0.0, posinf=0.0, neginf=0.0)
    return out, names

def ecfp_bitset(smiles: str, radius: int = 2, n_bits: int = 2048):
    """RDKit ExplicitBitVect for fast Tanimoto similarity (used in leakage/cliff audits)."""
    m = _mol(smiles)
    if m is None:
        return None
    gen = AllChem.GetMorganGenerator(radius=radius, fpSize=n_bits)
    return gen.GetFingerprint(m)

def tanimoto(fp_a, fp_b) -> float:
    return float(DataStructs.TanimotoSimilarity(fp_a, fp_b))
