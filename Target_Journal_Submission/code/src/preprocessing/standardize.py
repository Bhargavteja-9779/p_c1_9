"""Molecular structure standardization: the orbit operations, named end-to-end pipelines,
and deterministic orbit sampling used throughout the study.

A *standardization operation* transforms a molecule along one axis (tautomer, protonation,
fragment/salt, stereochemistry, aromaticity). A *pipeline* is a fixed composition of operation
settings producing one canonical SMILES. The *standardization orbit* of a molecule is the set of
distinct SMILES reachable by varying the operation settings; SORB trains and predicts over it.

All functions are deterministic given their arguments. RDKit warnings are silenced.
Dimorphite-DL is optional; if unavailable, the pH-7.4 protonation axis degrades to RDKit neutralization
and this is recorded in `PROTONATION_BACKEND`.
"""
from __future__ import annotations
from dataclasses import dataclass, field
from functools import lru_cache
from typing import Sequence
from rdkit import Chem, RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize

RDLogger.DisableLog("rdApp.*")

# ---- optional pH-based protonation backend -------------------------------------------------
try:
    from dimorphite_dl import protonate_smiles as _dimorphite_protonate  # type: ignore
    PROTONATION_BACKEND = "dimorphite_dl"
except Exception:  # pragma: no cover
    _dimorphite_protonate = None
    PROTONATION_BACKEND = "rdkit_uncharger"

# ---- shared RDKit engines (stateless, reusable) --------------------------------------------
_LFC = rdMolStandardize.LargestFragmentChooser()
_UNCHARGER = rdMolStandardize.Uncharger()
_TE = rdMolStandardize.TautomerEnumerator()

_AROM_MODELS = {
    "rdkit": Chem.AromaticityModel.AROMATICITY_RDKIT,
    "mdl": Chem.AromaticityModel.AROMATICITY_MDL,
    "simple": Chem.AromaticityModel.AROMATICITY_SIMPLE,
}


def _parse(smiles: str):
    return Chem.MolFromSmiles(smiles)


def _to_smiles(mol, isomeric: bool = True) -> str | None:
    if mol is None:
        return None
    try:
        return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=isomeric)
    except Exception:
        return None


# ---- individual operations -----------------------------------------------------------------
def op_largest_fragment(mol):
    return _LFC.choose(mol)


def op_neutralize(mol):
    return _UNCHARGER.uncharge(mol)


def op_canonical_tautomer(mol):
    try:
        return _TE.Canonicalize(mol)
    except Exception:
        return mol


def op_strip_stereo(mol):
    m = Chem.Mol(mol)
    Chem.RemoveStereochemistry(m)
    return m


def op_set_aromaticity(mol, model: str = "rdkit"):
    m = Chem.Mol(mol)
    try:
        Chem.Kekulize(m, clearAromaticFlags=True)
        Chem.SetAromaticity(m, _AROM_MODELS[model])
    except Exception:
        return mol
    return m


def op_protonate_ph74(smiles: str) -> str:
    """Return the pH-7.4 dominant protonation SMILES (Dimorphite-DL) or RDKit-neutralized fallback."""
    if _dimorphite_protonate is not None:
        try:
            res = _dimorphite_protonate(smiles, ph_min=7.4, ph_max=7.4, pka_precision=1.0)
            if res and _parse(res[0]) is not None:
                return res[0]
        except Exception:
            pass
    m = _parse(smiles)
    if m is None:
        return smiles
    out = _to_smiles(op_neutralize(m))
    return out if out else smiles


def tautomer_score(mol) -> float:
    """RDKit's default tautomer plausibility score (ScoreRings + ScoreSubstructs + ScoreHeteroHs).
    This is exactly the scoring RDKit uses to pick the canonical tautomer; higher = more plausible.
    Used to keep only chemically reasonable tautomers and discard artifacts (e.g. carboxylic-acid enols)."""
    try:
        return float(rdMolStandardize.ScoreRings(mol) + rdMolStandardize.ScoreSubstructs(mol)
                     + rdMolStandardize.ScoreHeteroHs(mol))
    except Exception:
        return float("-inf")


def enumerate_tautomers(mol, max_forms: int = 4) -> list[str]:
    """The `max_forms` most plausible distinct tautomer SMILES, ranked by RDKit's tautomer score
    (descending). The canonical tautomer is always rank 1. Implausible high-energy tautomers
    (acid/amide enols etc.) are naturally downranked and excluded when max_forms is small."""
    try:
        scored = []
        seen = set()
        for t in _TE.Enumerate(mol):
            s = _to_smiles(t)
            if s in seen:
                continue
            seen.add(s)
            scored.append((tautomer_score(t), s))
        if not scored:
            return [_to_smiles(mol)]
        scored.sort(key=lambda x: -x[0])
        return [s for _, s in scored[:max_forms]]
    except Exception:
        return [_to_smiles(mol)]


# ---- named end-to-end pipelines ------------------------------------------------------------
@dataclass(frozen=True)
class Pipeline:
    name: str
    largest_fragment: bool = True
    neutralize: bool = False
    canonical_tautomer: bool = False
    strip_stereo: bool = False
    aromaticity: str = "rdkit"
    protonate_ph74: bool = False

    def apply(self, smiles: str) -> str | None:
        """Apply the pipeline to a raw SMILES; return canonical SMILES or None if unparseable."""
        if self.protonate_ph74:
            smiles = op_protonate_ph74(smiles)
        mol = _parse(smiles)
        if mol is None:
            return None
        try:
            if self.largest_fragment:
                mol = op_largest_fragment(mol)
            if self.neutralize:
                mol = op_neutralize(mol)
            if self.canonical_tautomer:
                mol = op_canonical_tautomer(mol)
            if self.strip_stereo:
                mol = op_strip_stereo(mol)
            if self.aromaticity != "rdkit":
                mol = op_set_aromaticity(mol, self.aromaticity)
            return _to_smiles(mol, isomeric=not self.strip_stereo)
        except Exception:
            return None


# The six realistic pipelines used for the sensitivity/ranking study (E1-E3).
NAMED_PIPELINES: dict[str, Pipeline] = {
    "P_raw": Pipeline("P_raw", largest_fragment=False),
    "P_rdkit": Pipeline("P_rdkit", largest_fragment=True, neutralize=True),
    "P_taut": Pipeline("P_taut", largest_fragment=True, neutralize=True, canonical_tautomer=True),
    "P_ph74": Pipeline("P_ph74", largest_fragment=True, protonate_ph74=True),
    "P_mdlarom": Pipeline("P_mdlarom", largest_fragment=True, neutralize=True, aromaticity="mdl"),
    "P_aggressive": Pipeline("P_aggressive", largest_fragment=True, neutralize=True,
                             canonical_tautomer=True, strip_stereo=True),
}

# For one-operation ablation (RQ3): each toggles a single axis off a common base (P_rdkit).
ABLATION_PIPELINES: dict[str, Pipeline] = {
    "base": NAMED_PIPELINES["P_rdkit"],
    "+tautomer": Pipeline("+tautomer", neutralize=True, canonical_tautomer=True),
    "+ph74": Pipeline("+ph74", protonate_ph74=True),
    "+stripstereo": Pipeline("+stripstereo", neutralize=True, strip_stereo=True),
    "+mdlarom": Pipeline("+mdlarom", neutralize=True, aromaticity="mdl"),
    "-neutralize": Pipeline("-neutralize", neutralize=False),
}


# ---- canonical identity (for leakage-safe splitting / dedup) --------------------------------
@lru_cache(maxsize=200_000)
def canonical_identity(smiles: str) -> str | None:
    """A single standardization-invariant key for a molecule: largest fragment, neutralized,
    canonical tautomer, stereo-stripped, InChIKey. Two orbit members of one molecule map here
    to the SAME key, so splitting on it prevents standardization-variant leakage across train/test."""
    mol = _parse(smiles)
    if mol is None:
        return None
    try:
        mol = op_largest_fragment(mol)
        mol = op_neutralize(mol)
        mol = op_canonical_tautomer(mol)
        mol = op_strip_stereo(mol)
        return Chem.MolToInchiKey(mol) or _to_smiles(mol, isomeric=False)
    except Exception:
        return None


# ---- standardization orbit (for SORB) ------------------------------------------------------
@dataclass
class OrbitConfig:
    max_tautomers: int = 4
    include_ph74: bool = True
    include_neutral: bool = True
    include_stripstereo: bool = True
    include_mdlarom: bool = True
    max_members: int = 8


def standardization_orbit(smiles: str, cfg: OrbitConfig | None = None) -> list[str]:
    """Deterministic set of distinct canonical SMILES reachable by varying standardization axes.
    The first element is always the canonicalize-once (P_rdkit) form, so orbit[0] == baseline input.
    Bounded to cfg.max_members."""
    cfg = cfg or OrbitConfig()
    base = NAMED_PIPELINES["P_rdkit"].apply(smiles)
    if base is None:
        base = smiles  # fall back to the raw SMILES as the single orbit member
    members: list[str] = [base]
    bmol = _parse(base)
    if bmol is None:
        return members

    def add(s):
        if s and s not in members:
            members.append(s)

    def safe(fn):
        """Run a producer returning a SMILES/None, swallowing any RDKit error."""
        try:
            return fn()
        except Exception:
            return None

    # tautomers
    try:
        for s in enumerate_tautomers(bmol, cfg.max_tautomers):
            add(s)
    except Exception:
        pass
    # protonation / charge
    if cfg.include_ph74:
        add(safe(lambda: op_protonate_ph74(smiles)))
    if cfg.include_neutral:
        add(safe(lambda: _to_smiles(op_neutralize(bmol))))
    # stereo stripped
    if cfg.include_stripstereo:
        add(safe(lambda: _to_smiles(op_strip_stereo(bmol), isomeric=False)))
    # aromaticity model
    if cfg.include_mdlarom:
        add(safe(lambda: _to_smiles(op_set_aromaticity(bmol, "mdl"))))
    return members[: cfg.max_members]


def sample_orbit(smiles: str, n: int, seed: int, cfg: OrbitConfig | None = None) -> list[str]:
    """Deterministically sample n members (with replacement if orbit smaller) for train augmentation."""
    import numpy as np
    orb = standardization_orbit(smiles, cfg)
    if not orb:
        return []
    rng = np.random.default_rng(seed)
    if len(orb) >= n:
        idx = rng.choice(len(orb), size=n, replace=False)
    else:
        idx = rng.choice(len(orb), size=n, replace=True)
    return [orb[i] for i in idx]
