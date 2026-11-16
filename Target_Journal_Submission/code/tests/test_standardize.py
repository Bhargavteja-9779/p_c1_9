import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
from preprocessing.standardize import (NAMED_PIPELINES, standardization_orbit, canonical_identity,
                                       enumerate_tautomers, tautomer_score, OrbitConfig)
from rdkit import Chem


def test_pipelines_valid_or_none():
    for s in ["CC(=O)Oc1ccccc1C(=O)O", "[Na+].[Cl-].CCO", "Cc1cc(=O)[nH]c(=O)[nH]1"]:
        for p in NAMED_PIPELINES.values():
            out = p.apply(s)
            assert out is None or Chem.MolFromSmiles(out) is not None


def test_salt_stripping():
    assert NAMED_PIPELINES["P_rdkit"].apply("[Na+].[Cl-].CCO") == "CCO"


def test_orbit_nonempty_and_valid():
    for s in ["Cc1cc(=O)[nH]c(=O)[nH]1", "OC(=O)CCCCN"]:
        orb = standardization_orbit(s, OrbitConfig(max_tautomers=3, max_members=6))
        assert len(orb) >= 1
        assert all(Chem.MolFromSmiles(m) is not None for m in orb)


def test_orbit_handles_pathological_molecules():
    # the molecule that previously crashed the pool (carbanion salt) + junk
    for s in ["CC(=O)[c-]1c(=O)cc(C)oc1=O.[Na+]", "", "notasmiles", "[Na+].[Cl-]"]:
        orb = standardization_orbit(s, OrbitConfig())
        assert isinstance(orb, list)


def test_canonical_identity_collapses_orbit():
    s = "Cc1cc(=O)[nH]c(=O)[nH]1"
    ids = {canonical_identity(m) for m in standardization_orbit(s)}
    assert len([i for i in ids if i]) == 1  # all orbit members share one identity


def test_canonical_identity_distinguishes_molecules():
    assert canonical_identity("CCO") != canonical_identity("c1ccccc1")


def test_tautomer_score_ranks_canonical_first():
    m = Chem.MolFromSmiles("Oc1ccccn1")  # 2-hydroxypyridine; pyridone is canonical
    tauts = enumerate_tautomers(m, 4)
    # the highest-scored (first returned) should equal RDKit's canonical tautomer
    from rdkit.Chem.MolStandardize import rdMolStandardize
    canon = Chem.MolToSmiles(rdMolStandardize.TautomerEnumerator().Canonicalize(m))
    assert tauts[0] == canon
