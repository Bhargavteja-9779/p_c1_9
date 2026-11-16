import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
from data.splits import scaffold_split, random_split, murcko_scaffold
from preprocessing.featurize import ecfp_matrix, rdkit_descriptors, standardize_smiles
from preprocessing.standardize import canonical_identity
from models.zoo import build_model, predict_scores, MODEL_NAMES
from stat_tests.compare import paired_wilcoxon, bootstrap_ci, friedman, nemenyi_cd


def test_scaffold_split_disjoint():
    smiles = ["CCO", "CCN", "c1ccccc1", "c1ccccc1O", "c1ccncc1", "CC(=O)O", "CCCC", "c1ccc2ccccc2c1"] * 6
    tr, te = scaffold_split(smiles, 0.8, 0)
    assert len(set(tr) & set(te)) == 0
    strain = {murcko_scaffold(smiles[i]) for i in tr}
    stest = {murcko_scaffold(smiles[i]) for i in te}
    assert len(strain & stest) == 0


def test_featurization_shapes():
    s = ["CCO", "c1ccccc1O", "CC(=O)Oc1ccccc1C(=O)O"]
    X = ecfp_matrix(s); assert X.shape == (3, 2048)
    D, names = rdkit_descriptors(s); assert D.shape[0] == 3 and len(names) == D.shape[1]
    assert np.isfinite(X).all() and np.isfinite(D).all()


def test_models_fit_predict():
    rng = np.random.default_rng(0); X = rng.random((120, 32)); y = X @ rng.random(32)
    for n in MODEL_NAMES:
        m = build_model(n, "reg", 0); m.fit(X[:100], y[:100])
        p = predict_scores(m, X[100:], "reg"); assert p.shape == (20,) and np.isfinite(p).all()
        yc = (y > np.median(y)).astype(int)
        mc = build_model(n, "clf", 0); mc.fit(X[:100], yc[:100])
        pc = predict_scores(mc, X[100:], "clf"); assert ((pc >= 0) & (pc <= 1)).all()


def test_stats_sane():
    a = np.array([1.0, 2, 3, 4, 5]); b = a - 0.5
    w = paired_wilcoxon(a, b); assert 0 <= w["p"] <= 1
    lo = bootstrap_ci(a)[1]; hi = bootstrap_ci(a)[2]; assert lo <= np.mean(a) <= hi
    M = np.random.default_rng(0).random((10, 3)); f = friedman(M); assert 0 <= f["p"] <= 1
    assert nemenyi_cd(3, 10) > 0
