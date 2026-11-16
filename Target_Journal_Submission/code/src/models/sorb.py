"""SORB: Standardization-Orbit Robust learning.

Wraps any base estimator with (1) train-time augmentation over each molecule's standardization
orbit and (2) test-time prediction averaging over the orbit. The consistency-regularized MLP
variant lives in src/models/mlp_consistency.py.
"""
from __future__ import annotations
import numpy as np
from preprocessing.standardize import standardization_orbit, OrbitConfig, NAMED_PIPELINES
from preprocessing.featurize import ecfp_matrix, rdkit_descriptors
from models.zoo import build_model, predict_scores

ECFP_RADIUS, ECFP_BITS = 2, 2048


def _featurize(smiles_list, feat):
    if feat == "ecfp":
        return ecfp_matrix(smiles_list, ECFP_RADIUS, ECFP_BITS)
    X, _ = rdkit_descriptors(smiles_list)
    return X


class SORB:
    """Standardization-Orbit Robust estimator.

    Parameters
    ----------
    model_name, task, feat : base learner / task / featurizer
    train_aug : int  number of orbit members sampled per training molecule (1 = no augmentation)
    test_orbit : int number of orbit members averaged at test time (1 = canonicalize-once inference)
    seed : int
    orbit_cfg : OrbitConfig
    """

    def __init__(self, model_name="lgbm", task="reg", feat="ecfp", train_aug=4, test_orbit=8,
                 seed=0, orbit_cfg: OrbitConfig | None = None):
        self.model_name, self.task, self.feat = model_name, task, feat
        self.train_aug, self.test_orbit, self.seed = train_aug, test_orbit, seed
        self.orbit_cfg = orbit_cfg or OrbitConfig()
        self.model = None

    def _orbit_for(self, s, i, orbits):
        if orbits is not None:
            orb = list(orbits[i]) if orbits[i] else [NAMED_PIPELINES["P_rdkit"].apply(s) or s]
            return orb
        orb = standardization_orbit(s, self.orbit_cfg)
        return orb or [NAMED_PIPELINES["P_rdkit"].apply(s) or s]

    def fit(self, smiles, y, orbits=None):
        """orbits: optional precomputed list-of-lists (aligned to smiles) to avoid recomputing."""
        y = np.asarray(y, float)
        aug_smiles, aug_y = [], []
        for i, s in enumerate(smiles):
            orb = self._orbit_for(s, i, orbits)
            if self.train_aug <= 1:
                aug_smiles.append(orb[0]); aug_y.append(y[i])
            else:
                rng = np.random.default_rng(self.seed + i)
                k = self.train_aug
                idx = rng.choice(len(orb), size=k, replace=len(orb) < k)
                for j in idx:
                    aug_smiles.append(orb[int(j)]); aug_y.append(y[i])
        X = _featurize(aug_smiles, self.feat)
        self.model = build_model(self.model_name, self.task, self.seed)
        self.model.fit(X, np.asarray(aug_y))
        return self

    def predict(self, smiles, orbits=None):
        """Orbit-averaged prediction (values for reg, P(y=1) for clf)."""
        preds = np.zeros(len(smiles), float)
        counts = np.zeros(len(smiles), float)
        member_lists = []
        for i, s in enumerate(smiles):
            orb = self._orbit_for(s, i, orbits)
            member_lists.append(orb[0:1] if self.test_orbit <= 1 else orb[: self.test_orbit])
        # flatten, featurize once, scatter-average
        flat, owner = [], []
        for i, ml in enumerate(member_lists):
            for m in ml:
                flat.append(m); owner.append(i)
        X = _featurize(flat, self.feat)
        sc = predict_scores(self.model, X, self.task)
        for o, v in zip(owner, sc):
            preds[o] += v; counts[o] += 1
        counts[counts == 0] = 1
        return preds / counts
