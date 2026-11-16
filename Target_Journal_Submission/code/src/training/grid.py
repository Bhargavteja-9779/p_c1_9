"""Experiment grid definitions for quick and full modes."""
from __future__ import annotations
from data import datasets as D
from preprocessing.standardize import NAMED_PIPELINES, ABLATION_PIPELINES

MODELS = ["lgbm", "rf", "svm", "linear"]
SENSITIVITY_PIPELINES = list(NAMED_PIPELINES.keys())      # P_raw, P_rdkit, P_taut, P_ph74, P_mdlarom, P_aggressive
ABLATION_PIPE_NAMES = list(ABLATION_PIPELINES.keys())     # base, +tautomer, +ph74, +stripstereo, +mdlarom, -neutralize


def all_specs():
    reg = D.regression_registry()
    clf = [D.DatasetSpec(f"tdc:{t}", "tdc", t, "clf") for t in D.TDC_CLF]
    return reg + clf


def grid(mode: str):
    specs = all_specs()
    if mode == "quick":
        names = ["molnet:freesolv", "molnet:esol", "mace:CHEMBL1871_Ki", "tdc:bbb_martins"]
        specs = [s for s in specs if s.name in names]
        return dict(specs=specs, models=["lgbm", "rf", "linear"],
                    sens_pipes=["P_rdkit", "P_taut", "P_ph74"],
                    abl_pipes=["base", "+tautomer", "+ph74"],
                    seeds=[0, 1], method_models=["lgbm", "linear"], method_feats=["ecfp"],
                    method_seeds=[0, 1], split="scaffold")
    return dict(specs=specs, models=MODELS,
                sens_pipes=SENSITIVITY_PIPELINES, abl_pipes=ABLATION_PIPE_NAMES,
                seeds=[0, 1, 2, 3, 4], method_models=MODELS, method_feats=["ecfp", "desc"],
                method_seeds=[0, 1, 2, 3, 4], split="scaffold")
