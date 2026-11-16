"""Dataset registry and loaders for cached public benchmarks.

All loaders return a pandas DataFrame with canonical columns:
  smiles : str (raw SMILES as distributed)
  y      : float (regression target) or int (0/1 classification)
and, where available:
  relation : str in {'=','>','<'} for (future) censored-label experiments.

Data is expected under code/data_cache/ (see README for download instructions):
  - admet_group/<task>/{train_val.csv,test.csv}   (Therapeutics Data Commons ADMET benchmark group)
  - moleculeace_repo/MoleculeACE/Data/benchmark_data/<CHEMBLID>_<type>.csv
  - moleculenet/<name>.csv                         (optional; auto-downloaded on first use)
"""
from __future__ import annotations
import os, glob, io, urllib.request
from dataclasses import dataclass
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.abspath(os.path.join(HERE, "..", "..", "data_cache"))

# --- Therapeutics Data Commons ADMET benchmark group (22 datasets) ---------
TDC_REG = ["caco2_wang", "lipophilicity_astrazeneca", "solubility_aqsoldb", "ppbr_az",
           "vdss_lombardo", "half_life_obach", "clearance_hepatocyte_az",
           "clearance_microsome_az", "ld50_zhu"]
TDC_CLF = ["ames", "bbb_martins", "bioavailability_ma", "cyp2c9_veith", "cyp2d6_veith",
           "cyp3a4_veith", "cyp2c9_substrate_carbonmangels", "cyp2d6_substrate_carbonmangels",
           "cyp3a4_substrate_carbonmangels", "dili", "herg", "hia_hou", "pgp_broccatelli"]

def load_tdc(task: str) -> pd.DataFrame:
    base = os.path.join(CACHE, "admet_group", task)
    tv = pd.read_csv(os.path.join(base, "train_val.csv"))
    te = pd.read_csv(os.path.join(base, "test.csv"))
    tv["tdc_split"] = "train_val"; te["tdc_split"] = "test"
    df = pd.concat([tv, te], ignore_index=True)
    df = df.rename(columns={"Drug": "smiles", "Y": "y"})[["smiles", "y", "tdc_split"]]
    df["relation"] = "="
    return df

# --- MoleculeACE (30 ChEMBL regression datasets with activity cliffs) -------
def _mace_files():
    return sorted(glob.glob(os.path.join(CACHE, "moleculeace_repo", "MoleculeACE", "Data", "benchmark_data", "CHEMBL*_*.csv")))

def moleculeace_datasets() -> list[str]:
    return [os.path.basename(f)[:-4] for f in _mace_files()]

def load_moleculeace(name: str) -> pd.DataFrame:
    f = os.path.join(CACHE, "moleculeace_repo", "MoleculeACE", "Data", "benchmark_data", name + ".csv")
    df = pd.read_csv(f)
    # columns: smiles, exp_mean [nM], y (pEC50/pKi), cliff_mol (0/1), split (train/test)
    out = df.rename(columns={"y": "y_raw"})
    ycol = [c for c in df.columns if c.startswith("y [")]
    out["y"] = df[ycol[0]] if ycol else df["y"]
    out = out.rename(columns={"cliff_mol": "cliff", "split": "mace_split"})
    out["relation"] = "="
    keep = ["smiles", "y", "cliff", "mace_split"]
    return out[[c for c in keep if c in out.columns]]

# --- MoleculeNet physicochemical regression sets (auto-download) ------------
_MOLNET = {
    "esol": ("https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/delaney-processed.csv",
             "smiles", "measured log solubility in mols per litre"),
    "freesolv": ("https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/SAMPL.csv", "smiles", "expt"),
    "lipophilicity": ("https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/Lipophilicity.csv", "smiles", "exp"),
}

def load_moleculenet(name: str) -> pd.DataFrame:
    url, scol, ycol = _MOLNET[name]
    os.makedirs(os.path.join(CACHE, "moleculenet"), exist_ok=True)
    path = os.path.join(CACHE, "moleculenet", name + ".csv")
    if not os.path.exists(path):
        with urllib.request.urlopen(url, timeout=120) as r:
            open(path, "wb").write(r.read())
    df = pd.read_csv(path)
    out = df.rename(columns={scol: "smiles", ycol: "y"})[["smiles", "y"]].copy()
    out["relation"] = "="
    return out

@dataclass
class DatasetSpec:
    name: str
    loader: str      # 'tdc' | 'moleculeace' | 'moleculenet'
    key: str
    task: str        # 'reg' | 'clf'

def regression_registry() -> list[DatasetSpec]:
    specs = [DatasetSpec(f"tdc:{t}", "tdc", t, "reg") for t in TDC_REG]
    specs += [DatasetSpec(f"molnet:{m}", "moleculenet", m, "reg") for m in _MOLNET]
    specs += [DatasetSpec(f"mace:{d}", "moleculeace", d, "reg") for d in moleculeace_datasets()]
    return specs

def load(spec: DatasetSpec) -> pd.DataFrame:
    if spec.loader == "tdc":
        return load_tdc(spec.key)
    if spec.loader == "moleculeace":
        return load_moleculeace(spec.key)
    if spec.loader == "moleculenet":
        return load_moleculenet(spec.key)
    raise ValueError(spec.loader)
