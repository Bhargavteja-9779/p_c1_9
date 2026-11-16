# Standardization-Invariance in Molecular Property Prediction

Reproducible code for the study of **standardization-induced prediction instability** in molecular
property prediction and the **SORB** (Standardization-Orbit Robust) remedy.

> The same molecule, prepared with different but chemically-equivalent standardization choices
> (tautomer, protonation, salt/fragment, stereochemistry, aromaticity), can receive materially
> different model predictions — a per-prediction reproducibility failure that aggregate benchmark
> metrics hide. SORB enforces/learns invariance to this choice.

## 1. Environment
- Python 3.13, CPU-only (no GPU required). Tested on 4 cores / 15 GB RAM.
```bash
pip install -r requirements.txt
# PyTorch CPU wheel (only needed for the optional consistency-MLP variant):
pip install torch --index-url https://download.pytorch.org/whl/cpu
```
Or `conda env create -f environment.yml`.

## 2. Data (public, no credentials)
All datasets are public and are cached under `data_cache/` (git-ignored). To obtain them:
```bash
# TDC ADMET benchmark group (22 datasets) — Harvard Dataverse file 4426004
curl -L -o data_cache/admet_group.zip "https://dataverse.harvard.edu/api/access/datafile/4426004"
python -m zipfile -e data_cache/admet_group.zip data_cache/
# MoleculeACE (30 ChEMBL activity-cliff regression sets)
git clone --depth 1 https://github.com/molML/MoleculeACE.git data_cache/moleculeace_repo
# MoleculeNet ESOL/FreeSolv/Lipophilicity auto-download on first use (DeepChem S3).
```

## 3. Reproduce
```bash
python run_all.py --mode quick    # pipeline check (~3 min): small subset
python run_all.py --mode full     # full study (~2-3 h wall on 4 cores), checkpointed & resumable
python -m src.visualization.analyze    # -> results/metrics.json, statistics.json, tables/*.csv
python -m src.visualization.figures    # -> figures/*.{png,pdf}
python -m src.visualization.tables     # -> tables/*.{md,tex}
python scripts/leakage_audit.py        # -> results/leakage_audit.csv
pytest -q tests                        # unit tests
```
Re-running skips completed cells (`results/_cells/`). All randomness is seeded; the full study uses
seeds 0–4. Hardware auto-detected (CPU/CUDA/MPS).

## 4. Layout
```
src/
  data/          dataset loaders + leakage-safe splits (scaffold/random/three-way)
  preprocessing/ standardization orbit & pipelines (standardize.py), featurization, caching
  models/        model zoo (LightGBM/RF/Nystroem-SVM/linear) + SORB
  training/      experiment harness (sensitivity grid, ablation, invariance+method cell)
  evaluation/    metrics (regression/classification/calibration/interval)
  stat_tests/    Wilcoxon, bootstrap CIs, Cohen's dz, Friedman+Nemenyi, Holm-Bonferroni
  visualization/ analyze.py, figures.py, tables.py
experiments/  scripts/  tests/  results/  figures/  tables/
run_all.py    one-command reproduction
```

## 5. Key design points
- **Leakage-safe splits:** splitting is done on a standardization-invariant canonical identity
  (largest fragment + neutralized + canonical tautomer + stereo-stripped InChIKey), so no two
  standardization variants of one molecule cross the train/test boundary. Scaffold (Bemis–Murcko)
  splits are group-disjoint. See `scripts/leakage_audit.py`.
- **Standardization orbit:** RDKit `rdMolStandardize` (tautomers ranked by RDKit's tautomer score to
  keep only plausible forms) + Dimorphite-DL (pH 7.4). Defined in `src/preprocessing/standardize.py`.
- **Determinism:** fixed seeds; results cached to disk; figures/tables regenerated only from
  machine-readable result files (no hand-entered numbers).

## 6. Citation / license
Code released under the MIT License (see LICENSE). Datasets retain their original licenses
(TDC, MoleculeACE, MoleculeNet — all public, MIT/CC-BY).
