import sys; sys.path.insert(0, 'src')
import numpy as np
from data import datasets as D
from preprocessing.standardize import NAMED_PIPELINES, standardization_orbit
from preprocessing.featurize import ecfp

PIPES = list(NAMED_PIPELINES.values())
reg = D.regression_registry()
pick = ([s for s in reg if s.loader == 'tdc'][:4] + [s for s in reg if s.loader == 'moleculenet'] +
        [s for s in reg if s.loader == 'moleculeace'][:4])
print(f"{'dataset':28} {'n':>5} {'%SMILESvary':>11} {'%ECFPvary':>10} {'meanOrbit':>9}")
tot_smi, tot_ecfp, tot_n = 0, 0, 0
for spec in pick:
    df = D.load(spec)
    sm = df.smiles.tolist()
    rng = np.random.default_rng(0)
    idx = rng.choice(len(sm), size=min(300, len(sm)), replace=False)
    smi_vary = ecfp_vary = orbit_sizes = 0
    for i in idx:
        s = sm[int(i)]
        outs = [p.apply(s) for p in PIPES]
        outs = [o for o in outs if o]
        if len(set(outs)) > 1:
            smi_vary += 1
        # ECFP variation across the distinct pipeline outputs
        fps = [ecfp(o) for o in set(outs)]
        fps = [f for f in fps if f is not None]
        if len(fps) > 1 and any(not np.array_equal(fps[0], f) for f in fps[1:]):
            ecfp_vary += 1
        orbit_sizes += len(standardization_orbit(s))
    n = len(idx)
    print(f"{spec.name:28} {n:5d} {100*smi_vary/n:10.1f}% {100*ecfp_vary/n:9.1f}% {orbit_sizes/n:9.2f}")
    tot_smi += smi_vary; tot_ecfp += ecfp_vary; tot_n += n
print(f"\nTOTAL across sampled molecules: n={tot_n}  SMILES-vary={100*tot_smi/tot_n:.1f}%  ECFP-vary={100*tot_ecfp/tot_n:.1f}%")
