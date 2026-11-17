"""Figure 1: conceptual schematic — a molecule's standardization orbit maps to one true property but
divergent model predictions (the invariance violation). Uses a real molecule + its real orbit, and a
real trained model's predictions, so the figure is data-grounded, not a cartoon."""
from __future__ import annotations
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec
from rdkit import Chem
from rdkit.Chem import Draw
from data import datasets as D
from data.splits import scaffold_split
from preprocessing.cache import features, feature_spec_to_pipe, dataset_orbits
from preprocessing.standardize import OrbitConfig
from preprocessing.featurize import ecfp_matrix
from models.zoo import build_model, predict_scores

FIG = os.path.join(os.path.dirname(__file__), "..", "..", "figures")


def main(dataset="mace:CHEMBL1871_Ki"):
    spec = [s for s in D.regression_registry() if s.name == dataset][0]
    df = D.load(spec); raw = df.smiles.tolist(); y = df.y.to_numpy().astype(float)
    tr, te = scaffold_split(raw, 0.8, 0)
    X = features(spec.name, raw, feature_spec_to_pipe("P_rdkit"), "ecfp")
    model = build_model("linear", "reg", 0); model.fit(X[tr], y[tr])
    orbits = dataset_orbits(spec.name, raw, OrbitConfig(max_tautomers=3, max_members=6))
    # pick a test molecule with a diverse orbit and non-trivial prediction spread
    cand = sorted(te, key=lambda i: -len(orbits[i]))
    pick = None
    for i in cand:
        orb = orbits[i]
        if len(orb) >= 3:
            Xo = ecfp_matrix(orb, 2, 2048); pr = predict_scores(model, Xo, "reg")
            if pr.std() > 0.1:
                pick = (i, orb, pr); break
    if pick is None:
        i = cand[0]; orb = orbits[i]; pr = predict_scores(model, ecfp_matrix(orb, 2, 2048), "reg")
        pick = (i, orb, pr)
    i, orb, pr = pick
    orb = orb[:4]; pr = pr[:4]
    mols = [Chem.MolFromSmiles(s) for s in orb]
    fig = plt.figure(figsize=(9, 4.2))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1.3, 1])
    # left: orbit molecule images
    img = Draw.MolsToGridImage(mols, molsPerRow=2, subImgSize=(200, 150),
                               legends=[f"form {k+1}" for k in range(len(mols))])
    axL = fig.add_subplot(gs[0]); axL.imshow(img); axL.axis("off")
    axL.set_title("One molecule, standardization orbit\n(chemically equivalent forms)", fontsize=10)
    # right: predictions diverge, true label fixed
    axR = fig.add_subplot(gs[1])
    axR.bar(range(len(pr)), pr, color="#56B4E9", edgecolor="#0072B2")
    axR.axhline(y[i], color="#D55E00", ls="--", lw=1.5, label=f"true value ({y[i]:.2f})")
    axR.set_xticks(range(len(pr))); axR.set_xticklabels([f"form {k+1}" for k in range(len(pr))])
    axR.set_ylabel("model prediction"); axR.set_title("Predictions differ\n(spread is spurious)", fontsize=10)
    axR.legend(fontsize=8)
    fig.suptitle("Standardization-invariance violation", y=1.0, fontsize=12)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(FIG, f"fig1_concept.{ext}"), dpi=600, bbox_inches="tight")
    plt.close(fig)
    print("fig1 written; example molecule orbit size", len(orb), "pred spread", round(float(np.std(pr)), 3))


if __name__ == "__main__":
    main()
