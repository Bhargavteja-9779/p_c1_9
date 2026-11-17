"""Publication figures generated from the result CSVs. Colorblind-safe palette; PNG (600 dpi) + PDF.
Run after run_all.py and analyze.py. Figures are regenerated solely from machine-readable results.
"""
from __future__ import annotations
import os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.join(os.path.dirname(__file__), "..", "..")
RES = os.path.join(ROOT, "results"); TAB = os.path.join(ROOT, "tables")
FIG = os.path.join(ROOT, "figures"); os.makedirs(FIG, exist_ok=True)

plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.dpi": 120, "savefig.bbox": "tight", "axes.grid": True,
                     "grid.alpha": 0.25, "grid.linewidth": 0.5})
CB = {"blue": "#0072B2", "orange": "#E69F00", "green": "#009E73", "red": "#D55E00",
      "purple": "#CC79A7", "yellow": "#F0E442", "sky": "#56B4E9", "grey": "#999999"}


def _save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(FIG, f"{name}.{ext}"), dpi=600)
    plt.close(fig)


DROP_DATASETS = {"molnet:lipophilicity"}
def _load(n):
    p = os.path.join(RES, n)
    if not os.path.exists(p): return None
    df = pd.read_csv(p)
    return df[~df.dataset.isin(DROP_DATASETS)].copy() if "dataset" in df.columns else df


def fig_invariance_by_representation(meth):
    """Headline: baseline standardization-invariance violation by model x featurization."""
    base = meth[(meth.method == "baseline") & (meth.task == "reg")].copy()
    base["mf"] = base.model + "\n" + base.feat
    order = base.groupby("mf").orbit_mean_rel_sd.median().sort_values().index
    fig, ax = plt.subplots(figsize=(7, 4))
    data = [base[base.mf == m].orbit_mean_rel_sd.dropna().values for m in order]
    bp = ax.boxplot(data, patch_artist=True, showfliers=False, widths=0.6)
    for patch in bp["boxes"]:
        patch.set_facecolor(CB["sky"]); patch.set_alpha(0.7)
    for i, m in enumerate(order, 1):
        v = base[base.mf == m].orbit_mean_rel_sd.dropna().values
        ax.scatter(np.random.default_rng(0).normal(i, 0.05, len(v)), v, s=8, color=CB["blue"], alpha=0.5, zorder=3)
    ax.set_xticks(range(1, len(order) + 1)); ax.set_xticklabels(order)
    ax.set_ylabel("Orbit prediction SD / label SD\n(spurious standardization sensitivity)")
    ax.set_xlabel("Model / featurization")
    ax.set_title("Standardization-invariance violation is representation-dependent")
    ax.axhline(0, color="k", lw=0.5)
    _save(fig, "fig_invariance_by_representation")


def fig_silent_contrast(sens, meth):
    """Aggregate metric is ~robust (small spread) but per-prediction invariance is violated."""
    rq1 = pd.read_csv(os.path.join(TAB, "rq1_sensitivity.csv")) if os.path.exists(os.path.join(TAB, "rq1_sensitivity.csv")) else None
    base = meth[(meth.method == "baseline") & (meth.task == "reg")]
    fig, axes = plt.subplots(1, 2, figsize=(9, 4))
    if rq1 is not None:
        axes[0].hist(rq1.rel_spread_pct.clip(upper=15), bins=30, color=CB["grey"], alpha=0.8)
        axes[0].axvline(rq1.rel_spread_pct.median(), color=CB["red"], ls="--",
                        label=f"median {rq1.rel_spread_pct.median():.2f}%")
        axes[0].set_xlabel("Aggregate metric spread across pipelines (%)")
        axes[0].set_ylabel("Count (dataset x model x feat)")
        axes[0].set_title("(a) Aggregate benchmark metric:\nnearly unchanged"); axes[0].legend()
    axes[1].hist(base.orbit_mean_rel_sd.dropna().clip(upper=0.6), bins=30, color=CB["orange"], alpha=0.85)
    axes[1].axvline(base.orbit_mean_rel_sd.median(), color=CB["red"], ls="--",
                    label=f"median {base.orbit_mean_rel_sd.median():.3f}")
    axes[1].set_xlabel("Per-prediction orbit SD / label SD")
    axes[1].set_title("(b) Individual predictions:\nmaterially irreproducible"); axes[1].legend()
    fig.suptitle("The silent confounder: aggregate metrics hide per-prediction instability", y=1.02)
    _save(fig, "fig_silent_contrast")


def fig_operation_ablation():
    p = os.path.join(TAB, "rq3_ablation.csv")
    if not os.path.exists(p):
        return
    df = pd.read_csv(p)
    med = df.groupby("operation").abs_delta.median().sort_values(ascending=False)
    fig, ax = plt.subplots(figsize=(6, 3.5))
    ax.barh(range(len(med)), med.values, color=CB["green"], alpha=0.8)
    ax.set_yticks(range(len(med))); ax.set_yticklabels(med.index); ax.invert_yaxis()
    ax.set_xlabel("Median |Δ primary metric| vs base")
    ax.set_title("Which standardization operation drives instability")
    _save(fig, "fig_operation_ablation")


def fig_method_effect(meth):
    """SORB: invariance reduction (train-aug) and accuracy effect, paired across datasets."""
    fig, axes = plt.subplots(1, 2, figsize=(9, 4))
    # invariance reduction: baseline vs trainaug orbit SD, per dataset, ecfp
    for ax, (model, color) in zip([axes[0]], [("lgbm", CB["blue"])]):
        pass
    b = meth[(meth.method == "baseline")].groupby(["dataset", "model", "feat"]).orbit_mean_rel_sd.mean()
    t = meth[(meth.method == "sorb_trainaug")].groupby(["dataset", "model", "feat"]).orbit_mean_rel_sd.mean()
    j = pd.concat([b.rename("baseline"), t.rename("trainaug")], axis=1).dropna().reset_index()
    axes[0].scatter(j.baseline, j.trainaug, s=12, c=[CB["blue"] if f == "ecfp" else CB["orange"] for f in j.feat], alpha=0.6)
    lim = [0, max(j.baseline.max(), j.trainaug.max()) * 1.05]
    axes[0].plot(lim, lim, "k--", lw=0.6); axes[0].set_xlim(lim); axes[0].set_ylim(lim)
    axes[0].set_xlabel("Baseline orbit SD / label SD"); axes[0].set_ylabel("SORB train-aug orbit SD / label SD")
    axes[0].set_title("(a) SORB reduces invariance violation\n(points below diagonal = improved)")
    # accuracy effect: distribution of improvement for each variant (reg, primary)
    reg = meth[meth.task == "reg"]
    improvements = {}
    base = reg[reg.method == "baseline"].groupby(["dataset", "model", "feat"]).rmse.mean()
    for v in ["sorb_trainaug", "sorb_testorbit", "sorb_full"]:
        vv = reg[reg.method == v].groupby(["dataset", "model", "feat"]).rmse.mean()
        d = (base - vv).dropna()  # positive = RMSE reduced = improvement
        improvements[v] = d.values
    axes[1].boxplot([improvements[v] for v in improvements], showfliers=False,
                    tick_labels=["train-aug", "test-orbit", "full"], patch_artist=True)
    axes[1].axhline(0, color=CB["red"], ls="--", lw=0.8)
    axes[1].set_ylabel("RMSE improvement vs baseline\n(>0 = better)")
    axes[1].set_title("(b) Accuracy effect of SORB variants")
    _save(fig, "fig_method_effect")


def main():
    meth = _load("method_results.csv"); sens = _load("raw_results.csv")
    if meth is not None:
        fig_invariance_by_representation(meth)
        fig_silent_contrast(sens, meth)
        fig_method_effect(meth)
    fig_operation_ablation()
    print("figures written to", FIG, ":", sorted(os.listdir(FIG)))


if __name__ == "__main__":
    main()
