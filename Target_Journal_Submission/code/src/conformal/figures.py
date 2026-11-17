"""Figures for the AD-conditional conformal study."""
from __future__ import annotations
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = os.path.join(os.path.dirname(__file__), "..", ".."); RES = os.path.join(ROOT, "results"); FIG = os.path.join(ROOT, "figures")
os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False, "savefig.bbox": "tight", "axes.grid": True, "grid.alpha": 0.25})
CB = {"marginal": "#D55E00", "normalized": "#999999", "knn_normalized": "#CC79A7", "mondrian_scaffold": "#E69F00", "smoothed_ad": "#56B4E9", "mondrian_ad": "#0072B2"}
LAB = {"marginal": "marginal (standard)", "normalized": "normalized (learned)", "knn_normalized": "distance-normalized", "mondrian_scaffold": "Mondrian-scaffold", "smoothed_ad": "smoothed-AD (ours)", "mondrian_ad": "Mondrian-AD (ours)"}

def _save(fig, n):
    for e in ("png", "pdf"): fig.savefig(os.path.join(FIG, f"{n}.{e}"), dpi=600)
    plt.close(fig)

def fig_nearfar(d):
    ds = d[(d.split == "scaffold") & (d.alpha == 0.1)]
    methods = [x for x in ["marginal", "normalized", "knn_normalized", "mondrian_scaffold", "smoothed_ad", "mondrian_ad"] if x in ds.method.unique()]
    near = ds.groupby("method").near_cov.mean(); far = ds.groupby("method").far_cov.mean()
    fig, ax = plt.subplots(figsize=(7, 4)); x = np.arange(len(methods)); w = 0.38
    ax.bar(x - w/2, [near[m] for m in methods], w, label="familiar (near AD)", color="#56B4E9")
    ax.bar(x + w/2, [far[m] for m in methods], w, label="novel (far AD)", color="#D55E00")
    ax.axhline(0.9, color="k", ls="--", lw=1, label="nominal 0.90")
    ax.set_xticks(x); ax.set_xticklabels([LAB[m] for m in methods], rotation=20, ha="right")
    ax.set_ylabel("coverage"); ax.set_ylim(0.5, 1.0)
    ax.set_title("Coverage for familiar vs novel chemotypes (scaffold split, mean over datasets)"); ax.legend(fontsize=8)
    _save(fig, "cfig_nearfar")

def fig_tradeoff(d):
    ds = d[(d.split == "scaffold") & (d.alpha == 0.1)]
    agg = ds.groupby("method")[["worst_bin_cov", "rel_width"]].mean()
    fig, ax = plt.subplots(figsize=(6, 4.5))
    for mth, r in agg.iterrows():
        ax.scatter(r.rel_width, r.worst_bin_cov, s=90, color=CB.get(mth, "#333"), zorder=3)
        ax.annotate(LAB.get(mth, mth), (r.rel_width, r.worst_bin_cov), fontsize=8, xytext=(5, 3), textcoords="offset points")
    ax.axhline(0.9, color="k", ls="--", lw=1); ax.set_xlabel("relative interval width (/ label SD)")
    ax.set_ylabel("worst-AD-bin coverage"); ax.set_title("Reliability–efficiency trade-off\n(up = more reliable for novel compounds)")
    _save(fig, "cfig_tradeoff")

def fig_perdataset(d):
    ds = d[(d.split == "scaffold") & (d.alpha == 0.1)]
    piv = ds.groupby(["dataset", "method"]).worst_bin_cov.mean().unstack()
    if {"marginal", "mondrian_ad"} <= set(piv.columns):
        p = piv[["marginal", "mondrian_ad"]].dropna()
        fig, ax = plt.subplots(figsize=(5, 5))
        ax.scatter(p["marginal"], p["mondrian_ad"], s=18, color="#0072B2", alpha=0.6)
        lim = [min(p.min()) - 0.02, 1.0]; ax.plot(lim, lim, "k--", lw=0.6); ax.axhline(0.9, color="#D55E00", ls=":", lw=0.8); ax.axvline(0.9, color="#D55E00", ls=":", lw=0.8)
        ax.set_xlim(lim); ax.set_ylim(lim); ax.set_xlabel("marginal conformal: worst-bin coverage"); ax.set_ylabel("Mondrian-AD: worst-bin coverage")
        ax.set_title("Per-dataset worst-AD-bin coverage\n(points above diagonal = improved)")
        _save(fig, "cfig_perdataset")

def fig_profile(path):
    """Coverage across all 5 AD strata (familiar->novel) for marginal vs AD-conditional methods."""
    if not os.path.exists(path):
        return
    p = pd.read_csv(path)
    methods = [m for m in ["marginal", "normalized", "smoothed_ad", "mondrian_ad"] if m in p.method.unique()]
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    for mth in methods:
        g = p[p.method == mth].groupby("stratum").coverage.mean()
        ax.plot(g.index + 1, g.values, marker="o", color=CB.get(mth, "#333"), label=LAB.get(mth, mth), lw=1.8)
    ax.axhline(0.9, color="k", ls="--", lw=1, label="nominal 0.90")
    ax.set_xticks(range(1, 6)); ax.set_xticklabels(["1\n(familiar)", "2", "3", "4", "5\n(novel)"])
    ax.set_xlabel("applicability-domain stratum"); ax.set_ylabel("coverage"); ax.set_ylim(0.6, 1.0)
    ax.set_title("Coverage across the applicability domain (scaffold split, mean over datasets)")
    ax.legend(fontsize=8, loc="lower left")
    _save(fig, "cfig_profile")


def fig_selective(path):
    """Does interval width track novelty? Mean relative interval width by AD stratum, per method.
    Marginal conformal is flat (width carries no triage signal); AD-conditional width rises with novelty."""
    if not os.path.exists(path):
        return
    d = pd.read_csv(path)
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    for mth in [m for m in ["marginal", "smoothed_ad", "mondrian_ad"] if m in d.method.unique()]:
        g = d[d.method == mth].groupby("stratum").rel_width.mean()
        ax.plot(g.index + 1, g.values, marker="o", color=CB.get(mth, "#333"), label=LAB.get(mth, mth), lw=1.8)
    ax.set_xticks(range(1, 6)); ax.set_xticklabels(["1\n(familiar)", "2", "3", "4", "5\n(novel)"])
    ax.set_xlabel("applicability-domain stratum"); ax.set_ylabel("relative interval width (/ label SD)")
    ax.set_title("Interval width as a novelty signal\n(marginal is flat; AD-conditional rises with novelty)")
    ax.legend(fontsize=8, loc="upper left")
    _save(fig, "cfig_selective")


def main():
    d = pd.read_csv(os.path.join(RES, "conformal_results.csv"))
    fig_nearfar(d); fig_tradeoff(d); fig_perdataset(d)
    fig_profile(os.path.join(RES, "conformal_profile_results.csv"))
    fig_selective(os.path.join(RES, "conformal_selective_results.csv"))
    print("conformal figures written")

if __name__ == "__main__":
    main()
