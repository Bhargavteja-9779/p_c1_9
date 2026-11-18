"""Post-process conformal_manuscript_filled.md for pandoc:
  - resolve @fig-*/@tbl-* crossrefs to 'Figure N'/'Table N'
  - embed the generated figures inline (after first reference) with numbered captions
  - emit conformal_manuscript_build.md ready for pandoc -> DOCX
Run after fill_conformal.py and src/conformal/figures.py.
"""
from __future__ import annotations
import os, re

HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.join(HERE, "..", "code", "figures")

# figure number, image file, caption
FIGS = [
    ("fig-nearfar", "cfig_nearfar.png",
     "Coverage for familiar (near applicability domain) versus novel (far) chemotypes under a 90% "
     "target (scaffold split, mean over datasets). Standard marginal conformal under-covers novel "
     "compounds; the normalized variant breaks even marginal coverage; Mondrian-AD and smoothed-AD "
     "restore novel-compound coverage to nominal."),
    ("fig-tradeoff", "cfig_tradeoff.png",
     "Reliability-efficiency trade-off. Worst applicability-domain-bin coverage (vertical, closer to "
     "0.90 is better) against relative interval width. Applicability-domain-conditional methods buy "
     "reliable novel-compound coverage with modestly wider intervals, placed where the model "
     "extrapolates."),
    ("fig-perdataset", "cfig_perdataset.png",
     "Per-dataset worst applicability-domain-bin coverage: Mondrian-AD (vertical) versus marginal "
     "conformal (horizontal). Points above the diagonal are datasets where conditioning on the "
     "applicability domain improves worst-bin coverage; dotted lines mark the 0.90 nominal target."),
    ("fig-mechanism", "cfig_mechanism.png",
     "Why learned difficulty fails to adapt to novelty. Growth of three quantities from the familiar to "
     "the novel applicability-domain stratum, normalized to the familiar stratum (median over datasets). "
     "The true prediction error grows substantially out of domain and the applicability-domain distance "
     "grows more, but the learned difficulty sigma(x) barely grows — so a normalized score built on it "
     "cannot widen novel-compound intervals enough."),
    ("fig-profile", "cfig_profile.png",
     "Coverage across all five applicability-domain strata, from familiar (1) to novel (5), scaffold "
     "split, mean over datasets. Marginal conformal slopes downward out of domain; the "
     "applicability-domain-conditional methods stay near the nominal 0.90 line."),
    ("fig-selective", "cfig_selective.png",
     "Interval width as a novelty signal: mean relative interval width by applicability-domain stratum "
     "(familiar to novel). Marginal conformal assigns a constant width (flat line), so width cannot flag "
     "novel compounds; Mondrian-AD and smoothed-AD widen with novelty, making width an honest triage "
     "signal."),
]
TBLS = [("tbl-summary", None)]


def main():
    src = open(os.path.join(HERE, "conformal_manuscript_filled.md")).read()
    # any unfilled tokens -> fail loudly
    leftover = re.findall(r"\{\{[A-Z0-9_]+\}\}", src)
    if leftover:
        raise SystemExit(f"UNFILLED TOKENS remain: {sorted(set(leftover))}; run fill_conformal.py first")
    fignum = {k: i + 1 for i, (k, _, _) in enumerate(FIGS)}
    tblnum = {k: i + 1 for i, (k, _) in enumerate(TBLS)}
    # resolve crossrefs ("Figure @fig-x" / "@fig-x" / "Table @tbl-x")
    for k, n in fignum.items():
        src = src.replace(f"Figure @{k}", f"Figure {n}").replace(f"@{k}", f"Figure {n}")
    for k, n in tblnum.items():
        src = src.replace(f"Table @{k}", f"Table {n}").replace(f"@{k}", f"Table {n}")
    # embed each figure inline after the first paragraph that references it
    for k, img, cap in FIGS:
        n = fignum[k]
        path = os.path.join(FIGDIR, img)
        block = f"\n\n![**Figure {n}.** {cap}]({path}){{width=85%}}\n"
        marker = f"Figure {n}"
        idx = src.find(marker)
        if idx != -1:
            para_end = src.find("\n\n", idx)
            if para_end == -1:
                para_end = len(src)
            src = src[:para_end] + block + src[para_end:]
    open(os.path.join(HERE, "conformal_manuscript_build.md"), "w").write(src)
    print(f"build md written; {len(FIGS)} figures embedded")


if __name__ == "__main__":
    main()
