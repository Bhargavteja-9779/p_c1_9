"""Generate manuscript tables (markdown + LaTeX) from the auto-computed table CSVs in tables/.
All values trace to results/*.csv via src/visualization/analyze.py."""
from __future__ import annotations
import os
import pandas as pd

ROOT = os.path.join(os.path.dirname(__file__), "..", "..")
TAB = os.path.join(ROOT, "tables")


def _emit(df, name, caption, float_fmt="%.3f"):
    df.to_csv(os.path.join(TAB, f"{name}.csv"), index=False)
    with open(os.path.join(TAB, f"{name}.md"), "w") as f:
        f.write(f"**{caption}**\n\n{df.to_markdown(index=False, floatfmt=float_fmt.replace('%',''))}\n")
    with open(os.path.join(TAB, f"{name}.tex"), "w") as f:
        f.write(df.to_latex(index=False, float_format=lambda x: float_fmt % x, escape=True,
                            caption=caption, label=f"tab:{name}"))


def table_dataset_summary():
    """Table 1: datasets used."""
    import sys
    sys.path.insert(0, os.path.join(ROOT, "code", "src"))
    sys.path.insert(0, os.path.join(ROOT, "src"))
    from data import datasets as D
    rows = []
    for spec in D.regression_registry() + [D.DatasetSpec(f"tdc:{t}", "tdc", t, "clf") for t in D.TDC_CLF]:
        df = D.load(spec)
        rows.append({"Dataset": spec.name, "Task": spec.task, "N": len(df),
                     "Source": spec.loader})
    t = pd.DataFrame(rows)
    _emit(t, "table1_datasets", "Datasets used in this study.", float_fmt="%d")


def table_invariance_by_representation():
    p = os.path.join(TAB, "rq5_baseline_invariance_by_representation.csv")
    if os.path.exists(p):
        t = pd.read_csv(p).rename(columns={"orbit_mean_rel_sd": "Orbit SD / label SD",
                                           "pipe_mean_rel_sd": "Pipeline SD / label SD",
                                           "model": "Model", "feat": "Featurization"})
        _emit(t, "table2_invariance", "Baseline standardization-invariance violation by representation "
                                      "(median over datasets).")


def table_method():
    p = os.path.join(TAB, "rq4_method_accuracy.csv")
    if os.path.exists(p):
        t = pd.read_csv(p)
        _emit(t, "table3_method", "SORB variants vs baseline: paired accuracy improvement across datasets.")


def main():
    table_dataset_summary()
    table_invariance_by_representation()
    table_method()
    print("tables written to", TAB)


if __name__ == "__main__":
    main()
