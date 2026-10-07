"""    luad-run --data-dir data --out results --n-splits 30"""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter

from .data import build_dataset
from .evaluate import bootstrap_ci, repeated_splits, run_split


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--out", default="results")
    ap.add_argument("--n-splits", type=int, default=30)
    ap.add_argument("--n-genes", type=int, default=15)
    a = ap.parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)

    X, clin = build_dataset(a.data_dir)
    print(f"patients={len(clin)} genes={X.shape[1]} deaths={int(clin['OS'].sum())}")

    main_split = run_split(X, clin, seed=42, n_genes=a.n_genes)
    te = main_split["_test"]
    lo, hi = bootstrap_ci(te["OS.time"], te["OS"], te["sig"])
    summary = {k: v for k, v in main_split.items() if not k.startswith("_")}
    summary["c_signature_ci95"] = [lo, hi]
    print(json.dumps(summary, indent=2))

    rep = repeated_splits(X, clin, n_splits=a.n_splits, n_genes=a.n_genes)
    rep.to_csv(out / "repeated_splits.csv", index=False)
    summary["repeated_splits"] = {c: {"mean": float(rep[c].mean()), "sd": float(rep[c].std()),
                                      "min": float(rep[c].min()), "max": float(rep[c].max())}
                                  for c in ("c_signature", "c_clinical", "c_clinical_plus_signature")}
    summary["repeated_splits"]["frac_splits_logrank_p_lt_0.05"] = float((rep["logrank_p_high_vs_low"] < 0.05).mean())
    json.dump(summary, open(out / "metrics.json", "w"), indent=2)

    kmf = KaplanMeierFitter()
    plt.figure(figsize=(6, 5))
    for g in ("Low risk", "High risk"):
        ix = main_split["_group"] == g
        kmf.fit(te.loc[ix.values, "OS.time"], te.loc[ix.values, "OS"], label=g).plot_survival_function()
    plt.title(f"Held-out test (seed 42): log-rank p={summary['logrank_p_high_vs_low']:.3f}")
    plt.xlabel("Days"); plt.ylabel("Survival probability"); plt.tight_layout()
    plt.savefig(out / "km_test.png", dpi=150)

    plt.figure(figsize=(6, 4))
    for c, lab in (("c_clinical", "age + stage"), ("c_signature", "signature"), ("c_clinical_plus_signature", "age + stage + signature")):
        plt.hist(rep[c], bins=12, alpha=0.5, label=lab)
    plt.axvline(0.5, color="k", lw=1); plt.xlabel("Test C-index"); plt.legend(); plt.title(f"{a.n_splits} random 70/30 splits")
    plt.tight_layout(); plt.savefig(out / "c_index_distribution.png", dpi=150)


if __name__ == "__main__":
    main()
