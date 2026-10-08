"""    luad-pathways --data-dir data --out results      (needs internet and: pip install -e ".[gsea]")"""
import argparse
from pathlib import Path

from . import plots
from .annotate import ensembl_to_symbol
from .data import build_dataset
from .evaluate import run_split
from .pathways import ranking_by_symbol, run_gsea, run_ora, survival_ranking


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--out", default="results")
    ap.add_argument("--permutations", type=int, default=1000)
    a = ap.parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)

    X, clin = build_dataset(a.data_dir)
    sp = run_split(X, clin, seed=42)           # same training patients as luad-run (seed 42)
    sig = sp["_signature"]
    ids = list(sig.screening.dropna().index)
    ann = ensembl_to_symbol(ids)
    print(f"gene symbols found for {len(ann)} of {len(ids)} screened genes")
    rank = ranking_by_symbol(survival_ranking(sig.screening), ann)
    gsea = run_gsea(rank, permutations=a.permutations)
    gsea.to_csv(out / "gsea_hallmark.csv", index=False)
    plots.gsea_bar(gsea, out / "gsea_hallmark.png")
    print(gsea.head(10)[["Term", "NES", "FDR q-val"]].to_string(index=False))
    sym = [ann[g]["symbol"] for g in sig.genes if g in ann]
    ora = run_ora(sym)
    ora.to_csv(out / "ora_signature.csv", index=False)
    print("signature genes:", sym)
    print("ORA best adjusted p-value:", ora["Adjusted P-value"].min())


if __name__ == "__main__":
    main()
