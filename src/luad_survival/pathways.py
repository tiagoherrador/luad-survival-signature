"""Pathway-level interpretation (needs internet and `pip install -e ".[gsea]"`):
preranked GSEA of the training-set survival ranking, and over-representation analysis (ORA) of the signature genes."""
import numpy as np
import pandas as pd

from .annotate import ensembl_to_symbol


def survival_ranking(screening: pd.DataFrame) -> pd.Series:
    """Score per gene = sign(log hazard ratio) x -log10(p), clipped to +-50. Positive = higher expression, worse survival."""
    s = screening.dropna(subset=["beta", "p"]).copy()
    return (np.sign(s["beta"]) * (-np.log10(s["p"].clip(lower=1e-300)))).clip(-50, 50)


def ranking_by_symbol(ranking: pd.Series, annotation: dict) -> pd.Series:
    """Map Ensembl IDs to symbols; if several IDs share a symbol keep the one with the largest |score|."""
    df = pd.DataFrame({"score": ranking})
    df["symbol"] = [annotation.get(i, {}).get("symbol") for i in df.index]
    df = df.dropna(subset=["symbol"])
    df["abs"] = df["score"].abs()
    df = df.sort_values("abs", ascending=False).drop_duplicates("symbol")
    return df.set_index("symbol")["score"].sort_values(ascending=False)


def run_gsea(rank_symbols: pd.Series, gene_sets: str = "MSigDB_Hallmark_2020", permutations: int = 1000, seed: int = 42) -> pd.DataFrame:
    import gseapy as gp
    res = gp.prerank(rnk=rank_symbols,
                     gene_sets=gene_sets, min_size=15, max_size=500, permutation_num=permutations, seed=seed, outdir=None, verbose=False)
    out = res.res2d[["Term", "NES", "NOM p-val", "FDR q-val", "Lead_genes"]].copy()
    for c in ("NES", "NOM p-val", "FDR q-val"):
        out[c] = pd.to_numeric(out[c], errors="coerce")
    return out.sort_values("FDR q-val").reset_index(drop=True)


def run_ora(symbols: list[str], background: list[str] | None = None,
            libraries=("MSigDB_Hallmark_2020", "GO_Biological_Process_2021", "KEGG_2021_Human")) -> pd.DataFrame:
    """Over-representation of the signature genes. `background` should be the genes that could have been chosen (here the
    screened pool); without it Enrichr compares against the whole genome, which overstates significance."""
    import gseapy as gp
    kw = {"background": list(background)} if background else {}
    res = gp.enrichr(gene_list=list(symbols), gene_sets=list(libraries), organism="human", outdir=None, **kw).results
    if "Overlap" not in res.columns:  # with a custom background Enrichr no longer reports it; count the shared genes instead
        res = res.assign(Overlap=res["Genes"].str.count(";") + 1)
    return res[["Gene_set", "Term", "Overlap", "P-value", "Adjusted P-value", "Genes"]].sort_values("Adjusted P-value").reset_index(drop=True)
