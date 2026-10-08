"""A plain cell-proliferation score, used as a competitor and as an adjustment for the signature.

Most gene-expression signatures partly measure how fast the tumour divides, so a new signature is only interesting
if it adds something beyond such a score. The genes are widely used cell-cycle / mitosis markers; Ensembl IDs are given
without version because the Xena files carry versioned IDs (``ENSG00000148773.15``).
"""
import pandas as pd

PROLIFERATION_GENES = {
    "MKI67": "ENSG00000148773", "TOP2A": "ENSG00000131747", "BIRC5": "ENSG00000089685", "CCNB1": "ENSG00000134057",
    "CDC20": "ENSG00000117399", "UBE2C": "ENSG00000175063", "AURKB": "ENSG00000178999", "TPX2": "ENSG00000088325",
    "CENPF": "ENSG00000117724", "MYBL2": "ENSG00000101057", "CCNA2": "ENSG00000145386", "PLK1": "ENSG00000166851",
}
MIN_GENES = 3


def proliferation_columns(X: pd.DataFrame) -> list[str]:
    """Columns of X (versioned or not) that correspond to the proliferation genes."""
    wanted = set(PROLIFERATION_GENES.values())
    return [c for c in X.columns if str(c).split(".")[0] in wanted]


def proliferation_score(X_train: pd.DataFrame, X: pd.DataFrame, cols: list[str]) -> pd.Series:
    """Mean z-score of the proliferation genes, standardised with the TRAINING mean / SD."""
    mu, sd = X_train[cols].mean(), X_train[cols].std(ddof=0).replace(0, 1.0)
    return ((X[cols] - mu) / sd).mean(axis=1)
