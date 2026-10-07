"""Sign-aware gene-expression risk score fitted on a training set only."""
from dataclasses import dataclass

import numpy as np
import pandas as pd

from .cox import univariate_cox


def select_features(X_train: pd.DataFrame, min_mean: float = 1.0, n_top_var: int = 2000) -> list[str]:
    """Unsupervised filter computed on the TRAINING samples only: expressed genes, then most variable."""
    expressed = X_train.loc[:, X_train.mean(axis=0) > min_mean]
    return expressed.var(axis=0).sort_values(ascending=False).head(n_top_var).index.tolist()


@dataclass
class Signature:
    genes: list[str]
    signs: np.ndarray  # +1 risk gene, -1 protective gene
    mean: pd.Series
    sd: pd.Series
    threshold: float = 0.0
    screening: pd.DataFrame | None = None

    def score(self, X: pd.DataFrame) -> pd.Series:
        """Mean of z-scored expression (train mean/SD) multiplied by each gene's direction. Higher = worse."""
        z = (X[self.genes] - self.mean) / self.sd
        return (z * self.signs).mean(axis=1)

    def group(self, X: pd.DataFrame) -> pd.Series:
        return pd.Series(np.where(self.score(X) >= self.threshold, "High risk", "Low risk"), index=X.index)


def fit_signature(X_train: pd.DataFrame, time, event, n_genes: int = 15, n_top_var: int = 2000) -> Signature:
    """Univariate Cox screening on train; keep the `n_genes` smallest p-values; threshold = train median score."""
    feats = select_features(X_train, n_top_var=n_top_var)
    screen = univariate_cox(X_train[feats].values, np.asarray(time), np.asarray(event))
    screen.index = feats
    top = screen.dropna().sort_values("p").head(n_genes)
    genes = top.index.tolist()
    sig = Signature(genes=genes, signs=np.sign(top["beta"].values),
                    mean=X_train[genes].mean(), sd=X_train[genes].std(ddof=0), screening=screen)
    sig.threshold = float(sig.score(X_train).median())
    return sig
