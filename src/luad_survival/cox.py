"""Vectorised univariate Cox screening (Breslow ties), one Newton-Raphson fit per gene, all genes at once.

Used to rank thousands of genes in milliseconds, which makes repeated train/test splits affordable.
Validated against `lifelines.CoxPHFitter` in the test-suite.
"""
import numpy as np
import pandas as pd
from scipy.stats import norm


def univariate_cox(X: np.ndarray, time: np.ndarray, event: np.ndarray, n_iter: int = 30) -> pd.DataFrame:
    """X [n, g] (any scale; standardised internally), time [n], event [n] in {0,1}.

    Returns a frame with `beta` (log hazard ratio **per SD** of expression), `hr`, `se`, `z`, `p`.
    Genes with zero variance get NaN.
    """
    X = np.asarray(X, dtype=float)
    time, event = np.asarray(time, float), np.asarray(event).astype(bool)
    sd = X.std(axis=0, ddof=0)
    ok = sd > 0
    Z = np.zeros_like(X)
    Z[:, ok] = (X[:, ok] - X[:, ok].mean(axis=0)) / sd[ok]

    order = np.argsort(-time, kind="stable")  # descending time -> cumulative sums are risk sets
    t, e, Zo = time[order], event[order], Z[order]
    last = np.searchsorted(-t, -t, side="right") - 1  # last row of each tie group (Breslow)
    beta = np.zeros(X.shape[1])
    info = np.ones_like(beta)
    for _ in range(n_iter):
        eta = Zo * beta
        w = np.exp(eta - eta.max(axis=0))  # per-gene shift cancels in the ratios below
        s0, s1, s2 = (np.cumsum(a, axis=0)[last] for a in (w, w * Zo, w * Zo**2))
        m = s1 / s0
        score = (Zo[e] - m[e]).sum(axis=0)
        info = (s2 / s0 - m**2)[e].sum(axis=0)
        step = np.clip(score / np.where(info > 0, info, np.inf), -1.0, 1.0)
        beta += step
        if np.abs(step).max() < 1e-9:
            break
    se = 1 / np.sqrt(np.where(info > 0, info, np.nan))
    z = beta / se
    out = pd.DataFrame({"beta": beta, "hr": np.exp(beta), "se": se, "z": z, "p": 2 * norm.sf(np.abs(z))})
    out.loc[~ok] = np.nan
    return out
