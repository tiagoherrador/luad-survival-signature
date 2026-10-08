"""Held-out evaluation: C-index, log-rank, bootstrap CIs, clinical baseline and repeated random splits."""
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test
from lifelines.utils import concordance_index
from sklearn.model_selection import train_test_split

from .signature import fit_signature


def c_index(time, event, risk) -> float:
    """Harrell's C where a higher `risk` should mean shorter survival."""
    return float(concordance_index(time, -np.asarray(risk), event))


def bootstrap_ci(time, event, risk, n: int = 1000, seed: int = 0):
    rng = np.random.default_rng(seed)
    time, event, risk = map(np.asarray, (time, event, risk))
    vals = []
    for _ in range(n):
        i = rng.integers(0, len(time), len(time))
        if event[i].sum() > 0:
            vals.append(c_index(time[i], event[i], risk[i]))
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def _cox_risk(train: pd.DataFrame, test: pd.DataFrame, cols: list[str]) -> np.ndarray:
    cph = CoxPHFitter(penalizer=1e-3).fit(train[cols + ["OS.time", "OS"]], "OS.time", "OS")
    return cph.predict_partial_hazard(test[cols]).values


def random_signature_c_indices(X_train, X_test, screening, time_test, event_test, n_random, n_genes, seed=0):
    """Test C-index of `n_random` signatures made of RANDOM genes (drawn from the same screened pool).
    They get exactly the same treatment as the real signature (direction of each gene estimated on the training set,
    z-scoring with training mean/SD), so the comparison asks: is the *selection* of genes informative?"""
    rng = np.random.default_rng(seed)
    pool = screening.dropna().index.to_numpy()
    out = []
    for _ in range(n_random):
        g = list(rng.choice(pool, n_genes, replace=False))
        signs = np.sign(screening.loc[g, "beta"].values)
        z = (X_test[g] - X_train[g].mean()) / X_train[g].std(ddof=0)
        out.append(c_index(time_test, event_test, (z * signs).mean(axis=1)))
    return np.array(out)


def run_split(X: pd.DataFrame, clin: pd.DataFrame, seed: int = 42, test_size: float = 0.3,
              n_genes: int = 15, n_top_var: int = 2000, n_random: int = 0) -> dict:
    """One stratified split. Every data-dependent step (gene filter, screening, score normalisation,
    risk threshold, Cox baselines) is fitted on the training patients only."""
    tr, te = train_test_split(clin.index, test_size=test_size, random_state=seed, stratify=clin["OS"])
    sig = fit_signature(X.loc[tr], clin.loc[tr, "OS.time"], clin.loc[tr, "OS"], n_genes, n_top_var)
    ctr, cte = clin.loc[tr].copy(), clin.loc[te].copy()
    ctr["sig"], cte["sig"] = sig.score(X.loc[tr]), sig.score(X.loc[te])
    t, e = cte["OS.time"].values, cte["OS"].values
    grp = sig.group(X.loc[te])
    lr = logrank_test(t[grp == "High risk"], t[grp == "Low risk"], e[grp == "High risk"], e[grp == "Low risk"])
    clinical = ["age_years", "stage_advanced"]
    rand = {}
    if n_random:
        rc = random_signature_c_indices(X.loc[tr], X.loc[te], sig.screening, t, e, n_random, n_genes, seed=seed)
        c_sig = c_index(t, e, cte["sig"])
        rand = {"random_c_mean": float(rc.mean()), "random_c_p95": float(np.percentile(rc, 95)),
                "signature_percentile_vs_random": float((rc < c_sig).mean() + 0.5 * (rc == c_sig).mean()), "_random_c": rc}
    return {**rand, 
        "seed": seed, "n_train": len(tr), "n_test": len(te), "events_test": int(e.sum()),
        "c_signature": c_index(t, e, cte["sig"]),
        "c_clinical": c_index(t, e, _cox_risk(ctr, cte, clinical)),
        "c_clinical_plus_signature": c_index(t, e, _cox_risk(ctr, cte, clinical + ["sig"])),
        "logrank_p_high_vs_low": float(lr.p_value),
        "signature_genes": sig.genes, "_signature": sig, "_test": cte, "_group": grp,
    }


def repeated_splits(X, clin, n_splits: int = 30, **kw) -> pd.DataFrame:
    rows = []
    for sd in range(n_splits):
        r = run_split(X, clin, seed=sd, **kw)
        rows.append({k: v for k, v in r.items() if not k.startswith("_") and k != "signature_genes"})
    return pd.DataFrame(rows)
