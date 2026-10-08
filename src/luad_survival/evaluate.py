"""Held-out evaluation: C-index, log-rank, bootstrap CIs, clinical and proliferation baselines, repeated random splits,
an outcome-permutation null for the whole pipeline and gene-selection stability."""
from collections import Counter

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test
from lifelines.utils import concordance_index
from sklearn.model_selection import train_test_split

from .proliferation import MIN_GENES, proliferation_columns, proliferation_score
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


def clinical_columns(clin: pd.DataFrame) -> list[str]:
    """Clinical baseline: age, sex and AJCC stage in four levels (I, II, III, IV). A binary early/advanced stage throws
    information away and makes the baseline, and so the apparent value of any gene score, misleading."""
    stage = [c for c in ("stage_II", "stage_III", "stage_IV") if c in clin.columns]
    return ["age_years"] + (stage if stage else ["stage_advanced"]) + (["female"] if "female" in clin.columns else [])


def split_ids(clin: pd.DataFrame, seed: int, test_size: float = 0.3):
    return train_test_split(clin.index, test_size=test_size, random_state=seed, stratify=clin["OS"])


def _fdr_count(p: pd.Series, q: float = 0.05) -> int:
    """Benjamini-Hochberg: number of genes with FDR < q."""
    p = np.sort(p.dropna().values)
    m = len(p)
    bh = np.minimum.accumulate((p * m / np.arange(1, m + 1))[::-1])[::-1]
    return int((bh < q).sum())


def run_split(X: pd.DataFrame, clin: pd.DataFrame, seed: int = 42, test_size: float = 0.3,
              n_genes: int = 15, n_top_var: int = 2000, n_random: int = 0) -> dict:
    """One stratified split. Every data-dependent step (gene filter, screening, score normalisation,
    risk threshold, proliferation score, Cox baselines) is fitted on the training patients only."""
    tr, te = split_ids(clin, seed, test_size)
    sig = fit_signature(X.loc[tr], clin.loc[tr, "OS.time"], clin.loc[tr, "OS"], n_genes, n_top_var)
    ctr, cte = clin.loc[tr].copy(), clin.loc[te].copy()
    ctr["sig"], cte["sig"] = sig.score(X.loc[tr]), sig.score(X.loc[te])
    t, e = cte["OS.time"].values, cte["OS"].values
    grp = sig.group(X.loc[te])
    lr = logrank_test(t[grp == "High risk"], t[grp == "Low risk"], e[grp == "High risk"], e[grp == "Low risk"])
    clinical = clinical_columns(clin)
    rand = {}
    if n_random:
        rc = random_signature_c_indices(X.loc[tr], X.loc[te], sig.screening, t, e, n_random, n_genes, seed=seed)
        c_sig = c_index(t, e, cte["sig"])
        rand = {"random_c_mean": float(rc.mean()), "random_c_p95": float(np.percentile(rc, 95)),
                "signature_percentile_vs_random": float((rc < c_sig).mean() + 0.5 * (rc == c_sig).mean()), "_random_c": rc}
    pro = {}
    pcols = proliferation_columns(X)
    if len(pcols) >= MIN_GENES:
        ctr["pro"] = proliferation_score(X.loc[tr], X.loc[tr], pcols)
        cte["pro"] = proliferation_score(X.loc[tr], X.loc[te], pcols)
        pro = {"c_proliferation": c_index(t, e, cte["pro"]),
               "c_clinical_plus_proliferation": c_index(t, e, _cox_risk(ctr, cte, clinical + ["pro"])),
               "c_clinical_plus_proliferation_plus_signature": c_index(t, e, _cox_risk(ctr, cte, clinical + ["pro", "sig"])),
               "corr_signature_proliferation": float(np.corrcoef(cte["sig"], cte["pro"])[0, 1])}
    return {**rand, **pro,
        "seed": seed, "n_train": len(tr), "n_test": len(te), "events_test": int(e.sum()),
        "c_signature": c_index(t, e, cte["sig"]),
        "c_clinical": c_index(t, e, _cox_risk(ctr, cte, clinical)),
        "c_clinical_plus_signature": c_index(t, e, _cox_risk(ctr, cte, clinical + ["sig"])),
        "c_clinical_binary_stage": c_index(t, e, _cox_risk(ctr, cte, ["age_years", "stage_advanced"])),
        "logrank_p_high_vs_low": float(lr.p_value),
        "n_genes_fdr_lt_0.05_in_training": _fdr_count(sig.screening["p"]),
        "signature_genes": sig.genes, "_signature": sig, "_test": cte, "_group": grp,
    }


def repeated_splits(X, clin, n_splits: int = 30, **kw) -> pd.DataFrame:
    rows = []
    for sd in range(n_splits):
        r = run_split(X, clin, seed=sd, **kw)
        rows.append({k: v for k, v in r.items() if not k.startswith("_") and k != "signature_genes"})
    return pd.DataFrame(rows)


def paired_gain(rep: pd.DataFrame, better: str, baseline: str) -> dict:
    """Split-by-split difference in held-out C-index between two predictors."""
    d = rep[better] - rep[baseline]
    return {"mean": float(d.mean()), "sd": float(d.std()), "fraction_of_splits_better": float((d > 0).mean())}


def permutation_null(X, clin, n_perm: int = 100, n_splits: int = 5, seed: int = 0, **kw) -> dict:
    """Does the WHOLE pipeline (gene filter, screening, score, split) produce a good-looking C-index from noise?
    Survival time and event are shuffled together among patients (breaking any gene-outcome link, keeping the
    censoring structure) and the pipeline is re-run on the first `n_splits` splits. The observed mean C-index is compared
    with this null distribution; the p-value is (1 + #null >= observed) / (1 + n_perm)."""
    obs = float(np.mean([run_split(X, clin, seed=s, **kw)["c_signature"] for s in range(n_splits)]))
    rng = np.random.default_rng(seed)
    null = []
    for _ in range(n_perm):
        shuffled = clin.copy()
        shuffled[["OS.time", "OS"]] = clin[["OS.time", "OS"]].values[rng.permutation(len(clin))]
        null.append(np.mean([run_split(X, shuffled, seed=s, **kw)["c_signature"] for s in range(n_splits)]))
    null = np.array(null)
    return {"observed_mean_c": obs, "null_mean_c": float(null.mean()), "null_sd_c": float(null.std()),
            "null_p95_c": float(np.percentile(null, 95)), "n_permutations": n_perm, "n_splits_per_permutation": n_splits,
            "p_value": float((1 + (null >= obs).sum()) / (1 + n_perm)), "_null": null}


def selection_frequency(X, clin, n_splits: int = 30, n_genes: int = 15, n_top_var: int = 2000) -> pd.DataFrame:
    """How often each gene lands in the signature when only the training patients change (one row per gene ever chosen)."""
    count = Counter()
    for sd in range(n_splits):
        tr, _ = split_ids(clin, sd)
        count.update(fit_signature(X.loc[tr], clin.loc[tr, "OS.time"], clin.loc[tr, "OS"], n_genes, n_top_var).genes)
    out = pd.DataFrame({"times_selected": pd.Series(count)}).sort_values("times_selected", ascending=False)
    out["fraction_of_splits"] = out["times_selected"] / n_splits
    out.index.name = "ensembl_id"
    return out
