"""Cohort description and clinical analysis (patient level): the tables and figures of the exploratory notebook."""
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test, multivariate_logrank_test
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

STAGE_ORDER = ["Stage IA", "Stage IB", "Stage IIA", "Stage IIB", "Stage IIIA", "Stage IIIB", "Stage IV"]


def _km_median(time, event) -> float:
    m = KaplanMeierFitter().fit(time, event).median_survival_time_
    return float(m) if np.isfinite(m) else float("nan")


def cohort_overview(X: pd.DataFrame, clin: pd.DataFrame) -> pd.Series:
    """Headline numbers of the analysed cohort (one row per patient)."""
    return pd.Series({
        "patients": len(clin), "genes_measured": X.shape[1], "deaths": int(clin["OS"].sum()),
        "alive_or_censored": int((clin["OS"] == 0).sum()),
        "median_follow_up_days": float(clin["OS.time"].median()),
        "km_median_survival_days": _km_median(clin["OS.time"], clin["OS"]),
        "median_age_years": float(clin["age_years"].median()),
        "female_pct": float(100 * clin["female"].mean()),
    })


def stage_summary(clin: pd.DataFrame) -> pd.DataFrame:
    """Early (I-II) vs Advanced (III-IV): patients, median age, deaths, Kaplan-Meier median survival."""
    g = clin.groupby("stage_group")
    out = g.agg(patients=("OS", "size"), median_age_years=("age_years", "median"), deaths=("OS", "sum"),
                median_follow_up_days=("OS.time", "median"))
    out["km_median_survival_days"] = [_km_median(clin.loc[clin.stage_group == k, "OS.time"], clin.loc[clin.stage_group == k, "OS"]) for k in out.index]
    out["deaths_pct"] = (100 * out["deaths"] / out["patients"]).round(1)
    return out.round(1)


def stage_detail(clin: pd.DataFrame) -> pd.DataFrame:
    """Same, for each AJCC stage (IA ... IV)."""
    d = clin.dropna(subset=["stage_detail"]).copy()
    d["stage_detail"] = pd.Categorical(d["stage_detail"], STAGE_ORDER, ordered=True)
    g = d.groupby(["stage_group", "stage_detail"], observed=True)
    out = g.agg(patients=("OS", "size"), deaths=("OS", "sum"), median_follow_up_days=("OS.time", "median"))
    out["deaths_pct"] = (100 * out["deaths"] / out["patients"]).round(1)
    return out


def cox_table(df: pd.DataFrame, covariates: list[str], labels: dict | None = None, duration: str = "OS.time", event: str = "OS"):
    """Multivariable Cox model -> (table with hazard ratio, 95 % CI, p; dict with n, events, concordance, LR-test p)."""
    d = df[covariates + [duration, event]].dropna()
    cph = CoxPHFitter().fit(d, duration, event)
    s = cph.summary
    tbl = pd.DataFrame({"hazard_ratio": s["exp(coef)"], "ci95_low": s["exp(coef) lower 95%"], "ci95_high": s["exp(coef) upper 95%"], "p_value": s["p"]})
    if labels:
        tbl.index = [labels.get(i, i) for i in tbl.index]
    info = {"patients": int(len(d)), "events": int(d[event].sum()), "concordance": float(cph.concordance_index_),
            "likelihood_ratio_test_p": float(cph.log_likelihood_ratio_test().p_value)}
    return tbl, info


def clinical_cox(clin: pd.DataFrame):
    return cox_table(clin, ["age_years", "stage_advanced", "female"],
                     {"age_years": "Age (per year)", "stage_advanced": "Advanced vs early stage", "female": "Female vs male"})


def logrank_by_stage(clin: pd.DataFrame) -> float:
    a, b = clin[clin.stage_group == "Early"], clin[clin.stage_group == "Advanced"]
    return float(logrank_test(a["OS.time"], b["OS.time"], a["OS"], b["OS"]).p_value)


def cohens_d(a, b) -> float:
    a, b = np.asarray(a, float), np.asarray(b, float)
    sp = np.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) / (len(a) + len(b) - 2))
    return float((a.mean() - b.mean()) / sp)


def pca_by_stage(X: pd.DataFrame, clin: pd.DataFrame, n_genes: int = 2000):
    """Unsupervised PCA on the most variable genes (exploration only, nothing is fitted to survival).
    Returns (scores DataFrame with PC1, PC2, stage_group; explained variance; Cohen's d Early - Advanced for PC1 and PC2)."""
    Xf = X.loc[clin.index]
    top = Xf.var(axis=0).sort_values(ascending=False).head(n_genes).index
    Z = StandardScaler().fit_transform(Xf[top])
    pca = PCA(n_components=2, random_state=0)
    P = pca.fit_transform(Z)
    sc = pd.DataFrame(P, index=clin.index, columns=["PC1", "PC2"]); sc["stage_group"] = clin["stage_group"].values
    e, a = sc[sc.stage_group == "Early"], sc[sc.stage_group == "Advanced"]
    return sc, pca.explained_variance_ratio_, {"PC1": cohens_d(e.PC1, a.PC1), "PC2": cohens_d(e.PC2, a.PC2)}
