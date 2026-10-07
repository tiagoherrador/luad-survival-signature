import numpy as np
import pandas as pd
import pytest
from lifelines import CoxPHFitter

from luad_survival.cox import univariate_cox
from luad_survival.data import is_primary_tumour, patient_id, simplify_stage
from luad_survival.evaluate import c_index, run_split
from luad_survival.signature import fit_signature, select_features


def synthetic(n=300, g=60, n_signal=3, seed=0, effect=0.9):
    """Gene 0..n_signal-1 drive hazard (genes 0 and 1 risk, gene 2 protective); the rest is noise."""
    rng = np.random.default_rng(seed)
    X = rng.normal(5, 1, (n, g))
    lp = effect * (X[:, 0] - 5) + effect * (X[:, 1] - 5) - effect * (X[:, 2] - 5)
    t = rng.exponential(1000 * np.exp(-lp))
    c = rng.exponential(1500, n)
    time, event = np.minimum(t, c), (t <= c).astype(int)
    cols = [f"G{i}" for i in range(g)]
    X = pd.DataFrame(X, columns=cols, index=[f"P{i}" for i in range(n)])
    clin = pd.DataFrame({"OS.time": time, "OS": event, "age_years": rng.normal(65, 8, n),
                         "stage_advanced": rng.integers(0, 2, n).astype(float)}, index=X.index)
    return X, clin


def test_fast_cox_matches_lifelines():
    X, clin = synthetic(n=200, g=8)
    ours = univariate_cox(X.values, clin["OS.time"].values, clin["OS"].values)
    for j in (0, 2, 5):
        z = (X.iloc[:, j] - X.iloc[:, j].mean()) / X.iloc[:, j].std(ddof=0)
        df = pd.DataFrame({"x": z, "t": clin["OS.time"], "e": clin["OS"]})
        ref = CoxPHFitter().fit(df, "t", "e")
        assert ours.loc[j, "beta"] == pytest.approx(ref.params_["x"], abs=2e-3)
        assert ours.loc[j, "p"] == pytest.approx(ref.summary.loc["x", "p"], abs=1e-3, rel=1e-2)


def test_zero_variance_gene_is_nan():
    X = np.column_stack([np.ones(50), np.random.default_rng(1).normal(size=50)])
    out = univariate_cox(X, np.arange(1, 51, dtype=float), np.ones(50))
    assert np.isnan(out.loc[0, "beta"]) and not np.isnan(out.loc[1, "beta"])


def test_screening_recovers_planted_genes_with_right_signs():
    X, clin = synthetic()
    sig = fit_signature(X, clin["OS.time"], clin["OS"], n_genes=3, n_top_var=60)
    assert set(sig.genes) == {"G0", "G1", "G2"}
    signs = dict(zip(sig.genes, sig.signs))
    assert signs["G0"] == 1 and signs["G1"] == 1 and signs["G2"] == -1


def test_signature_generalises_on_signal_and_not_on_noise():
    X, clin = synthetic(seed=3)
    assert run_split(X, clin, seed=0, n_genes=3, n_top_var=60)["c_signature"] > 0.65
    # null calibration: pure-noise expression, many INDEPENDENT datasets -> mean held-out C-index ~ 0.5
    cs = []
    for k in range(25):
        Xk, ck = synthetic(seed=100 + k)
        Xn = pd.DataFrame(np.random.default_rng(k).normal(5, 1, Xk.shape), index=Xk.index, columns=Xk.columns)
        cs.append(run_split(Xn, ck, seed=0, n_genes=5, n_top_var=60)["c_signature"])
    assert abs(np.mean(cs) - 0.5) < 0.03


def test_test_set_never_influences_the_fitted_signature():
    X, clin = synthetic(seed=5)
    tr = X.index[:200]
    a = fit_signature(X.loc[tr], clin.loc[tr, "OS.time"], clin.loc[tr, "OS"], 5, 60)
    X2 = X.copy(); X2.loc[X.index[200:]] = 123.0  # corrupt the "test" rows
    b = fit_signature(X2.loc[tr], clin.loc[tr, "OS.time"], clin.loc[tr, "OS"], 5, 60)
    assert a.genes == b.genes and a.threshold == pytest.approx(b.threshold)


def test_c_index_direction():
    t = np.array([1, 2, 3, 4.]); e = np.ones(4)
    assert c_index(t, e, np.array([4, 3, 2, 1.])) == 1.0  # highest risk dies first
    assert c_index(t, e, np.array([1, 2, 3, 4.])) == 0.0


def test_select_features_drops_unexpressed_genes():
    X = pd.DataFrame({"low": np.zeros(10), "hi": np.arange(10, dtype=float) + 2, "flat": np.full(10, 3.0)})
    assert select_features(X, min_mean=1.0, n_top_var=5)[0] == "hi"
    assert "low" not in select_features(X, min_mean=1.0, n_top_var=5)


def test_stage_and_barcode_helpers():
    assert simplify_stage("Stage IIB") == "Early" and simplify_stage("Stage IIIA") == "Advanced"
    assert simplify_stage("Stage IV") == "Advanced" and np.isnan(simplify_stage(None))
    assert patient_id("TCGA-38-7271-01A") == "TCGA-38-7271"
    assert is_primary_tumour("TCGA-38-7271-01A") and not is_primary_tumour("TCGA-38-7271-11A")
