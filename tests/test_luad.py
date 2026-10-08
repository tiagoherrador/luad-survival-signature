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


def test_random_signatures_do_not_beat_a_real_signal_but_look_like_noise_on_noise():
    from luad_survival.evaluate import run_split
    X, clin = synthetic(seed=7)
    r = run_split(X, clin, seed=0, n_genes=3, n_top_var=60, n_random=60)
    assert r["c_signature"] > r["random_c_mean"] and r["signature_percentile_vs_random"] > 0.9
    # on pure noise the real "signature" is just another random draw: its percentile is not systematically extreme
    pc = []
    for k in range(10):
        Xk, ck = synthetic(seed=200 + k)
        Xn = pd.DataFrame(np.random.default_rng(k).normal(5, 1, Xk.shape), index=Xk.index, columns=Xk.columns)
        pc.append(run_split(Xn, ck, seed=0, n_genes=5, n_top_var=60, n_random=40)["signature_percentile_vs_random"])
    assert 0.3 < np.mean(pc) < 0.7


def test_ensembl_to_symbol_strips_versions_and_survives_failures():
    from luad_survival.annotate import ensembl_to_symbol, strip_version
    assert strip_version("ENSG00000107984.10") == "ENSG00000107984"
    fake = lambda ids: [{"query": "ENSG00000107984", "symbol": "DKK1", "name": "dickkopf WNT signaling pathway inhibitor 1"}]
    out = ensembl_to_symbol(["ENSG00000107984.10", "ENSG00000000000.1"], fetch=fake)
    assert out["ENSG00000107984.10"]["symbol"] == "DKK1" and "ENSG00000000000.1" not in out
    def broken(ids):
        raise RuntimeError("offline")
    with pytest.warns(UserWarning):
        assert ensembl_to_symbol(["ENSG00000107984"], fetch=broken) == {}


def _toy_clin(n=200, seed=0):
    rng = np.random.default_rng(seed)
    adv = rng.integers(0, 2, n).astype(float)
    t = rng.exponential(1000 * np.exp(-0.9 * adv)); c = rng.exponential(1500, n)
    clin = pd.DataFrame({"OS.time": np.minimum(t, c), "OS": (t <= c).astype(int), "age_years": rng.normal(65, 8, n),
                         "stage_advanced": adv, "female": rng.integers(0, 2, n).astype(float)})
    clin["stage_group"] = np.where(adv == 1, "Advanced", "Early")
    clin["stage_detail"] = np.where(adv == 1, "Stage IIIA", "Stage IA")
    clin.index = [f"P{i}" for i in range(n)]
    return clin


def test_stage_tables_count_patients_once_and_cox_recovers_the_stage_effect():
    from luad_survival.eda import clinical_cox, stage_detail, stage_summary
    clin = _toy_clin()
    s = stage_summary(clin)
    assert s["patients"].sum() == len(clin) and set(s.index) == {"Early", "Advanced"}
    assert (s["deaths"] <= s["patients"]).all()
    assert stage_detail(clin)["patients"].sum() == len(clin)
    tbl, info = clinical_cox(clin)
    hr = tbl.loc["Advanced vs early stage"]
    assert 1.8 < hr["hazard_ratio"] < 3.5 and hr["p_value"] < 0.001 and info["patients"] == len(clin)


def test_pca_effect_size_has_the_right_sign_and_cohens_d():
    from luad_survival.eda import cohens_d, pca_by_stage
    assert cohens_d([1, 2, 3, 4], [3, 4, 5, 6]) < 0  # first group lower -> negative d
    clin = _toy_clin(120)
    rng = np.random.default_rng(1)
    X = pd.DataFrame(rng.normal(size=(120, 50)), index=clin.index, columns=[f"G{i}" for i in range(50)])
    X.iloc[:, :10] += (clin["stage_advanced"].values[:, None] * 3)  # strong stage signal
    sc, var, d = pca_by_stage(X, clin, n_genes=50)
    assert abs(d["PC1"]) > 1.0 and len(var) == 2


def test_survival_ranking_sign_clipping_and_symbol_mapping():
    from luad_survival.pathways import ranking_by_symbol, survival_ranking
    scr = pd.DataFrame({"beta": [0.5, -0.5, 0.1, np.nan], "p": [1e-4, 1e-100, 0.5, 0.01]}, index=["a", "b", "c", "d"])
    r = survival_ranking(scr)
    assert r["a"] > 0 and r["b"] == -50 and "d" not in r.index
    ann = {"a": {"symbol": "AAA"}, "b": {"symbol": "AAA"}, "c": {"symbol": "CCC"}}
    out = ranking_by_symbol(r, ann)
    assert out["AAA"] == -50 and out["CCC"] > 0  # duplicate symbol keeps the largest |score|
