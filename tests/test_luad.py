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


def test_clinical_baseline_uses_four_stage_levels_when_available():
    from luad_survival.evaluate import clinical_columns
    clin = _toy_clin(20)
    assert clinical_columns(clin) == ["age_years", "stage_advanced", "female"]          # fallback
    for k in ("II", "III", "IV"):
        clin[f"stage_{k}"] = 0.0
    assert clinical_columns(clin) == ["age_years", "stage_II", "stage_III", "stage_IV", "female"]


def test_proliferation_score_uses_training_statistics_and_needs_enough_genes():
    from luad_survival.proliferation import PROLIFERATION_GENES, proliferation_columns, proliferation_score
    ids = list(PROLIFERATION_GENES.values())[:4]
    rng = np.random.default_rng(0)
    X = pd.DataFrame(rng.normal(5, 1, (30, 5)), columns=[f"{i}.7" for i in ids] + ["ENSG_other"])
    cols = proliferation_columns(X)
    assert cols == [f"{i}.7" for i in ids]                                              # versions are ignored
    s = proliferation_score(X.iloc[:20], X.iloc[20:], cols)
    assert len(s) == 10 and abs(proliferation_score(X, X, cols).mean()) < 1e-9          # centred on its own training set
    assert proliferation_columns(X[["ENSG_other"]]) == []


def test_run_split_reports_proliferation_and_binary_baselines_when_possible():
    from luad_survival.proliferation import PROLIFERATION_GENES
    X, clin = synthetic(g=60)
    X = X.rename(columns=dict(zip(X.columns[:6], [f"{i}.1" for i in list(PROLIFERATION_GENES.values())[:6]])))
    r = run_split(X, clin, seed=0, n_genes=3, n_top_var=60)
    for k in ("c_proliferation", "c_clinical_plus_proliferation", "c_clinical_plus_proliferation_plus_signature",
              "c_clinical_binary_stage", "n_genes_fdr_lt_0.05_in_training"):
        assert k in r
    assert "c_proliferation" not in run_split(synthetic(g=60)[0], clin, seed=0, n_genes=3, n_top_var=60)


def test_permutation_null_is_centred_on_one_half_and_real_signal_beats_it():
    from luad_survival.evaluate import permutation_null
    X, clin = synthetic(seed=11)
    pn = permutation_null(X, clin, n_perm=20, n_splits=2, n_genes=3, n_top_var=60)
    assert abs(pn["null_mean_c"] - 0.5) < 0.05 and pn["observed_mean_c"] > pn["null_p95_c"] and pn["p_value"] < 0.1
    Xn = pd.DataFrame(np.random.default_rng(0).normal(5, 1, X.shape), index=X.index, columns=X.columns)
    assert permutation_null(Xn, clin, n_perm=20, n_splits=2, n_genes=3, n_top_var=60)["p_value"] > 0.05


def test_selection_frequency_finds_planted_genes_and_paired_gain():
    from luad_survival.evaluate import paired_gain, selection_frequency
    X, clin = synthetic(seed=2)
    f = selection_frequency(X, clin, n_splits=6, n_genes=3, n_top_var=60)
    assert set(f.index[:3]) == {"G0", "G1", "G2"} and f["fraction_of_splits"].iloc[0] == 1.0
    g = paired_gain(pd.DataFrame({"a": [0.7, 0.6, 0.65], "b": [0.6, 0.6, 0.7]}), "a", "b")
    assert g["mean"] == pytest.approx(0.05 / 3) and g["fraction_of_splits_better"] == pytest.approx(1 / 3)


def test_fdr_count_and_stage_level_cox():
    from luad_survival.evaluate import _fdr_count
    from luad_survival.eda import clinical_cox_stage_levels
    assert _fdr_count(pd.Series([1e-8, 1e-6, 0.5, 0.9, np.nan])) == 2
    clin = _toy_clin(300)
    level = np.random.default_rng(4).choice(["I", "II", "III", "IV"], 300)
    for k in ("II", "III", "IV"):
        clin[f"stage_{k}"] = (level == k).astype(float)
    tbl, info = clinical_cox_stage_levels(clin)
    assert {"Stage II vs I", "Stage III vs I", "Stage IV vs I"} <= set(tbl.index) and info["patients"] == 300


def test_loading_keeps_primary_tumours_one_row_per_patient_and_four_stage_levels(tmp_path):
    genes = ["ENSG1.1", "ENSG2.1"]
    expr = pd.DataFrame({"TCGA-AA-0001-01A": [1.0, 2.0], "TCGA-AA-0001-11A": [9.0, 9.0], "TCGA-AA-0002-01A": [3.0, 4.0],
                         "TCGA-AA-0002-01B": [7.0, 7.0], "TCGA-AA-0003-01A": [5.0, 6.0]}, index=genes)
    expr.index.name = "Ensembl_ID"
    expr.to_csv(tmp_path / "TCGA-LUAD.star_tpm.tsv.gz", sep="\t")
    pd.DataFrame({"sample": ["TCGA-AA-0001-01A", "TCGA-AA-0001-11A", "TCGA-AA-0002-01A", "TCGA-AA-0003-01A"],
                  "OS": [1, 1, 0, 1], "OS.time": [100, 100, 400, 0]}).to_csv(tmp_path / "TCGA-LUAD.survival.tsv.gz", sep="\t", index=False)
    pd.DataFrame({"sample": ["TCGA-AA-0001-01A", "TCGA-AA-0001-01A", "TCGA-AA-0001-11A", "TCGA-AA-0002-01A", "TCGA-AA-0003-01A"],
                  "age_at_diagnosis.diagnoses": [365.25 * 60] * 5, "gender.demographic": ["female", "female", "female", "male", "male"],
                  "ajcc_pathologic_stage.diagnoses": ["Stage IIIA", "Stage IIIA", "Stage IIIA", "Stage IB", "Stage IV"]}
                 ).to_csv(tmp_path / "TCGA-LUAD.clinical.tsv.gz", sep="\t", index=False)
    from luad_survival.data import build_dataset
    X, clin = build_dataset(tmp_path)
    assert list(clin.index) == ["TCGA-AA-0001", "TCGA-AA-0002"]               # normal tissue, repeated rows and OS.time = 0 are out
    assert X.loc["TCGA-AA-0001"].tolist() == [1.0, 2.0] and X.loc["TCGA-AA-0002"].tolist() == [3.0, 4.0]  # tumour values, first sample
    assert clin.loc["TCGA-AA-0001", ["stage_II", "stage_III", "stage_IV"]].tolist() == [0.0, 1.0, 0.0]
    assert clin.loc["TCGA-AA-0002", ["stage_II", "stage_III", "stage_IV"]].tolist() == [0.0, 0.0, 0.0]
