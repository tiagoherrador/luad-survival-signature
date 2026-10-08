"""    luad-run --data-dir data --out results --n-splits 30"""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from .annotate import ensembl_to_symbol
from . import eda, plots
from .data import build_dataset
from .evaluate import bootstrap_ci, paired_gain, permutation_null, repeated_splits, run_split, selection_frequency


def _round_keep_p(tbl, nd=4):
    """Round estimates but keep p-values at full precision (they can be far below 1e-4)."""
    t = tbl.copy()
    for c in t.columns:
        if not c.startswith("p_value"):
            t[c] = t[c].round(nd) if t[c].dtype.kind == "f" else t[c]
    return t


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--out", default="results")
    ap.add_argument("--n-splits", type=int, default=30)
    ap.add_argument("--n-genes", type=int, default=15)
    ap.add_argument("--n-random", type=int, default=200, help="random gene signatures per split for the comparison (0 = skip)")
    ap.add_argument("--n-perm", type=int, default=100, help="outcome permutations for the whole-pipeline null (0 = skip)")
    ap.add_argument("--gene-symbols", action="store_true", help="look up gene symbols on MyGene.info (needs internet; only gene IDs are sent)")
    a = ap.parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)

    X, clin = build_dataset(a.data_dir)
    print(f"patients={len(clin)} genes={X.shape[1]} deaths={int(clin['OS'].sum())}")

    # --- clinical / exploratory outputs, at patient level
    eda.cohort_overview(X, clin).to_frame("value").to_csv(out / "cohort_overview.csv")
    eda.stage_summary(clin).to_csv(out / "stage_summary.csv")
    eda.stage_detail(clin).to_csv(out / "stage_detail.csv")
    cox_tbl, cox_info = eda.clinical_cox(clin)
    _round_keep_p(cox_tbl).to_csv(out / "clinical_cox.csv")
    cox4, cox4_info = eda.clinical_cox_stage_levels(clin)
    _round_keep_p(cox4).to_csv(out / "clinical_cox_stage_levels.csv")
    lr_p = eda.logrank_by_stage(clin)
    plots.survival_time_hist(clin, out / "survival_time_hist.png")
    plots.age_hist(clin, out / "age_hist.png")
    plots.km_by_stage(clin, out / "km_by_stage.png", lr_p)
    pc_scores, pc_var, pc_d = eda.pca_by_stage(X, clin)
    plots.pca_stage(pc_scores, pc_var, pc_d, out / "pca_by_stage.png")
    clinical_summary = {"cox_clinical": cox_info, "cox_clinical_stage_levels": cox4_info, "logrank_early_vs_advanced_p": lr_p, "pca_cohens_d_early_minus_advanced": pc_d}

    main_split = run_split(X, clin, seed=42, n_genes=a.n_genes, n_random=a.n_random)
    te = main_split["_test"]
    lo, hi = bootstrap_ci(te["OS.time"], te["OS"], te["sig"])
    summary = {k: v for k, v in main_split.items() if not k.startswith("_")}
    summary["c_signature_ci95"] = [lo, hi]
    print(json.dumps(summary, indent=2))

    rep = repeated_splits(X, clin, n_splits=a.n_splits, n_genes=a.n_genes, n_random=a.n_random)
    rep.to_csv(out / "repeated_splits.csv", index=False)
    summary["repeated_splits"] = {c: {"mean": float(rep[c].mean()), "sd": float(rep[c].std()),
                                      "min": float(rep[c].min()), "max": float(rep[c].max())}
                                  for c in ("c_signature", "c_clinical", "c_clinical_binary_stage", "c_clinical_plus_signature",
                                            "c_proliferation", "c_clinical_plus_proliferation",
                                            "c_clinical_plus_proliferation_plus_signature") if c in rep.columns}
    summary["gains"] = {"signature_over_clinical": paired_gain(rep, "c_clinical_plus_signature", "c_clinical")}
    if "c_proliferation" in rep.columns:
        summary["gains"]["proliferation_over_clinical"] = paired_gain(rep, "c_clinical_plus_proliferation", "c_clinical")
        summary["gains"]["signature_over_clinical_plus_proliferation"] = paired_gain(
            rep, "c_clinical_plus_proliferation_plus_signature", "c_clinical_plus_proliferation")
        summary["repeated_splits"]["mean_correlation_signature_proliferation"] = float(rep["corr_signature_proliferation"].mean())
    summary["repeated_splits"]["median_genes_fdr_lt_0.05_in_training"] = float(rep["n_genes_fdr_lt_0.05_in_training"].median())
    summary["repeated_splits"]["frac_splits_logrank_p_lt_0.05"] = float((rep["logrank_p_high_vs_low"] < 0.05).mean())
    if a.n_random:
        summary["repeated_splits"]["random_signatures"] = {
            "mean_c_of_random_signatures": float(rep["random_c_mean"].mean()),
            "mean_percentile_of_our_signature": float(rep["signature_percentile_vs_random"].mean()),
            "fraction_of_splits_above_95th_percentile": float((rep["c_signature"] > rep["random_c_p95"]).mean())}
    plots.c_index_by_predictor(rep, out / "c_index_by_predictor.png")
    freq = selection_frequency(X, clin, n_splits=a.n_splits, n_genes=a.n_genes)
    summary["gene_selection_stability"] = {"distinct_genes_chosen": int(len(freq)), "chosen_in_every_split": int((freq["times_selected"] == a.n_splits).sum()),
                                           "chosen_in_at_least_half": int((freq["fraction_of_splits"] >= 0.5).sum())}
    if a.gene_symbols:
        ann_f = ensembl_to_symbol(list(freq.index))
        freq.insert(0, "gene_symbol", [ann_f.get(i, {}).get("symbol", "") for i in freq.index])
    freq.round(3).to_csv(out / "gene_selection_stability.csv")
    if a.n_perm:
        pn = permutation_null(X, clin, n_perm=a.n_perm, n_genes=a.n_genes)
        plots.permutation_null(pn.pop("_null"), pn["observed_mean_c"], out / "permutation_null.png")
        summary["permutation_null"] = pn
    summary["clinical"] = clinical_summary
    json.dump(summary, open(out / "metrics.json", "w"), indent=2)

    # --- tables: the signature genes and every held-out patient's risk score
    sig, te_tbl = main_split["_signature"], main_split["_test"]
    tbl = sig.screening.loc[sig.genes, ["hr", "p"]].rename(columns={"hr": "hazard_ratio_per_SD", "p": "p_value_train"}).copy()
    tbl.insert(0, "direction", ["risk (worse survival)" if sg > 0 else "protective" for sg in sig.signs])
    tbl.index.name = "ensembl_id"
    if a.gene_symbols:
        ann = ensembl_to_symbol(list(tbl.index))
        tbl.insert(0, "gene_symbol", [ann.get(i, {}).get("symbol", "") for i in tbl.index])
        tbl.insert(1, "gene_name", [ann.get(i, {}).get("name", "") for i in tbl.index])
    _round_keep_p(tbl, 3).to_csv(out / "signature_genes.csv")
    risk = te_tbl[["OS.time", "OS", "age_years", "stage_group", "sig"]].rename(columns={"sig": "risk_score"}).copy()
    risk["risk_group"] = main_split["_group"].values
    risk.index.name = "patient"
    risk.round(4).sort_values("risk_score", ascending=False).to_csv(out / "test_risk_scores.csv")
    tt = te_tbl.copy(); tt["risk_score_per_SD"] = (tt["sig"] - tt["sig"].mean()) / tt["sig"].std()
    ctab, cinfo = eda.cox_table(tt, ["age_years", "female", "stage_II", "stage_III", "stage_IV", "risk_score_per_SD"],
                                {"age_years": "Age (per year)", "female": "Female vs male", "stage_II": "Stage II vs I", "stage_III": "Stage III vs I",
                                 "stage_IV": "Stage IV vs I", "risk_score_per_SD": "Signature risk score (per SD)"})
    _round_keep_p(ctab).to_csv(out / "test_cox_with_signature.csv")
    summary["test_cox_with_signature"] = cinfo
    json.dump(summary, open(out / "metrics.json", "w"), indent=2)
    plots.km_test_groups(te_tbl, main_split["_group"], out / "km_test_signature.png", summary["logrank_p_high_vs_low"])
    plots.signature_genes(tbl, out / "signature_genes.png")

    if a.n_random:
        plt.figure(figsize=(6, 4))
        plt.hist(main_split["_random_c"], bins=20, color="grey", alpha=0.7, label=f"{a.n_random} random {a.n_genes}-gene signatures")
        plt.axvline(summary["c_signature"], color="crimson", lw=2, label=f"our signature ({summary['c_signature']:.3f})")
        plt.axvline(0.5, color="k", lw=1); plt.xlabel("Test C-index"); plt.legend(fontsize=8)
        plt.title("Is the signature better than random genes? (seed-42 split)"); plt.tight_layout()
        plt.savefig(out / "random_signature_comparison.png", dpi=150)


if __name__ == "__main__":
    main()
