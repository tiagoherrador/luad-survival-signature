"""Figures (English labels, same style for every plot)."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from lifelines import KaplanMeierFitter

RISK, PROTECT, NEUTRAL = "#d62728", "#2ca02c", "#4c72b0"


def _save(fig, path):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout(); fig.savefig(path, dpi=140); plt.close(fig)


def survival_time_hist(clin, path):
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(clin["OS.time"] / 365.25, bins=40, color=NEUTRAL)
    ax.set_xlabel("Follow-up time (years)"); ax.set_ylabel("Number of patients")
    ax.set_title("How long patients were followed (deaths and censored)")
    _save(fig, path)


def age_hist(clin, path):
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(clin["age_years"], bins=30, color=NEUTRAL)
    ax.set_xlabel("Age at diagnosis (years)"); ax.set_ylabel("Number of patients"); ax.set_title("Age at diagnosis")
    _save(fig, path)


def km_by_stage(clin, path, logrank_p):
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    for g, col in (("Early", PROTECT), ("Advanced", RISK)):
        ix = clin["stage_group"] == g
        KaplanMeierFitter().fit(clin.loc[ix, "OS.time"] / 365.25, clin.loc[ix, "OS"], label=f"{g} (n={int(ix.sum())})").plot_survival_function(ax=ax, color=col)
    ax.set_xlabel("Years since diagnosis"); ax.set_ylabel("Probability of being alive")
    ax.set_title(f"Survival by stage (log-rank p = {logrank_p:.1e})")
    _save(fig, path)


def pca_stage(scores, var, d, path):
    fig, ax = plt.subplots(figsize=(6, 4.8))
    for g, col in (("Early", PROTECT), ("Advanced", RISK)):
        s = scores[scores.stage_group == g]
        ax.scatter(s.PC1, s.PC2, s=14, alpha=0.6, color=col, label=g)
    ax.set_xlabel(f"PC1 ({100*var[0]:.0f} % of variance)"); ax.set_ylabel(f"PC2 ({100*var[1]:.0f} % of variance)")
    ax.set_title(f"Gene-expression PCA by stage (effect size d: PC1 {d['PC1']:.2f}, PC2 {d['PC2']:.2f})", fontsize=9); ax.legend()
    _save(fig, path)


def km_test_groups(test_tbl, groups, path, p):
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    for g, col in (("Low risk", PROTECT), ("High risk", RISK)):
        ix = (groups == g).values
        KaplanMeierFitter().fit(test_tbl.loc[ix, "OS.time"] / 365.25, test_tbl.loc[ix, "OS"], label=f"{g} (n={int(ix.sum())})").plot_survival_function(ax=ax, color=col)
    ax.set_xlabel("Years since diagnosis"); ax.set_ylabel("Probability of being alive")
    ax.set_title(f"Held-out patients split by the signature (log-rank p = {p:.3f})")
    _save(fig, path)


def signature_genes(tbl, path):
    t = tbl.copy().iloc[::-1]
    label = t["gene_symbol"].where(t.get("gene_symbol", "") != "", t.index) if "gene_symbol" in t else t.index
    fig, ax = plt.subplots(figsize=(6.5, 5))
    ax.barh(label, np.log(t["hazard_ratio_per_SD"]), color=[RISK if h > 1 else PROTECT for h in t["hazard_ratio_per_SD"]])
    ax.axvline(0, color="k", lw=1)
    ax.set_xlabel("log(hazard ratio) per 1 SD of expression (training set)\nred = worse survival, green = better survival")
    ax.set_title("The 15 genes of the signature")
    _save(fig, path)


def gsea_bar(gsea, path, top=12):
    t = gsea.head(top).sort_values("NES")
    fig, ax = plt.subplots(figsize=(7.5, 5))
    ax.barh(t["Term"], t["NES"], color=[RISK if v > 0 else PROTECT for v in t["NES"]])
    ax.axvline(0, color="k", lw=1)
    ax.set_xlabel("Enrichment score (higher = pathway more active in patients with worse survival)", fontsize=9)
    ax.set_title("Hallmark pathways linked to survival (GSEA, training patients)", fontsize=10)
    _save(fig, path)


PREDICTORS = [("c_clinical_binary_stage", "age + early/advanced stage"), ("c_clinical", "age + sex + stage I-IV"),
              ("c_signature", "15-gene signature"), ("c_proliferation", "proliferation score"),
              ("c_clinical_plus_signature", "clinical + signature"), ("c_clinical_plus_proliferation", "clinical + proliferation"),
              ("c_clinical_plus_proliferation_plus_signature", "clinical + proliferation + signature")]


def c_index_by_predictor(rep, path):
    cols = [(c, lab) for c, lab in PREDICTORS if c in rep.columns]
    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    ax.boxplot([rep[c] for c, _ in cols], vert=False, tick_labels=[lab for _, lab in cols], widths=0.6)
    ax.axvline(0.5, color="k", lw=1, ls="--")
    ax.invert_yaxis()
    ax.set_xlabel("C-index on held-out patients (30 random splits; 0.5 = coin flip)")
    ax.set_title("What predicts survival?")
    _save(fig, path)


def permutation_null(null, observed, path):
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(null, bins=20, color="grey", alpha=0.7, label="survival shuffled among patients")
    ax.axvline(observed, color=RISK, lw=2, label=f"real data ({observed:.3f})")
    ax.axvline(0.5, color="k", lw=1)
    ax.set_xlabel("Mean held-out C-index of the signature"); ax.set_ylabel("Number of shuffles"); ax.legend(fontsize=8)
    ax.set_title("Would the pipeline look good on pure noise?")
    _save(fig, path)
