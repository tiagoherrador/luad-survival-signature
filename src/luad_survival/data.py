"""Loading and harmonising the three UCSC Xena (GDC TCGA-LUAD) files.

* expression  : `TCGA-LUAD.star_tpm.tsv.gz`   genes x samples, already log2(TPM + 1)
* survival    : `TCGA-LUAD.survival.tsv.gz`   OS.time (days), OS (1 = death)
* clinical    : `TCGA-LUAD.clinical.tsv.gz`

Only **primary tumour** samples (TCGA sample type 01) are kept, one per patient.
Trimming barcodes to 12 characters and averaging duplicates would also mix in matched normal tissue, and applying
`log1p` on top of the already log-transformed values would distort the scale, so neither is done here.
"""
from pathlib import Path

import numpy as np
import pandas as pd

CLIN_COLS = ["age_at_diagnosis.diagnoses", "gender.demographic", "ajcc_pathologic_stage.diagnoses"]


def patient_id(barcode: str) -> str:
    return barcode[:12]


def is_primary_tumour(barcode: str) -> bool:
    return barcode[13:15] == "01"


def simplify_stage(stage) -> float | str:
    """AJCC stage -> 'Early' (I-II) / 'Advanced' (III-IV); NaN if unknown."""
    if stage is None or (isinstance(stage, float) and np.isnan(stage)):
        return np.nan
    stage = str(stage).strip()
    if stage.startswith(("Stage III", "Stage IV")):
        return "Advanced"
    if stage.startswith(("Stage I", "Stage II")):
        return "Early"
    return np.nan


def load_expression(path) -> pd.DataFrame:
    """Samples x genes matrix (log2(TPM+1)), primary tumours only, indexed by patient ID."""
    expr = pd.read_csv(path, sep="\t", index_col=0)
    keep = [c for c in expr.columns if is_primary_tumour(c)]
    expr = expr[keep].T
    expr.index = [patient_id(c) for c in expr.index]
    return expr[~expr.index.duplicated(keep="first")]


def load_clinical(survival_path, clinical_path) -> pd.DataFrame:
    surv = pd.read_csv(survival_path, sep="\t")
    clin = pd.read_csv(clinical_path, sep="\t")
    surv = surv[surv["sample"].map(is_primary_tumour)].copy()
    clin = clin[clin["sample"].map(is_primary_tumour)].copy()
    df = surv.merge(clin[["sample"] + CLIN_COLS], on="sample", how="inner")
    df["patient"] = df["sample"].map(patient_id)
    df["age_years"] = df["age_at_diagnosis.diagnoses"] / 365.25
    df["stage_group"] = df["ajcc_pathologic_stage.diagnoses"].map(simplify_stage)
    df["stage_advanced"] = (df["stage_group"] == "Advanced").astype(float).where(df["stage_group"].notna())
    df["female"] = (df["gender.demographic"] == "female").astype(float).where(df["gender.demographic"].notna())
    df = df.rename(columns={"ajcc_pathologic_stage.diagnoses": "stage_detail"})
    level = df["stage_detail"].str.extract(r"Stage (IV|III|II|I)")[0]  # I (reference), II, III, IV
    for k in ("II", "III", "IV"):
        df[f"stage_{k}"] = (level == k).astype(float).where(level.notna())
    df = df.drop_duplicates("patient").set_index("patient")
    return df[["OS.time", "OS", "age_years", "stage_group", "stage_detail", "stage_advanced", "stage_II", "stage_III", "stage_IV",
               "female", "gender.demographic"]]


def build_dataset(data_dir) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Returns (X patients x genes, clin) aligned on patient ID, complete cases for time/event/age/stage."""
    d = Path(data_dir)
    X = load_expression(d / "TCGA-LUAD.star_tpm.tsv.gz")
    clin = load_clinical(d / "TCGA-LUAD.survival.tsv.gz", d / "TCGA-LUAD.clinical.tsv.gz")
    clin = clin.dropna(subset=["OS.time", "OS", "age_years", "stage_advanced", "stage_II", "female"])
    clin = clin[clin["OS.time"] > 0]
    ids = X.index.intersection(clin.index)
    return X.loc[ids], clin.loc[ids]
