# LUAD – can the activity of a tumour's genes help predict how long lung-cancer patients live?

[![tests](https://github.com/tiagoherrador/luad-survival-signature/actions/workflows/ci.yml/badge.svg)](https://github.com/tiagoherrador/luad-survival-signature/actions)

> **Research and portfolio project — not a medical tool.** It was built on public research data and has never been tested in a clinic. It must not be used to estimate any real person's prognosis.

## The idea in plain words

Lung adenocarcinoma (LUAD) is the most common type of lung cancer. Doctors already estimate a patient's outlook with the **stage** of the disease: how far the tumour has spread when it is found.
But two patients with the same stage can do very differently, and the reason is partly written in the tumour itself: **which of its ~60,000 genes are switched on, and how strongly.**

This project uses the public **TCGA-LUAD** project, which measured gene activity in the tumours of hundreds of patients and then followed them over the years. We ask three questions:

1. **How much does stage tell us?** (the clinical baseline)
2. **Is there a combination of genes whose activity tracks survival** — and does it still work on patients that were *not* used to build it?
3. **What biology do those genes point to?**

Everything is done at the level of the **patient**, with the strict rule that a gene is only chosen using some patients (the *training* group) and judged on the others (the *held-out* group).

## What the analysis gives you

Running it produces **figures** and **tables** (all in [`results/`](results)). For example, the signature that came out of the training patients — 15 genes, with the direction in which each one moves survival
(*hazard ratio per SD*: above 1 means patients with more activity of that gene tend to live shorter; below 1, longer):

![The 15 genes of the signature](results/signature_genes.png)

| Gene | What it is | Effect | Hazard ratio per SD |
|---|---|---|---|
| DKK1 | dickkopf Wnt signaling pathway inhibitor 1 | worse survival | 1.82 |
| SACK1A | scaffolding CK1 anchoring protein A | worse survival | 1.67 |
| PKP2 | plakophilin 2 | worse survival | 1.61 |
| KRT6A | keratin 6A | worse survival | 1.57 |
| SERPINB5 | serpin family B member 5 | worse survival | 1.56 |
| GJB3 | gap junction protein beta 3 | worse survival | 1.56 |
| RHOV | ras homolog family member V | worse survival | 1.53 |
| FSCN1 | fascin actin-bundling protein 1 | worse survival | 1.51 |
| SDCBP2 | syndecan binding protein 2 | worse survival | 1.51 |
| LOXL2 | lysyl oxidase like 2 | worse survival | 1.51 |
| BMAL2 | basic helix-loop-helix ARNT like 2 | worse survival | 1.48 |
| LYPD3 | LY6/PLAUR domain containing 3 | worse survival | 1.46 |
| ABCC2 | ATP binding cassette subfamily C member 2 | worse survival | 1.41 |
| KRT14 | keratin 14 | worse survival | 1.40 |
| IRX5 | iroquois homeobox 5 | **better** survival | 0.65 |

The same table, with every p-value, is saved as a spreadsheet (`results/signature_genes.csv`). The other spreadsheet is **one line per held-out patient** (`results/test_risk_scores.csv`): the signature's risk score, whether the patient fell in the "high" or "low" risk group, and what actually happened to them.

| patient | follow-up (days) | died? | age | stage | risk score | risk group |
|---|---|---|---|---|---|---|
| TCGA-62-A471 | 1246 | no | 64.9 | Early | 2.03 | High risk |
| TCGA-55-7726 | 652 | no | 72.4 | Early | 1.79 | High risk |
| … | … | … | … | … | … | … |

## Following the analysis step by step

The steps follow the order of the original exploratory notebook (kept in [`notebooks/exploratory/`](notebooks/exploratory) with all its outputs).

### 1. Who are the patients?

After keeping **primary tumours with complete data, one row per patient**, the cohort has **474 patients**: 171 died during follow-up and 303 were alive or lost to follow-up ("censored").
Half of the patients were followed for less than ~1.8 years, and the typical (Kaplan–Meier) survival time is about 4 years.

| Cohort | |
|---|---|
| Patients | 474 |
| Genes measured | 60,660 |
| Deaths / still alive or censored | 171 / 303 |
| Median age at diagnosis | 66.6 years |
| Women | 54.9 % |
| Median follow-up | 651 days (1.8 years) |
| Median survival time (Kaplan–Meier) | 1,454 days (4.0 years) |

| Follow-up time | Age at diagnosis |
|---|---|
| ![follow-up](results/survival_time_hist.png) | ![age](results/age_hist.png) |

### 2. Stage matters — a lot

Patients were split into **early** (stages I–II) and **advanced** (stages III–IV) disease.

| Group | Patients | Deaths | Median survival (Kaplan–Meier) | Median age |
|---|---|---|---|---|
| Early (I–II) | 370 | 111 (30 %) | 4.4 years | 66.7 |
| Advanced (III–IV) | 104 | 60 (58 %) | 2.3 years | 66.3 |

![Survival by stage](results/km_by_stage.png)

*How to read a survival curve (Kaplan–Meier):* each step down is a death; the line shows the share of patients still alive over time, and the shaded band is the uncertainty. The further apart two curves are, the more the groups differ.
Here the curves separate early and stay apart (log-rank test p = 2×10⁻⁹): advanced stage means a clearly worse outlook, as clinical experience says.

| Stage | Patients | Deaths | % who died |
|---|---|---|---|
| IA | 125 | 27 | 21.6 % |
| IB | 128 | 36 | 28.1 % |
| IIA | 48 | 19 | 39.6 % |
| IIB | 63 | 27 | 42.9 % |
| IIIA | 69 | 38 | 55.1 % |
| IIIB | 10 | 6 | 60.0 % |
| IV | 25 | 16 | 64.0 % |

A statistical model that looks at age, stage and sex *together* (Cox regression) confirms that **stage is the main factor**; age and sex add little:

| Factor | Hazard ratio (95 % CI) | p-value |
|---|---|---|
| Advanced vs early stage | **2.55** (1.86 – 3.51) | 6×10⁻⁹ |
| Age (per extra year) | 1.01 (0.99 – 1.02) | 0.30 |
| Female vs male | 0.98 (0.73 – 1.33) | 0.92 |

*Hazard ratio:* how many times higher the risk of dying at any moment is, compared with the reference group. 2.55 means roughly two and a half times the risk. Stage alone ranks patients with a **C-index of about 0.62** (see the glossary below).

### 3. Does the whole gene-activity profile separate early from advanced tumours?

A *PCA* squeezes the 2,000 most variable genes into two summary axes, so that each patient becomes one dot. Nothing here is fitted to survival; it is just a look at the data.

![PCA](results/pca_by_stage.png)

Early and advanced tumours overlap a lot. The shift between the groups is small on the first axis (effect size 0.19) and small-to-moderate on the second (0.42): **stage leaves only a weak imprint on the overall gene profile**, so any survival signal in the genes will have to come from something more specific.

### 4. Building and testing a 15-gene signature

The patients were split at random, **70 % for training (331) and 30 % held out (143 patients, 52 deaths)**. Using *only the training patients*:

1. every gene is tested one by one for its link with survival;
2. the 15 strongest are kept, each with its direction (risk or protective);
3. a patient's **risk score** is the average of those genes' activity, each one counted in its direction and expressed relative to the training group.

Then the score is computed for the held-out patients, who split into a **high-risk** and a **low-risk** group at the training median:

![Held-out patients split by the signature](results/km_test_signature.png)

The high-risk group does worse, but the two curves overlap a lot and the difference is borderline (log-rank p = 0.058). Taking age and stage into account in the same model, the signature keeps a modest, statistically detectable effect:

| Factor (held-out patients only) | Hazard ratio (95 % CI) | p-value |
|---|---|---|
| Advanced vs early stage | 3.39 (1.94 – 5.91) | 1.7×10⁻⁵ |
| **Signature risk score (per SD)** | **1.36 (1.04 – 1.78)** | **0.025** |
| Age (per extra year) | 1.02 (0.99 – 1.05) | 0.15 |

### 5. How much better is the signature — and is it better than random genes?

A single random split can be lucky, so the whole procedure was **repeated on 30 different random splits**. For each one, three predictors were scored on the held-out patients with the **C-index**: stage and age alone, the signature alone, and both together.

| Held-out C-index | Average over 30 splits (± variation) | The seed-42 split above |
|---|---|---|
| Age + stage (clinical baseline) | 0.610 ± 0.044 | 0.668 |
| 15-gene signature alone | 0.617 ± 0.038 | 0.631 (95 % CI 0.538 – 0.716) |
| **Age + stage + signature** | **0.650 ± 0.038** | 0.661 |

![C-index over 30 splits](results/c_index_distribution.png)

Adding the signature to age and stage improves the average C-index by about **+0.04**, and does so in 27 of the 30 splits. That is a modest gain.

**The important check.** A well-known pitfall of gene signatures is that *almost any* set of genes shows some link to cancer survival, largely because the activity of many genes simply follows how fast the tumour is growing (shown for breast cancer by [Venet et al. 2011](docs/references.md)). So we asked: *does our carefully selected set beat sets of 15 random genes treated exactly the same way?*

![Our signature vs 200 random signatures](results/random_signature_comparison.png)

**Not clearly.** Random 15-gene sets reach an average C-index of about 0.59; ours is 0.63 on this split, better than ~80 % of the random sets, but on average over the 30 splits it beats 95 % of the random sets in only **3 splits out of 30**.
In other words, much of what the signature "predicts" is a general proliferation-related signal that many gene sets share; the *specific* choice of these 15 genes adds only a little. This is an honest, important limit of the result.

### 6. What biology do the genes point to?

Looking at the whole genome (not only the 15 genes), we ranked every gene by how strongly its activity relates to survival in the training patients and asked which **biological programmes** (the *Hallmark* gene collections) are concentrated at the "worse survival" end — a method called *gene-set enrichment analysis* (GSEA).

![Hallmark pathways](results/gsea_hallmark.png)

| Pathway | In plain words | Enrichment score | Adjusted p-value |
|---|---|---|---|
| Epithelial–mesenchymal transition | cells loosening up and becoming able to move and invade | 2.21 | < 0.001 |
| G2-M checkpoint | cells preparing to divide | 2.18 | < 0.001 |
| E2F targets | the cell-division programme | 2.07 | < 0.001 |
| Glycolysis | tumour feeding on sugar | 2.07 | < 0.001 |
| mTORC1 signalling | growth and metabolism switch | 1.97 | 0.0003 |
| Hypoxia | adapting to low oxygen | 1.96 | 0.0002 |
| TNF-α signalling via NF-κB | inflammation and stress response | 1.90 | 0.0004 |
| p53 pathway | stress and DNA-damage response | 1.90 | 0.0003 |

17 of the 28 pathways tested were significant, all in the same direction. The picture is coherent: **fast cell division, the ability to invade, and metabolic adaptation** are the traits that go with a worse outlook in LUAD — and several of the signature's own genes sit inside these programmes (for example DKK1 and LOXL2 in the invasion programme).

A second, narrower test asks whether the 15 signature genes *themselves* are enriched in some biological process. After correcting for the number of tests, the clearest hit is **intermediate-filament organisation** (2 of the 7 genes of that process, adjusted p = 0.002), which comes from the two keratin genes (KRT6A and KRT14). With only 15 genes this test has little power, so it should be read as an observation, not a finding.

> **Careful with the story.** The pathway picture says *what kind of biology* tracks survival. Because proliferation-related signals show up in almost any gene set, it does not by itself prove that *these 15 genes* are special (see step 5).

## How well does it work? (in short)

| Question | Answer |
|---|---|
| Does stage predict survival? | **Yes, strongly.** Advanced vs early: hazard ratio 2.55; the curves separate clearly. |
| Do genes add information beyond age and stage? | **A little.** About +0.04 in C-index on average (0.61 → 0.65); a modest, positive gain in 27 of 30 splits. |
| Is the 15-gene signature special? | **Not clearly.** It does not consistently beat random gene sets treated the same way. |
| What biology is involved? | **Cell division, invasion (EMT) and metabolic adaptation** — consistent with aggressive tumours. |

**What the metrics mean**
- **C-index (0.5 to 1)** — take two patients, one who died earlier than the other: how often does the predictor give the one who died earlier the higher risk? 0.5 is a coin flip, 1.0 is always right, and **0.6 to 0.7 is modest**.
- **Hazard ratio** — how many times higher the risk of dying is in one group than in another (1 = no difference).
- **Kaplan–Meier curve** — the share of patients still alive over time. **Log-rank p-value** — how likely a gap this large between two curves would be by chance alone.
- **Enrichment score and adjusted p-value (GSEA)** — how concentrated a biological programme is at the "worse survival" end of the gene ranking, and how likely that is by chance after correcting for testing many programmes.

## Things to keep in mind (honest limits)

- **It is a modest effect.** C-index values around 0.6–0.65 mean the genes help a bit; they are far from a reliable predictor for an individual.
- **No external check.** All patients come from the same cohort (TCGA). The signature has not been tried on an independent group of patients, which is the real test of a biomarker.
- **The 30 splits are not 30 independent studies.** They reshuffle the same 474 patients.
- **The signature is not clearly better than random gene sets** (step 5), and its biology is largely the general "fast-growing tumour" signal.
- **Only some patients have all the data.** Of the patients in the files, 474 had primary-tumour samples and complete age, stage and sex.
- **The size of the signature (15 genes) was not tuned**; it was inherited from the first version of the project. Other sizes might behave differently.
- **Survival here is "overall survival"**: death from any cause, not only from lung cancer.

## Differences from the original notebook

This repository follows the original exploratory notebook step by step, but re-does it at **patient level** and fixes problems found along the way. The notebook is kept, with all its outputs, in `notebooks/exploratory/`.

| In the original notebook | In this repository |
|---|---|
| Tables and Cox models counted **rows**, not patients: 878 observations (the survival file has 645 rows for 507 patients, and the merge repeated them) | One row per patient: **474 patients** |
| Normal-tissue samples mixed in with the tumours | **Primary tumours only** |
| Gene values log-transformed a second time (they were already on a log scale) | Values used as provided |
| Risk score = plain average of the 15 genes, so risk and protective genes **cancelled each other** | Each gene counted **in its own direction** |
| Gene filtering done with all patients | Done **on the training patients only** |
| One split, one number; signature effect reported per raw unit (hazard ratio 8.95, CI 1.02 – 78.8) | **30 splits**, comparison with random genes, and effect **per standard deviation** (1.36, CI 1.04 – 1.78) |
| "No significant enrichment" for the 15 genes | A small enrichment for keratin-related processes (2 genes) |
| Conclusion: the signature captures "a real biological axis, not a spurious association" | The **pathways** do point to a real biological axis; the **signature** shows only a modest, not clearly specific, added value |

## Try it yourself

Everything runs on **public data**; there is no trained model to request — the "model" is just the 15 genes above.

```bash
pip install -e ".[gsea,dev]"
# download the three TCGA-LUAD files from UCSC Xena into data/ (see data/README.md)

luad-run      --data-dir data --out results --gene-symbols   # tables and figures of steps 1-5
luad-pathways --data-dir data --out results                   # step 6 (needs internet)
pytest -q                                                      # automated tests
```
`--gene-symbols` and `luad-pathways` contact public web services (MyGene.info and Enrichr); only gene identifiers are sent, never patient data.

<details>
<summary><b>For developers</b></summary>

**Code layout** (`src/luad_survival/`): `data.py` (loading and harmonising the three Xena files: primary tumours, one row per patient) · `eda.py` (cohort, stage tables, Cox, PCA) · `cox.py` (vectorised univariate Cox, validated against `lifelines`) ·
`signature.py` (sign-aware score fitted on training patients only) · `evaluate.py` (held-out C-index, bootstrap CI, repeated splits, random-signature baseline) · `pathways.py` (GSEA and ORA) · `annotate.py` (gene symbols) · `plots.py` · `cli.py`, `pathways_cli.py`.

**Protocol.** 70/30 stratified split by outcome (seed 42 for the figures above; seeds 0–29 for the repeated analysis). Every data-dependent step — expression filter, top-2,000-variance gene selection, univariate screening, score normalisation, risk threshold (training median), Cox baselines — is fitted on training patients only.
The random-signature baseline draws 15 genes at random from the same screened pool, gives them the training-estimated directions and the same normalisation, and scores them on the held-out patients (200 draws per split).

**Tests** (`pytest -q`): the fast Cox fit matches `lifelines`; planted risk and protective genes are recovered with the right signs; on pure-noise expression the held-out C-index is 0.50 and the real signature is not an outlier among random ones; changing test-patient expression does not change the fitted signature (no leakage); stage tables count each patient once.

**Design decisions and the literature behind them:** [docs/design_decisions.md](docs/design_decisions.md) · [docs/references.md](docs/references.md).

</details>

## Author
Tiago Alejo Herrador — MIT licensed.
