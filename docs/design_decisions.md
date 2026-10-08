# Design decisions

Why I built the analysis the way I did. Sources are in [references.md](references.md); where there is **no** literature behind a choice, I say so.

## 1. Everything data-dependent is fitted on training patients only
Gene filtering, screening, score normalisation, the risk threshold, the proliferation score and the clinical baselines all use only the training patients; the held-out patients are touched once, at the end.
*Why:* if genes are chosen using the same patients that are later tested, the error looks smaller than it really is (Ambroise & McLachlan 2002). A test (`test_test_set_never_influences_the_fitted_signature`) corrupts the held-out patients' expression and checks that the signature does not change.

## 2. Many random splits, not one
A single 70/30 split is one draw; a lucky or unlucky draw can mislead. I repeat the whole procedure on 30 random splits and report the average and spread.
*Why:* in a re-analysis of seven published studies, gene lists were highly unstable across training sets and five of seven did not classify better than chance; the authors recommend validation by repeated random sampling (Michiels et al. 2005).
*Limit:* the 30 splits reuse the same 474 patients, so they are not independent evidence. An external cohort is still missing.

## 3. C-index as the headline metric, with a bootstrap interval
The C-index ranks pairs of patients by who died first, handles censored follow-up, and is easy to read (0.5 = coin flip, 1 = perfect). Intervals come from resampling the held-out patients.
*Why:* Harrell et al. (1996) describe it as an interpretable index of predictive discrimination for survival data and recommend validating it by bootstrap or cross-validation.

## 4. Compare against a proper clinical baseline, a proliferation score *and* random genes
A signature is only useful if it adds something to what is already known, and a baseline that is too weak makes any gene score look good. With age and early/advanced stage only, the baseline is 0.61 and the signature seems to add about +0.04; with the four stage levels (I–IV) plus age and sex the baseline rises to 0.66 and the gain shrinks to about +0.01.
So I report: (a) the clinical baseline, (b) the signature alone, (c) both together, (d) a simple 12-gene proliferation score alone and added to the baseline, and (e) the C-index of 200 random 15-gene signatures per split, treated exactly like the real one.
*Why:* most random gene sets show some association with cancer outcome, largely through proliferation (Venet et al. 2011, in breast cancer); and molecular scores should be judged on what they add beyond stage and the other standard clinical factors (Subramanian & Simon 2010; Shedden et al. 2008 found that combining clinical and molecular information worked best).
*Result:* the signature adds ≈ +0.01 C-index to the proper baseline, the proliferation score adds ≈ +0.02, and the signature adds nothing (−0.007) on top of proliferation. It does not consistently beat random sets either (above the 95th percentile in 3 of 30 splits). I report this as it came out.
*Limit:* the proliferation genes are a standard cell-cycle marker list I chose, not an optimised or published signature.

## 5. Shuffle the survival data and rerun everything
To know what "good" looks like when there is nothing to find, I permute survival time and outcome among patients and rerun the full pipeline (filter, screening, score, split) 100 times. On shuffled data the held-out C-index is 0.50 ± 0.03; on the real data it is 0.64. This shows the pipeline does not create performance from nothing, which is something a train/test split alone cannot prove.
*Limit:* with 100 shuffles the smallest possible p-value is about 0.01.

## 6. Each gene counts in its own direction
The risk score is the mean of the z-scored expression of the 15 genes, each multiplied by +1 (risk) or −1 (protective), using training-set statistics.
*Why:* a plain average mixes genes that push survival in opposite directions, so their effects partly cancel. This is my own design choice; there is no specific paper behind it, but it is checked on synthetic data (the planted risk and protective genes are recovered with the right signs).

## 7. Primary tumours only, one row per patient
TCGA files also contain matched normal tissue and repeated samples. Mixing them in, or counting a patient several times, distorts both the expression analysis and the survival statistics (in an early version of this analysis that inflated a Cox model to 878 "observations" coming from only 475 patients, and averaged tumour and normal tissue for 59 of them).
*Why:* basic data hygiene (TCGA sample-type convention; see the note in references.md).

## 8. Fast univariate Cox screening
Ranking thousands of genes with a full Cox fit each is slow; a vectorised Newton–Raphson implementation (Breslow ties) fits all genes at once, which makes 30 splits × 200 random signatures affordable. It is checked against `lifelines` in the tests.

## 9. Check how stable the gene list is
I count how often each gene lands in the signature across the 30 training sets. Over 30 splits 106 different genes were chosen and none in all of them; about 100 genes pass a strict false-discovery cut-off in a typical training set. So the 15 genes are one sample of a larger group of correlated genes, and I give more weight to the biology than to the exact list.

## 10. Pathway analysis is interpretation, not validation
GSEA of the training-set ranking and an over-representation test of the 15 genes describe *what kind of biology* relates to survival. They are not evidence that the signature will predict well in new patients.
The ranking covers the ~1,900 most variable, well-expressed genes (the pool the signature was picked from). The over-representation test uses that same pool as its background; comparing 15 genes against the whole genome overstates significance (Wijesooriya et al. 2022). With the right background the best hit is weaker than it first looked (adjusted p = 0.009 instead of 0.002).

## 11. Reporting
I follow the spirit of transparent-reporting guidance (McShane & Hayes 2012; Collins et al. 2015): the cohort and inclusion criteria, the split, every result (including the unfavourable ones) and the limitations are stated. This repository is not a formal REMARK or TRIPOD submission.

## Not done (and why it matters)
- **External validation** on an independent cohort (for example from GEO). This is the most important missing piece.
- **More clinical and technical covariates:** smoking, driver mutations (EGFR, KRAS), tumour purity and lab batch.
- **Penalised models** (elastic-net Cox) or stability selection to get a more reliable gene list, and **time-dependent** performance curves.
- Tuning the signature size (15 genes was fixed in advance, not optimised).
