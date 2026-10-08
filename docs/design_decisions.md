# Design decisions

Why the analysis is built the way it is. Sources are in [references.md](references.md); where there is **no** literature behind a choice, it says so.

## 1. Everything data-dependent is fitted on training patients only
Gene filtering, screening, score normalisation, the risk threshold and the clinical baselines all use only the training patients; the held-out patients are touched once, at the end.
*Why:* if genes are chosen using the same patients that are later tested, the error looks smaller than it really is (Ambroise & McLachlan 2002). A test (`test_test_set_never_influences_the_fitted_signature`) corrupts the held-out patients' expression and checks that the signature does not change.

## 2. Many random splits, not one
A single 70/30 split is one draw; a lucky or unlucky draw can mislead. We repeat the whole procedure on 30 random splits and report the average and spread.
*Why:* in a re-analysis of seven published studies, gene lists were highly unstable across training sets and five of seven did not classify better than chance; the authors recommend validation by repeated random sampling (Michiels et al. 2005).
*Limit:* the 30 splits reuse the same 474 patients, so they are not independent evidence. An external cohort is still missing.

## 3. C-index as the headline metric, with a bootstrap interval
The C-index ranks pairs of patients by who died first, handles censored follow-up, and is easy to read (0.5 = coin flip, 1 = perfect). Intervals come from resampling the held-out patients.
*Why:* Harrell et al. (1996) describe it as an interpretable index of predictive discrimination for survival data and recommend validating it by bootstrap or cross-validation.

## 4. Compare against a clinical baseline *and* against random genes
A signature is only useful if it adds something to what is already known (age and stage), and if it is better than what chance gives. So we report (a) age + stage alone, (b) the signature alone, (c) both, and (d) the distribution of C-index values of 200 random 15-gene signatures per split, treated exactly like the real one (directions estimated on training patients, same normalisation).
*Why:* most random gene sets show some association with cancer outcome, largely through proliferation (Venet et al. 2011, in breast cancer).
*Result:* the signature adds ≈ +0.04 C-index to age and stage on average, but does not consistently beat random sets (above the 95th percentile of random sets in 3 of 30 splits). This is reported as it came out.

## 5. Each gene counts in its own direction
The risk score is the mean of the z-scored expression of the 15 genes, each multiplied by +1 (risk) or −1 (protective), using training-set statistics.
*Why:* a plain average mixes genes that push survival in opposite directions, so their effects partly cancel. This is our own design choice; there is no specific paper behind it, but it is checked on synthetic data (the planted risk and protective genes are recovered with the right signs).

## 6. Primary tumours only, one row per patient
TCGA files also contain matched normal tissue and repeated samples. Mixing them in, or counting a patient several times, distorts both the expression analysis and the survival statistics (the original notebook's Cox model counted 878 observations although the files hold only 507 patients with survival data).
*Why:* basic data hygiene (TCGA sample-type convention; see the note in references.md).

## 7. Fast univariate Cox screening
Ranking thousands of genes with a full Cox fit each is slow; a vectorised Newton–Raphson implementation (Breslow ties) fits all genes at once, which makes 30 splits × 200 random signatures affordable. It is checked against `lifelines` in the tests.

## 8. Pathway analysis is interpretation, not validation
GSEA of the training-set ranking and an over-representation test of the 15 genes describe *what kind of biology* relates to survival. They are not evidence that the signature will predict well in new patients.

## 9. Reporting
We follow the spirit of transparent-reporting guidance (McShane & Hayes 2012; Collins et al. 2015): the cohort and inclusion criteria, the split, every result (including the unfavourable comparison with random signatures) and the limitations are stated. This repository is not a formal REMARK or TRIPOD submission.

## Not done (and why it matters)
- **External validation** on an independent cohort (for example from GEO).
- **Adjusting for proliferation** (a "meta-proliferation" score) to see whether the signature carries anything beyond it (Venet et al. 2011).
- **Penalised models** (elastic-net Cox) and **time-dependent** performance curves.
- Tuning the signature size (15 genes is inherited, not optimised).
