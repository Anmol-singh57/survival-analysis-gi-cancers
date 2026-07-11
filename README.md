# Survival Analysis of Colorectal and Stomach Cancers

**A comparative study using TCGA clinical data.**

![Python](https://img.shields.io/badge/Python-3.10-blue)
![License](https://img.shields.io/badge/License-MIT-green)
![Status](https://img.shields.io/badge/MSc%20Dissertation-2026-lightgrey)

MSc Health Data Science dissertation project. A full survival analysis pipeline (Kaplan-Meier, penalised Cox regression, random survival forest with SHAP, and cross-cancer transfer learning) applied to three gastrointestinal cancer cohorts in The Cancer Genome Atlas (TCGA): colon (COAD), rectum (READ) and stomach (STAD). The central question is not only what predicts survival in each cancer, but which predictors behave the same way across cancers, and whether a model trained on one cancer transfers to another.

## Key results

**Cohorts** (after deduplication to one record per patient and exclusion of records with missing survival time):

| Cohort | Patients | Deaths | Death rate | Median KM survival |
|---|---|---|---|---|
| COAD (colon) | 437 | 95 | 21.7% | 2,821 days |
| READ (rectum) | 161 | 27 | 16.8% | 1,741 days |
| STAD (stomach) | 412 | 168 | 40.8% | 881 days |

**Discrimination** (concordance index):

| Model | COAD | READ | STAD |
|---|---|---|---|
| Penalised Cox (in-sample) | 0.76 | 0.83 | 0.70 |
| Random survival forest (5-fold CV) | 0.76 | 0.66 | 0.65 |

**Headline findings**

- **Stage is the dominant predictor** of survival across all three cancers, at every level of the pipeline (log-rank, Cox, permutation importance, SHAP). Residual disease after surgery is the second most consistent signal. Demographic variables carry no reliable prognostic signal after adjustment.
- **The prognostic structure is shared.** Testing each predictor for a difference in effect between colorectal and stomach, all four covariate-by-cancer interaction tests are non-significant (age p = 0.46, residual disease p = 0.48, sex p = 0.81, stage p = 0.10).
- **Cross-cancer transfer works.** A random survival forest trained only on colorectal patients, and never shown a stomach patient, scores a concordance index of **0.646** on stomach, against a within-stomach cross-validated ceiling of **0.650**. Splitting stomach patients by their colorectal-predicted risk separates survival sharply (log-rank p = 1.5e-7; integrated Brier score 0.192; mean time-dependent AUC 0.671).
- **What does not transfer is absolute risk.** After adjusting for age, stage and sex, colon has a hazard ratio of **0.36** versus stomach (95% CI 0.27 to 0.47), so stomach is roughly three times worse stage-for-stage. Its poor prognosis is intrinsic to the disease, not a consequence of later-stage diagnosis.

## Methods

Each cancer is analysed independently before any cross-cancer comparison is made. The pipeline runs in the following order:

1. **Exploratory data analysis** of observed times, age, sex and stage distributions.
2. **Kaplan-Meier estimation** with log-rank tests, stratified by stage, sex and age group.
3. **Penalised (ridge) Cox proportional hazards regression** on an eight-variable clinical subset.
4. **Random survival forest** on the full eleven-variable set (including T, N and M stage), evaluated by five-fold cross-validation, with **permutation importance** and **SHAP** explanations.
5. **Cross-cancer analyses**: an adjusted between-cancer comparison, transfer of a colorectal-trained forest to stomach, and an invariance test of which predictors keep the same effect across cancers.

## Repository structure

```
.
├── README.md
├── LICENSE
├── .gitignore
├── environment-thesis.yml         # conda env for KM + Cox (lifelines, scipy<1.12)
├── environment-thesis_rsf.yml     # conda env for RSF + SHAP + transfer (scikit-survival, shap)
├── Clean/                         # cleaned, analysis-ready data and preprocessing
├── Figures/                       # generated figures (all figures in the report)
├── Results/                       # metrics, tables and model outputs (CSV)
├── report/
│   ├── Survival_Analysis_of_Colorectal_and_Stomach_Cancers.pdf   # dissertation
│   └── Survival_Analysis_of_Colorectal_slides.pptx               # presentation
└── src/  (or notebooks/)          # analysis scripts / notebooks
```

> Adjust the `src/` names to match your actual scripts (for example `eda_vlab.py` and your KM/Cox/RSF/transfer files). Keep `Figures/` and `Results/` tracked so every figure and statistic in the report can be traced back to an output.

## Data

All clinical data come from TCGA via the NCI Genomic Data Commons (GDC) portal, `portal.gdc.cancer.gov`. Only **clinical** data are used (not genomic or imaging), and TCGA clinical data are **open-access**.

Raw GDC downloads (`.tsv`) are not redistributed in this repository (see `.gitignore`); to reproduce from scratch, download the clinical files for projects **TCGA-COAD**, **TCGA-READ** and **TCGA-STAD** from the GDC portal and place them where the preprocessing expects them. The cleaned, analysis-ready data and all results are included so the analysis can be followed and re-run.

Note on preprocessing: the GDC clinical file stores one row per diagnosis-and-follow-up combination, so the pipeline first collapses to one record per patient (keyed on the case identifier) before any analysis. This step is essential; treating rows as patients inflates each cohort roughly threefold.

## Environment setup

Two conda environments are used because the Cox and forest stacks have an incompatible dependency (`lifelines` requires `scipy<1.12`, which conflicts with `scikit-survival`).

```bash
# KM and Cox
conda env create -f environment-thesis.yml
conda activate thesis

# Random survival forest, SHAP and transfer
conda env create -f environment-thesis_rsf.yml
conda activate thesis_rsf
```

> The environment files capture the key pins. To reproduce exact versions from your own machine, run `conda env export --no-builds > environment-thesis.yml` (and the same for `thesis_rsf`) and commit the result.

## Reproducing the analysis

Run the stages in the order listed under [Methods](#methods). Use the `thesis` environment for exploratory analysis, Kaplan-Meier and Cox; switch to `thesis_rsf` for the random survival forest, SHAP and the cross-cancer transfer. Figures are written to `Figures/` and metrics to `Results/`; both match the values reported in the dissertation.

## Report and presentation

The full dissertation and the presentation slides are in [`report/`](report/).

## Author

**Anmol Singh** — MSc Health Data Science, University of Birmingham
Supervisor: Dr Sudip Mondal

## License

Code in this repository is released under the MIT License (see [LICENSE](LICENSE)). The dissertation text and figures in `report/` are the author's academic work.
