# Survival Analysis of Colorectal and Stomach Cancers

A comparative study using TCGA clinical data.

![Python](https://img.shields.io/badge/python-3.10-blue) ![License](https://img.shields.io/badge/license-MIT-green)

MSc Health Data Science dissertation project. A full survival analysis pipeline (Kaplan-Meier, penalised Cox regression, random survival forest with SHAP, and cross-cancer transfer learning) applied to three gastrointestinal cancer cohorts in The Cancer Genome Atlas (TCGA): colon (COAD), rectum (READ) and stomach (STAD). The central question is not only what predicts survival in each cancer, but which predictors behave the same way across cancers, and whether a model trained on one cancer transfers to another, in both directions.

## Key results

Cohorts (after deduplication to one record per patient and exclusion of records with missing survival time):

| Cohort | Patients | Deaths | Death rate | Median KM survival |
| --- | --- | --- | --- | --- |
| COAD (colon) | 437 | 95 | 21.7% | 2,821 days |
| READ (rectum) | 161 | 27 | 16.8% | 1,741 days |
| STAD (stomach) | 412 | 168 | 40.8% | 881 days |

Discrimination (concordance index):

| Model | COAD | READ | STAD |
| --- | --- | --- | --- |
| Penalised Cox (in-sample) | 0.76 | 0.83 | 0.70 |
| Random survival forest (5-fold CV) | 0.76 | 0.66 | 0.65 |

### Headline findings

- **Stage is the dominant predictor** of survival across all three cancers, at every level of the pipeline (log-rank, Cox, permutation importance, SHAP). Residual disease after surgery is the second most consistent signal. Demographic variables carry no reliable prognostic signal after adjustment.
- **The prognostic structure is shared.** Testing each predictor for a difference in effect between colorectal and stomach, all four covariate-by-cancer interaction tests are non-significant (age p = 0.46, residual disease p = 0.48, sex p = 0.81, stage p = 0.10).
- **Cross-cancer transfer works, in both directions.** A random survival forest trained only on colorectal patients, and never shown a stomach patient, scores a concordance index of 0.646 on stomach against a within-stomach ceiling of 0.648; run in reverse, a stomach-trained forest scores 0.763 on colorectal against a within-colorectal ceiling of 0.760. Transfer therefore reaches the within-cohort ceiling in each direction. A Cox model transferred in place of the forest performs almost identically (0.649 against 0.646), so the result is not specific to the forest. Splitting stomach patients by their colorectal-predicted risk separates survival sharply (log-rank p = 1.5e-7; integrated Brier score 0.192; mean time-dependent AUC 0.671, corroborated by time-dependent ROC-AUCs of 0.658, 0.674 and 0.678 at one, two and three years).
- **What does not transfer is absolute risk.** After adjusting for age, stage and sex, colon has a hazard ratio of 0.36 versus stomach (95% CI 0.27 to 0.47), so stomach is roughly three times worse stage-for-stage. Its poor prognosis is intrinsic to the disease, not a consequence of later-stage diagnosis. This is why discrimination transfers but calibration does not: the risk ordering is shared across cancers, but the absolute level is not. The concordance index is only comparable within a target cohort, not between cohorts, so each transfer direction is judged against its own within-cohort ceiling.

## Methods

Each cancer is analysed independently before any cross-cancer comparison is made. The pipeline runs in the following order:

1. Exploratory data analysis of observed times, age, sex and stage distributions.
2. Kaplan-Meier estimation with log-rank tests, stratified by stage, sex and age group.
3. Penalised (ridge) Cox proportional hazards regression on an eight-variable clinical subset.
4. Random survival forest on the full eleven-variable set (including T, N and M stage), evaluated by five-fold cross-validation, with permutation importance and SHAP explanations.
5. Cross-cancer analyses: an adjusted between-cancer comparison; transfer of a colorectal-trained forest to stomach, with ROC-AUC, a method comparison, and the reverse direction against each cohort's own ceiling; and an invariance test of which predictors keep the same effect across cancers.

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
├── notebooks/
│   ├── Thesis.ipynb               # EDA, Kaplan-Meier, Cox (thesis env)
│   ├── Thesis.html                # rendered run of Thesis.ipynb
│   ├── thesis_rsf.ipynb           # random survival forest, SHAP, transfer (thesis_rsf env)
│   ├── thesis_rsf.html            # rendered run of thesis_rsf.ipynb
│   └── rsf_common.py              # shared feature definitions and helpers, imported by thesis_rsf.ipynb
└── report/
    ├── Survival_Analysis_of_Colorectal_and_Stomach_Cancers.pdf   # dissertation
    ├── Survival_Analysis_of_Colorectal_slides.pptx               # presentation
    └── Survival_Analysis_Supplementary_Slides.pptx               # supplementary slides (backup detail)
```

`Figures/` and `Results/` are kept under version control so every figure and statistic in the report can be traced back to an output.

## Data

All clinical data come from TCGA via the NCI Genomic Data Commons (GDC) portal, portal.gdc.cancer.gov. Only clinical data are used (not genomic or imaging), and TCGA clinical data are open-access.

Raw GDC downloads (.tsv) are not redistributed in this repository (see `.gitignore`); to reproduce from scratch, download the clinical files for projects TCGA-COAD, TCGA-READ and TCGA-STAD from the GDC portal and place them where the preprocessing expects them. The cleaned, analysis-ready data and all results are included so the analysis can be followed and re-run.

**Note on preprocessing:** the GDC clinical file stores one row per diagnosis-and-follow-up combination, so the pipeline first collapses to one record per patient (keyed on the case identifier) before any analysis. This step is essential; treating rows as patients inflates each cohort roughly threefold.

## Environment setup

Two conda environments are used because the Cox and forest stacks have an incompatible dependency (lifelines requires scipy<1.12, which conflicts with scikit-survival).

```bash
# KM and Cox
conda env create -f environment-thesis.yml
conda activate thesis

# Random survival forest, SHAP and transfer
conda env create -f environment-thesis_rsf.yml
conda activate thesis_rsf
```

## Reproducing the analysis

Run the two notebooks in `notebooks/` in order. `Thesis.ipynb` (exploratory analysis, Kaplan-Meier, Cox) runs in the `thesis` environment and writes the cleaned per-patient data and the Cox outputs. `thesis_rsf.ipynb` (random survival forest, SHAP, cross-cancer transfer) runs in the `thesis_rsf` environment and reads those cleaned files, so run `Thesis.ipynb` first. `thesis_rsf.ipynb` imports shared feature definitions and helpers from `rsf_common.py`, which must sit alongside it. Figures are written to `Figures/` and metrics to `Results/`; both match the values reported in the dissertation. Rendered HTML of a full run of each notebook is included next to it in `notebooks/`.

## Report and presentation

The dissertation, the presentation slides and a supplementary slide deck (backup detail for questions) are in `report/`.

## Author

Anmol Singh — MSc Health Data Science, University of Birmingham
Supervisor: Dr Sudip Mondal
