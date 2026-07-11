"""
=============================================================
Shared helpers for the thesis_rsf kernel (Random Survival Forest + SHAP)
Student ID: 3068804
Supervisor: Dr Sudip Mondal
MSc Health Data Science, University of Birmingham
=============================================================

Both random_survival_forest.py and rsf_shap.py import from here so they encode the
feature set the *same* way. If they diverged, the SHAP explanation would describe a
different representation than the forest was trained on, and the two figures would
not agree. Keeping one encoder is the whole point of this module.

ENVIRONMENT: the thesis_rsf env (scikit-survival). Run the two scripts from the
project root (Thesis/) so the relative paths resolve.

FEATURE SET (expanded clinical, matches the clean CSVs written by kernel 1):
  age, gender, overall stage, T stage, N stage, M stage, prior malignancy,
  residual disease, treatment received, race, ethnicity.

Encoding: one integer column per variable (not one-hot), so the forest gives one
importance value per clinical variable and the SHAP beeswarm has one row per
variable. Ordered variables (stage, T, N, M, residual, treatment) get a meaningful
ordinal code; 'Unknown' is kept separable at -1 so the forest can split it out
rather than treating it as a real level. The genuinely nominal variables (race,
ethnicity, prior malignancy) are given fixed integer codes - a mild abuse for a
tree, acceptable at this low cardinality, and noted in the write-up.

Complete-case is on age + overall stage, exactly as the Cox model in kernel 1, so
n matches (426 / 153 / 393). Every other variable carries an 'Unknown' level, so no
further rows are dropped.
"""

import os
import numpy as np
import pandas as pd
from sksurv.util import Surv

# --- paths (relative, run from Thesis/ root) --------------------------------
CLEAN_DIR   = 'Clean'
FIGURE_DIR  = 'Figures'
RESULTS_DIR = 'Results'

# --- labels and house colours -----------------------------------------------
COLOURS = {'coad': '#1E6091', 'read': '#2A9D8F', 'stad': '#E63946'}
LABELS  = {'coad': 'COAD (Colon)', 'read': 'READ (Rectum)', 'stad': 'STAD (Stomach)'}

# --- forest / SHAP feature set ----------------------------------------------
FEATURES = ['age', 'gender', 'stage', 't_stage', 'n_stage', 'm_stage',
            'prior_malignancy', 'residual_disease', 'treatment_therapy',
            'race', 'ethnicity']

FEATURE_LABELS = {
    'age': 'Age', 'gender': 'Gender (Male)', 'stage': 'Overall stage',
    't_stage': 'T stage', 'n_stage': 'N stage', 'm_stage': 'M stage',
    'prior_malignancy': 'Prior malignancy', 'residual_disease': 'Residual disease',
    'treatment_therapy': 'Treatment received', 'race': 'Race', 'ethnicity': 'Ethnicity',
}

# integer encodings. Ordered levels get an ascending code; 'Unknown' -> -1 (kept
# separable). age and gender are handled directly in _encode_forest.
ENCODINGS = {
    'stage':             {'Stage I': 1, 'Stage II': 2, 'Stage III': 3, 'Stage IV': 4},
    't_stage':           {'T1': 1, 'T2': 2, 'T3': 3, 'T4': 4, 'Unknown': -1},
    'n_stage':           {'N0': 1, 'N1': 2, 'N2': 3, 'N3': 4, 'Unknown': -1},
    'm_stage':           {'M0': 1, 'M1': 2, 'Unknown': -1},
    'prior_malignancy':  {'No': 0, 'Yes': 1, 'Unknown': -1},
    'residual_disease':  {'R0': 0, 'R+': 1, 'Unknown': -1},
    'treatment_therapy': {'No': 0, 'Yes': 1, 'Unknown': -1},
    'race':              {'White': 0, 'Black': 1, 'Asian': 2, 'Other/Unknown': 3},
    'ethnicity':         {'Non-Hispanic': 0, 'Hispanic': 1, 'Unknown': -1},
}

# --- matched Cox design (mirrors kernel 1's survival_analysis.py) ------------
# Used only by the forest script, to score a Cox model on the same CV folds as the
# forest (a fair out-of-sample comparison). T/N/M are excluded (collinear with
# overall stage); prior_treatment is excluded (near-constant). One-hot, reference
# level dropped. Column names are internal only (never plotted here).
COX_CATEGORICALS = {
    'stage':             ('Stage I',      ['Stage II', 'Stage III', 'Stage IV']),
    'prior_malignancy':  ('No',           ['Yes', 'Unknown']),
    'residual_disease':  ('R0',           ['R+', 'Unknown']),
    'treatment_therapy': ('No',           ['Yes', 'Unknown']),
    'race':              ('White',        ['Black', 'Asian', 'Other/Unknown']),
    'ethnicity':         ('Non-Hispanic', ['Hispanic', 'Unknown']),
}


def _complete_case(df):
    """Drop rows missing overall stage or age (same rule as the Cox model)."""
    return df.dropna(subset=['stage', 'age']).copy()


def _encode_forest(d):
    """One integer column per feature, in FEATURES order."""
    X = pd.DataFrame(index=d.index)
    X['age']    = d['age'].astype(float)
    X['gender'] = (d['gender'] == 'Male').astype(int)
    for col, mapping in ENCODINGS.items():
        X[col] = d[col].map(mapping)
    X = X[FEATURES]

    if X.isna().any().any():   # a category outside the maps -> should not happen post-kernel-1
        offenders = {c: d[c][X[c].isna()].unique().tolist() for c in X.columns if X[c].isna().any()}
        print(f"  [warn] unmapped categories coded as -1: {offenders}")
        X = X.fillna(-1)
    return X


def _encode_cox(d):
    """One-hot Cox design (no duration/event columns; those live in y)."""
    X = pd.DataFrame(index=d.index)
    X['age']         = d['age'].astype(float)
    X['gender_male'] = (d['gender'] == 'Male').astype(int)
    for feat, (_ref, levels) in COX_CATEGORICALS.items():
        for lev in levels:
            name = lev if feat == 'stage' else f'{feat}={lev}'
            X[name] = (d[feat] == lev).astype(int)
    return X


def _make_y(d):
    """scikit-survival structured target (event bool + time)."""
    return Surv.from_arrays(event=(d['event'] == 1).values,
                            time=d['survival_time'].astype(float).values)


def build_xy(df):
    """Forest / SHAP design: (X features, y). Complete-case on age + stage."""
    d = _complete_case(df)
    return _encode_forest(d), _make_y(d)


def prepare_cohort(df):
    """Forest + matched-Cox designs on identical rows, plus the shared target.
    Returns (X_forest, X_cox, y), all aligned by position so the same CV fold
    indices apply to all three."""
    d = _complete_case(df)
    return _encode_forest(d), _encode_cox(d), _make_y(d)


def load_clean():
    """Read the three clean per-patient CSVs written by kernel 1."""
    frames = {}
    for key in ('coad', 'read', 'stad'):
        frames[key] = pd.read_csv(f'{CLEAN_DIR}/{key}_clean.csv')
    return frames


def read_cox_cindex():
    """cohort key -> Cox C-index, read from Results/cox_cindex_summary.csv (written
    by kernel 1). This is the lifelines apparent (in-sample) value reported in the
    main analysis; the forest figure shows it for reference. Returns {} if the file
    is not there yet, in which case the reference bar is simply omitted."""
    path = f'{RESULTS_DIR}/cox_cindex_summary.csv'
    if not os.path.exists(path):
        print(f"  [warn] {path} not found - run kernel 1 first. Cox reference bar omitted.")
        return {}
    dfc = pd.read_csv(path)
    label_to_key = {v: k for k, v in LABELS.items()}
    out = {}
    for _, r in dfc.iterrows():
        key = label_to_key.get(str(r['cohort']).strip())
        if key is not None:
            out[key] = float(r['Cox_C_index'])
    return out
