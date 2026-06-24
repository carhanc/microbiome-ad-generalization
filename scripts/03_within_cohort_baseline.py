import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import (GridSearchCV, StratifiedKFold,
                                      StratifiedGroupKFold, cross_val_predict)
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import LabelEncoder

try:
    import lightgbm as lgb
    HAS_LIGHTGBM = True
except ImportError:
    HAS_LIGHTGBM = False
    warnings.warn(
        "LightGBM not installed. Run: conda install -n microbiome-ad lightgbm -c conda-forge\n"
        "LightGBM results will be skipped.",
        stacklevel=2
    )

PROJECT_ROOT = Path(__file__).resolve().parent.parent
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
RESULTS_DIR   = PROJECT_ROOT / "results"
TABLES_DIR    = RESULTS_DIR / "tables"
MODEL_DIR     = RESULTS_DIR / "model_outputs" / "within_cohort_cv"

TABLES_DIR.mkdir(parents=True, exist_ok=True)
MODEL_DIR.mkdir(parents=True, exist_ok=True)

SUPERVISED_COHORTS = ["zhuang2018", "ling2020", "shanghai2022", "kazakhstan2022"]

RANDOM_STATE = 42

OUTER_FOLDS = 10
INNER_FOLDS = 5
N_BOOTSTRAP = 1000


def bootstrap_auc_ci(y_true: np.ndarray, y_score: np.ndarray,
                     n_boot: int = N_BOOTSTRAP,
                     alpha: float = 0.05,
                     rng: np.random.Generator = None) -> tuple[float, float]:
    if rng is None:
        rng = np.random.default_rng(RANDOM_STATE)
    n = len(y_true)
    auc_boot = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        y_t = y_true[idx]
        y_s = y_score[idx]
        if len(np.unique(y_t)) < 2:
            continue
        auc_boot.append(roc_auc_score(y_t, y_s))
    auc_boot = np.array(auc_boot)
    pct_lo = 100 * alpha / 2
    pct_hi = 100 * (1 - alpha / 2)
    return float(np.percentile(auc_boot, pct_lo)), \
           float(np.percentile(auc_boot, pct_hi))


def youden_threshold_metrics(y_true: np.ndarray,
                              y_score: np.ndarray) -> dict:
    from sklearn.metrics import roc_curve
    fpr, tpr, thresholds = roc_curve(y_true, y_score)
    j_scores = tpr - fpr
    best_idx  = np.argmax(j_scores)
    best_thr  = thresholds[best_idx]
    sensitivity = float(tpr[best_idx])
    specificity = float(1 - fpr[best_idx])
    return {
        "threshold": float(best_thr),
        "sensitivity": round(sensitivity, 4),
        "specificity": round(specificity, 4),
        "youden_j": round(j_scores[best_idx], 4),
    }


def make_logreg_pipeline():
    model = LogisticRegression(
        solver="saga",
        l1_ratio=0.5,
        max_iter=10000,
        random_state=RANDOM_STATE,
    )
    param_grid = {"C": [0.001, 0.01, 0.1, 1.0, 10.0]}
    return model, param_grid


def make_lgbm_pipeline():
    model = lgb.LGBMClassifier(
        objective="binary",
        n_jobs=4,
        random_state=RANDOM_STATE,
        verbose=-1,
    )
    param_grid = {
        "n_estimators": [50, 100, 200],
        "learning_rate": [0.05, 0.1],
        "max_depth":     [3, 5],
        "num_leaves":    [15, 31],
    }
    return model, param_grid


def run_nested_cv(X: np.ndarray, y: np.ndarray,
                  model, param_grid: dict,
                  model_name: str,
                  cohort: str,
                  outer_folds: int = OUTER_FOLDS,
                  inner_folds: int = INNER_FOLDS,
                  groups: np.ndarray = None) -> dict:
    if groups is not None:
        outer_cv = StratifiedGroupKFold(n_splits=outer_folds)
    else:
        outer_cv = StratifiedKFold(n_splits=outer_folds, shuffle=True,
                                   random_state=RANDOM_STATE)
    inner_cv = StratifiedKFold(n_splits=inner_folds, shuffle=True,
                               random_state=RANDOM_STATE)

    oof_proba  = np.zeros(len(y))
    fold_aucs  = []

    split_iter = (outer_cv.split(X, y, groups=groups) if groups is not None
                  else outer_cv.split(X, y))

    for fold_idx, (train_idx, test_idx) in enumerate(split_iter):
        X_train, X_test = X[train_idx], X[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]

        gs = GridSearchCV(
            estimator=model,
            param_grid=param_grid,
            cv=inner_cv,
            scoring="roc_auc",
            n_jobs=-1,
            refit=True,
        )
        gs.fit(X_train, y_train)
        best_model = gs.best_estimator_

        fold_proba = best_model.predict_proba(X_test)[:, 1]
        oof_proba[test_idx] = fold_proba

        fold_auc = roc_auc_score(y_test, fold_proba)
        fold_aucs.append(fold_auc)
        print(f"    fold {fold_idx+1:2d}/{outer_folds}: AUC={fold_auc:.4f}  "
              f"best_params={gs.best_params_}")

    oof_auc = roc_auc_score(y, oof_proba)

    rng = np.random.default_rng(RANDOM_STATE)
    ci_lo, ci_hi = bootstrap_auc_ci(y, oof_proba, rng=rng)

    thresh_metrics = youden_threshold_metrics(y, oof_proba)

    print(f"    pooled OOF AUC = {oof_auc:.4f}  "
          f"95% CI [{ci_lo:.4f}, {ci_hi:.4f}]  "
          f"sens={thresh_metrics['sensitivity']:.3f}  "
          f"spec={thresh_metrics['specificity']:.3f}")

    oof_df = pd.DataFrame({
        "y_true":  y,
        "y_score": oof_proba,
    })
    oof_path = MODEL_DIR / f"{cohort}_{model_name}_oof_predictions.csv"
    oof_df.to_csv(oof_path, index=False)

    return {
        "cohort":       cohort,
        "model":        model_name,
        "n_samples":    len(y),
        "n_AD":         int(y.sum()),
        "n_CN":         int((y == 0).sum()),
        "auc_oof":      round(oof_auc, 4),
        "auc_ci_lo":    round(ci_lo, 4),
        "auc_ci_hi":    round(ci_hi, 4),
        "auc_mean_folds": round(float(np.mean(fold_aucs)), 4),
        "auc_sd_folds":   round(float(np.std(fold_aucs)), 4),
        "sensitivity":  thresh_metrics["sensitivity"],
        "specificity":  thresh_metrics["specificity"],
        "youden_j":     thresh_metrics["youden_j"],
        "fold_aucs":    [round(a, 4) for a in fold_aucs],
    }


def main():
    print("03_within_cohort_baseline.py — nested CV baselines")

    clr_path = PROCESSED_DIR / "clr_matrix.csv"
    print(f"\nLoading {clr_path}")
    clr_df = pd.read_csv(clr_path, index_col=0)
    print(f"  {clr_df.shape[0]} samples × {clr_df.shape[1]} columns")

    genus_cols = [c for c in clr_df.columns if c != "cohort"]
    features = genus_cols  # same list, used below in loops
    print(f"  {len(features)} genus features")

    diag_path = PROCESSED_DIR / "diagnosis_labels.csv"
    print(f"\nLoading {diag_path}")
    diag_df = pd.read_csv(diag_path)

    clr_reset = clr_df.reset_index().rename(columns={"sample_id": "run_id"})
    merged = clr_reset.merge(diag_df[["run_id", "cohort", "diagnosis"]],
                             on="run_id", suffixes=("_clr", ""))
    if "cohort_clr" in merged.columns:
        merged = merged.drop(columns=["cohort_clr"])
    merged = merged.set_index("run_id")

    print(f"  merged: {merged.shape}")
    print("  diagnosis counts:")
    print(merged.groupby(["cohort", "diagnosis"]).size().to_string())

    all_results = []

    for cohort in SUPERVISED_COHORTS:
        print(f"\n{cohort}")

        cohort_data = merged[merged["cohort"] == cohort].copy()

        if "MCI" in cohort_data["diagnosis"].unique():
            mci_n = (cohort_data["diagnosis"] == "MCI").sum()
            print(f"  {mci_n} MCI samples — binary analysis uses AD vs CN only")
            cohort_binary = cohort_data[cohort_data["diagnosis"].isin(["AD", "CN"])].copy()
        else:
            cohort_binary = cohort_data.copy()

        cohort_binary = cohort_binary[cohort_binary["diagnosis"].notna()]

        n_ad = (cohort_binary["diagnosis"] == "AD").sum()
        n_cn = (cohort_binary["diagnosis"] == "CN").sum()
        print(f"  AD={n_ad}  CN={n_cn}  total={len(cohort_binary)}")

        if n_ad < 5 or n_cn < 5:
            print(f"  too few samples per class for {OUTER_FOLDS}-fold CV, skipping")
            continue

        X = cohort_binary[genus_cols].to_numpy()
        y = (cohort_binary["diagnosis"] == "AD").astype(int).to_numpy()

        groups = None

        print("\n  logistic regression (elastic-net)")
        logreg_model, logreg_grid = make_logreg_pipeline()
        logreg_result = run_nested_cv(
            X, y, logreg_model, logreg_grid,
            model_name="logreg",
            cohort=cohort,
            groups=groups,
        )
        all_results.append(logreg_result)

        if HAS_LIGHTGBM:
            print("\n  LightGBM")
            lgbm_model, lgbm_grid = make_lgbm_pipeline()
            lgbm_result = run_nested_cv(
                X, y, lgbm_model, lgbm_grid,
                model_name="lgbm",
                cohort=cohort,
                groups=groups,
            )
            all_results.append(lgbm_result)

    print("\nshanghai2022 — 3-class secondary (AD / MCI / CN)")

    sh_data = merged[merged["cohort"] == "shanghai2022"].copy()
    sh_data = sh_data[sh_data["diagnosis"].notna()]
    X_sh = sh_data[genus_cols].to_numpy()

    le = LabelEncoder()
    y_sh = le.fit_transform(sh_data["diagnosis"])
    class_names = le.classes_
    print(f"  classes: {class_names}  encoded: {dict(zip(class_names, range(len(class_names))))}")

    outer_cv_3 = StratifiedKFold(n_splits=OUTER_FOLDS, shuffle=True, random_state=RANDOM_STATE)
    inner_cv_3 = StratifiedKFold(n_splits=INNER_FOLDS, shuffle=True, random_state=RANDOM_STATE)

    oof_proba_3 = np.zeros((len(y_sh), len(class_names)))
    fold_aucs_3 = []

    for fold_idx, (train_idx, test_idx) in enumerate(outer_cv_3.split(X_sh, y_sh)):
        X_tr, X_te = X_sh[train_idx], X_sh[test_idx]
        y_tr, y_te = y_sh[train_idx], y_sh[test_idx]

        # sklearn 1.7 dropped multi_class — saga figures it out
        gs3 = GridSearchCV(
            estimator=LogisticRegression(
                solver="saga",
                l1_ratio=0.5, max_iter=10000,
                random_state=RANDOM_STATE
            ),
            param_grid={"C": [0.001, 0.01, 0.1, 1.0, 10.0]},
            cv=inner_cv_3,
            scoring="roc_auc_ovr_weighted",
            n_jobs=-1,
        )
        gs3.fit(X_tr, y_tr)
        fold_proba = gs3.predict_proba(X_te)
        oof_proba_3[test_idx] = fold_proba
        fold_auc_macro = roc_auc_score(y_te, fold_proba, multi_class="ovr",
                                        average="macro")
        fold_aucs_3.append(fold_auc_macro)
        print(f"    fold {fold_idx+1:2d}/{OUTER_FOLDS}: macro-AUC={fold_auc_macro:.4f}  "
              f"best_C={gs3.best_params_['C']}")

    oof_auc_3 = roc_auc_score(y_sh, oof_proba_3, multi_class="ovr", average="macro")
    print(f"    pooled OOF macro-AUC = {oof_auc_3:.4f}")

    all_results.append({
        "cohort":       "shanghai2022_3class",
        "model":        "logreg",
        "n_samples":    len(y_sh),
        "n_AD":         int((sh_data["diagnosis"] == "AD").sum()),
        "n_CN":         int((sh_data["diagnosis"] == "CN").sum()),
        "n_MCI":        int((sh_data["diagnosis"] == "MCI").sum()),
        "auc_oof":      round(oof_auc_3, 4),
        "auc_ci_lo":    None,
        "auc_ci_hi":    None,
        "auc_mean_folds": round(float(np.mean(fold_aucs_3)), 4),
        "auc_sd_folds":   round(float(np.std(fold_aucs_3)), 4),
        "sensitivity":  None,
        "specificity":  None,
        "youden_j":     None,
        "fold_aucs":    [round(a, 4) for a in fold_aucs_3],
        "note":         "3-class macro OVR AUC; secondary analysis",
    })

    results_df = pd.DataFrame(all_results)
    out_path = TABLES_DIR / "within_cohort_auc.csv"
    results_df.to_csv(out_path, index=False)

    print("\nResults summary")
    primary = results_df[~results_df["cohort"].str.endswith("3class")].copy()
    print(primary[["cohort", "model", "n_samples",
                   "auc_oof", "auc_ci_lo", "auc_ci_hi",
                   "sensitivity", "specificity"]].to_string(index=False))

    print(f"\n  {out_path}")
    print(f"  OOF predictions in {MODEL_DIR}/")
    print("\nFinished. Next: python scripts/04_cross_cohort_generalization.py")


if __name__ == "__main__":
    main()
