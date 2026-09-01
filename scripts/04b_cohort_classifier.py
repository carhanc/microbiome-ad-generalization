from pathlib import Path
import warnings

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GridSearchCV, StratifiedKFold
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import LabelEncoder, LabelBinarizer

PROJECT_ROOT  = Path(__file__).resolve().parent.parent
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
TABLES_DIR    = PROJECT_ROOT / "results" / "tables"
TABLES_DIR.mkdir(parents=True, exist_ok=True)

RANDOM_STATE = 42
OUTER_FOLDS  = 5
INNER_FOLDS  = 5

LOGREG_GRID = {"C": [0.001, 0.01, 0.1, 1.0, 10.0]}


def main():
    print("04b_cohort_classifier.py — cohort identity classifier")
    print("Quantifies residual cohort-specific signal in CLR-transformed data,\n"
          "complementing the PERMANOVA cohort-variance result.")

    clr_path = PROCESSED_DIR / "clr_matrix.csv"
    print(f"\nLoading {clr_path}")
    clr_df = pd.read_csv(clr_path, index_col=0)
    genus_cols = [c for c in clr_df.columns if c != "cohort"]
    X = clr_df[genus_cols].to_numpy()

    le = LabelEncoder()
    y = le.fit_transform(clr_df["cohort"])
    cohort_names = le.classes_
    n_classes = len(cohort_names)

    print(f"  {X.shape[0]} samples x {X.shape[1]} genera, {n_classes} cohorts")
    print("  class sizes:", dict(zip(cohort_names, np.bincount(y))))

    outer_cv = StratifiedKFold(n_splits=OUTER_FOLDS, shuffle=True,
                               random_state=RANDOM_STATE)
    inner_cv = StratifiedKFold(n_splits=INNER_FOLDS, shuffle=True,
                               random_state=RANDOM_STATE)

    oof_proba = np.zeros((len(y), n_classes))

    print(f"\nRunning {OUTER_FOLDS}-fold CV (5-class cohort prediction, "
          f"{INNER_FOLDS}-fold inner CV for C)...")
    for fold_idx, (train_idx, test_idx) in enumerate(outer_cv.split(X, y)):
        X_train, X_test = X[train_idx], X[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]

        base = LogisticRegression(
            solver="saga",
            l1_ratio=0.5,
            max_iter=10000,
            random_state=RANDOM_STATE,
        )
        gs = GridSearchCV(
            estimator=base,
            param_grid=LOGREG_GRID,
            cv=inner_cv,
            scoring="roc_auc_ovr_weighted",
            n_jobs=-1,
            refit=True,
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            gs.fit(X_train, y_train)

        fold_proba = gs.predict_proba(X_test)
        oof_proba[test_idx] = fold_proba

        fold_auc_macro = roc_auc_score(y_test, fold_proba, multi_class="ovr",
                                        average="macro")
        print(f"    fold {fold_idx+1:2d}/{OUTER_FOLDS}: macro-AUC={fold_auc_macro:.4f}  "
              f"best_C={gs.best_params_['C']}")

    macro_auc = roc_auc_score(y, oof_proba, multi_class="ovr", average="macro")
    print(f"\n  Pooled OOF macro-AUC (cohort identity, 5-class OvR) = {macro_auc:.4f}")

    lb = LabelBinarizer()
    y_bin = lb.fit_transform(y)  # (n_samples, n_classes) one-hot, class order = sorted unique y

    print("\n  Per-cohort one-vs-rest AUC:")
    rows = []
    for i, cohort in enumerate(cohort_names):
        auc_i = roc_auc_score(y_bin[:, i], oof_proba[:, i])
        n_i = int((y == i).sum())
        print(f"    {cohort:20s}  n={n_i:3d}  AUC={auc_i:.4f}")
        rows.append({
            "cohort":    cohort,
            "n_samples": n_i,
            "auc_ovr":   round(auc_i, 4),
        })

    rows.append({
        "cohort":    "macro_average",
        "n_samples": len(y),
        "auc_ovr":   round(macro_auc, 4),
    })

    result_df = pd.DataFrame(rows)
    out_path = TABLES_DIR / "cohort_identity_classifier.csv"
    result_df.to_csv(out_path, index=False)
    print(f"\n  wrote {out_path}")
    print(result_df.to_string(index=False))

    print("\nDone.")
    print(f"  {out_path}")


if __name__ == "__main__":
    main()
