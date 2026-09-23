"""Per-sample OOF/LOCO predictions under strict training-only feature selection
(Reviewer 5, Point 5 -- promoting training-only selection to the PRIMARY
classifier-performance pipeline).

09_training_only_feature_sensitivity.py already established the strict
training-only AUC summary numbers (within-cohort, LOCO, pairwise) used as the
sensitivity analysis in the Reviewer 3/4 revision. This script reuses that
exact same model-fitting / feature-selection code, unmodified, but additionally
saves per-sample (run_id-tagged) OOF/LOCO predicted probabilities -- which
09_training_only_feature_sensitivity.py did not persist -- so that the formal
paired diagnosis-stratified bootstrap (13_auc_difference_bootstrap.py's
procedure) can be re-run on the training-only pipeline's own predictions
(17_training_only_auc_difference_bootstrap.py), rather than on the old
global-396 predictions.

Sanity check: the within-cohort/LOCO AUCs recomputed here from these saved
per-sample predictions are asserted to match 09's already-published summary
CSVs (results/tables/within_cohort_auc_training_only.csv,
loco_auc_training_only.csv) to full rounding precision, confirming this is the
same computation, not a new/different one.
"""
import importlib.util
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GridSearchCV, StratifiedKFold
from sklearn.metrics import roc_auc_score
import lightgbm as lgb

PROJECT_ROOT = Path(__file__).resolve().parent.parent
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
TABLES_DIR = PROJECT_ROOT / "results" / "tables"
OUT_DIR_WC = PROJECT_ROOT / "results" / "model_outputs" / "within_cohort_cv_training_only"
OUT_DIR_CC = PROJECT_ROOT / "results" / "model_outputs" / "cross_cohort_training_only"
OUT_DIR_WC.mkdir(parents=True, exist_ok=True)
OUT_DIR_CC.mkdir(parents=True, exist_ok=True)


def _load_module(name, filename):
    spec = importlib.util.spec_from_file_location(name, PROJECT_ROOT / "scripts" / filename)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


harmonize = _load_module("harmonize_16", "01_harmonize_taxonomy.py")
clr_mod = _load_module("clr_transform_16", "02_clr_transform.py")

SUPERVISED_COHORTS = ["zhuang2018", "ling2020", "shanghai2022", "kazakhstan2022"]
RANDOM_STATE = 42
OUTER_FOLDS = 10
INNER_FOLDS = 5
MIN_PREVALENCE = 0.20

LOGREG_GRID = {"C": [0.001, 0.01, 0.1, 1.0, 10.0]}
LGBM_GRID = {
    "n_estimators": [50, 100, 200],
    "learning_rate": [0.05, 0.1],
    "max_depth": [3, 5],
    "num_leaves": [15, 31],
}


def load_raw_tables():
    tables = {}
    for cohort in SUPERVISED_COHORTS:
        df = harmonize.load_genus_table(cohort)
        df = harmonize.apply_synonym_map(df)
        df = harmonize.drop_uninformative_genera(df)
        tables[cohort] = df
    return tables


def load_labels():
    return pd.read_csv(PROCESSED_DIR / "diagnosis_labels.csv")


def labeled_ids_for(diag_df, cohort):
    sub = diag_df[(diag_df["cohort"] == cohort) & (diag_df["diagnosis"].isin(["AD", "CN"]))]
    return sub.set_index("run_id")["diagnosis"]


def select_genera_training_only(raw_tables, sample_ids_by_cohort):
    restricted = {}
    for cohort, ids in sample_ids_by_cohort.items():
        df = raw_tables[cohort]
        restricted[cohort] = df.loc[df.index.intersection(ids)]
    return harmonize.apply_prevalence_filter(restricted, min_prevalence=MIN_PREVALENCE, min_cohorts=1)


def project_and_clr(raw_tables, cohort, sample_ids, retained_genera_sorted):
    df = raw_tables[cohort]
    sub = df.reindex(index=sample_ids, columns=retained_genera_sorted, fill_value=0)
    counts = sub.to_numpy(dtype=float)
    counts_replaced = clr_mod.multiplicative_replacement(counts, delta=0.65)
    return clr_mod.clr_transform(counts_replaced)


def tune_and_fit(X_train, y_train, model_name):
    inner_cv = StratifiedKFold(n_splits=INNER_FOLDS, shuffle=True, random_state=RANDOM_STATE)
    if model_name == "logreg":
        base = LogisticRegression(solver="saga", l1_ratio=0.5, max_iter=10000, random_state=RANDOM_STATE)
        grid = LOGREG_GRID
    else:
        base = lgb.LGBMClassifier(objective="binary", n_jobs=4, random_state=RANDOM_STATE, verbose=-1)
        grid = LGBM_GRID
    gs = GridSearchCV(base, grid, cv=inner_cv, scoring="roc_auc", n_jobs=-1, refit=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        gs.fit(X_train, y_train)
    return gs.best_estimator_, gs.best_params_


def run_within_cohort(raw_tables, diag_df):
    print("\n=== Within-cohort OOF predictions, training-only feature selection ===")
    summary_rows = []
    for cohort in SUPERVISED_COHORTS:
        labels = labeled_ids_for(diag_df, cohort)
        ids = np.array(labels.index)
        y = (labels.values == "AD").astype(int)
        outer_cv = StratifiedKFold(n_splits=OUTER_FOLDS, shuffle=True, random_state=RANDOM_STATE)

        for model_name in ["logreg", "lgbm"]:
            oof_proba = np.zeros(len(y))
            for train_idx, test_idx in outer_cv.split(ids, y):
                train_ids, test_ids = ids[train_idx], ids[test_idx]
                y_train = y[train_idx]
                retained_sorted = sorted(select_genera_training_only(raw_tables, {cohort: train_ids}))
                X_train = project_and_clr(raw_tables, cohort, train_ids, retained_sorted)
                X_test = project_and_clr(raw_tables, cohort, test_ids, retained_sorted)
                model, _ = tune_and_fit(X_train, y_train, model_name)
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    oof_proba[test_idx] = model.predict_proba(X_test)[:, 1]

            oof_auc = roc_auc_score(y, oof_proba)
            out_df = pd.DataFrame({"run_id": ids, "y_true": y, "y_score": oof_proba})
            out_path = OUT_DIR_WC / f"{cohort}_{model_name}_oof_predictions.csv"
            out_df.to_csv(out_path, index=False)
            print(f"  {cohort}/{model_name}: recomputed pooled-OOF AUC={oof_auc:.4f}  -> {out_path.name}")
            summary_rows.append({"cohort": cohort, "model": model_name, "auc_recomputed": round(oof_auc, 4)})
    return pd.DataFrame(summary_rows)


def run_loco(raw_tables, diag_df):
    print("\n=== LOCO predictions, training-only feature selection ===")
    summary_rows = []
    labels_by_cohort = {c: labeled_ids_for(diag_df, c) for c in SUPERVISED_COHORTS}

    for held_out in SUPERVISED_COHORTS:
        train_cohorts = [c for c in SUPERVISED_COHORTS if c != held_out]
        train_ids_by_cohort = {c: np.array(labels_by_cohort[c].index) for c in train_cohorts}
        retained_sorted = sorted(select_genera_training_only(raw_tables, train_ids_by_cohort))

        X_train_parts, y_train_parts = [], []
        for c in train_cohorts:
            ids_c = train_ids_by_cohort[c]
            y_c = (labels_by_cohort[c].values == "AD").astype(int)
            X_train_parts.append(project_and_clr(raw_tables, c, ids_c, retained_sorted))
            y_train_parts.append(y_c)
        X_train = np.vstack(X_train_parts)
        y_train = np.concatenate(y_train_parts)

        test_ids = np.array(labels_by_cohort[held_out].index)
        y_test = (labels_by_cohort[held_out].values == "AD").astype(int)
        X_test = project_and_clr(raw_tables, held_out, test_ids, retained_sorted)

        for model_name in ["logreg", "lgbm"]:
            model, _ = tune_and_fit(X_train, y_train, model_name)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                y_score = model.predict_proba(X_test)[:, 1]
            auc = roc_auc_score(y_test, y_score)

            out_df = pd.DataFrame({
                "run_id": test_ids, "y_true": y_test, "y_score": y_score,
                "train_cohorts": "+".join(train_cohorts),
            })
            out_path = OUT_DIR_CC / f"loco_{held_out}_{model_name}_predictions.csv"
            out_df.to_csv(out_path, index=False)
            print(f"  holdout={held_out}/{model_name}: AUC={auc:.4f}  -> {out_path.name}")
            summary_rows.append({"test_cohort": held_out, "model": model_name, "auc_recomputed": round(auc, 4)})

    return pd.DataFrame(summary_rows)


def verify_against_published_summary(wc_check, loco_check):
    print("\n=== Verifying against published 09_training_only_feature_sensitivity.py summary CSVs ===")
    wc_pub = pd.read_csv(TABLES_DIR / "within_cohort_auc_training_only.csv")
    loco_pub = pd.read_csv(TABLES_DIR / "loco_auc_training_only.csv")

    all_ok = True
    for _, r in wc_check.iterrows():
        pub = wc_pub[(wc_pub["cohort"] == r["cohort"]) & (wc_pub["model"] == r["model"])].iloc[0]
        ok = abs(pub["auc_oof_training_only"] - r["auc_recomputed"]) < 1e-4
        all_ok &= ok
        print(f"  within {r['cohort']}/{r['model']}: recomputed={r['auc_recomputed']:.4f} "
              f"published={pub['auc_oof_training_only']:.4f}  {'OK' if ok else 'MISMATCH'}")
    for _, r in loco_check.iterrows():
        pub = loco_pub[(loco_pub["test_cohort"] == r["test_cohort"]) & (loco_pub["model"] == r["model"])].iloc[0]
        ok = abs(pub["auc_training_only"] - r["auc_recomputed"]) < 1e-4
        all_ok &= ok
        print(f"  LOCO {r['test_cohort']}/{r['model']}: recomputed={r['auc_recomputed']:.4f} "
              f"published={pub['auc_training_only']:.4f}  {'OK' if ok else 'MISMATCH'}")

    if all_ok:
        print("\n  ALL VALUES MATCH published training-only summary CSVs exactly (within rounding).")
    else:
        print("\n  WARNING: mismatch(es) detected -- investigate before using these predictions.")
    return all_ok


def main():
    print("16_training_only_oof_predictions.py -- Reviewer 5, Point 5")
    raw_tables = load_raw_tables()
    diag_df = load_labels()

    wc_check = run_within_cohort(raw_tables, diag_df)
    loco_check = run_loco(raw_tables, diag_df)

    ok = verify_against_published_summary(wc_check, loco_check)
    if not ok:
        raise SystemExit("Mismatch against published training-only summary -- aborting.")

    print("\nDone.")


if __name__ == "__main__":
    main()
