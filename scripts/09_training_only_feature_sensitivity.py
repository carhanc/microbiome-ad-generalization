"""Strict training-only genus-selection sensitivity analysis (Reviewer 3, Round 2, concern #1).

The original pipeline (01_harmonize_taxonomy.py) selects the retained 396-genus
feature universe ONCE, globally, using raw genus tables from all five cohorts
simultaneously (>=20% prevalence in >=1 cohort). That feature universe is then
reused, unchanged, by every downstream supervised analysis (03, 04, 07, 08) —
so for every LOCO/pairwise/nested-CV split, the held-out cohort's own raw
prevalence contributed to deciding which genera exist as features at all.

This script repeats the three supervised analyses (within-cohort nested CV,
LOCO, pairwise transfer) with genus selection re-derived independently inside
every single split, using ONLY the samples on the training side of that split.
The held-out/test side contributes nothing to feature selection. Test/held-out
samples are then projected onto the training-derived genus set (absent genera
filled with zero count), and zero-replacement (delta=0.65) + CLR are applied
per-sample, independently, AFTER projection — using the exact functions from
02_clr_transform.py, unmodified.

Genus universe used for prevalence computation in every split = the AD/CN-
labeled samples that would actually enter that split's classifier (i.e., the
exact X/y row set), not the full raw per-cohort table (which includes
MCI/other/unlabeled samples). This is a deliberate, documented design choice:
it means "training-only" is defined relative to the samples the classifier
actually trains on, which is the safest and least ambiguous reading of
"the held-out cohort must contribute nothing to genus selection." kbase2022
(unlabeled, not part of any supervised analysis) is excluded from this
script entirely, exactly as it already is from 03/04/07/08.

Same model families, random state, hyperparameter grids, fold counts, and
bootstrap CI procedure as the historical scripts (03_within_cohort_baseline.py,
04_cross_cohort_generalization.py) — the only change under test is the
feature-selection scoping. DADA2 is NOT rerun; this script reads the
git-tracked pre-filter per-cohort tables (data/processed/dada2/<cohort>/genus_table.csv)
and reuses 01_harmonize_taxonomy.py's synonym-map / uninformative-genus-drop /
prevalence-filter functions and 02_clr_transform.py's zero-replacement/CLR
functions unmodified, via direct import.
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
TABLES_DIR.mkdir(parents=True, exist_ok=True)


def _load_module(name, filename):
    spec = importlib.util.spec_from_file_location(name, PROJECT_ROOT / "scripts" / filename)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


harmonize = _load_module("harmonize_09", "01_harmonize_taxonomy.py")
clr_mod = _load_module("clr_transform_09", "02_clr_transform.py")

SUPERVISED_COHORTS = ["zhuang2018", "ling2020", "shanghai2022", "kazakhstan2022"]
RANDOM_STATE = 42
OUTER_FOLDS = 10
INNER_FOLDS = 5
N_BOOTSTRAP = 1000
MIN_PREVALENCE = 0.20

LOGREG_GRID = {"C": [0.001, 0.01, 0.1, 1.0, 10.0]}
LGBM_GRID = {
    "n_estimators": [50, 100, 200],
    "learning_rate": [0.05, 0.1],
    "max_depth": [3, 5],
    "num_leaves": [15, 31],
}

feature_count_rows = []
feature_list_rows = []


def bootstrap_auc_ci(y_true, y_score, n_boot=N_BOOTSTRAP, alpha=0.05, rng=None):
    if rng is None:
        rng = np.random.default_rng(RANDOM_STATE)
    n = len(y_true)
    auc_boot = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        yt, ys = y_true[idx], y_score[idx]
        if len(np.unique(yt)) < 2:
            continue
        auc_boot.append(roc_auc_score(yt, ys))
    auc_boot = np.array(auc_boot)
    lo = float(np.percentile(auc_boot, 100 * alpha / 2))
    hi = float(np.percentile(auc_boot, 100 * (1 - alpha / 2)))
    return lo, hi


def load_raw_tables():
    """Per-cohort raw genus tables, post synonym-map + uninformative-genus drop,
    PRE prevalence-filter. Indexed by run_id. Identical preprocessing to
    01_harmonize_taxonomy.py::harmonize_all up to (not including) the
    prevalence filter step."""
    tables = {}
    for cohort in SUPERVISED_COHORTS:
        df = harmonize.load_genus_table(cohort)
        df = harmonize.apply_synonym_map(df)
        df = harmonize.drop_uninformative_genera(df)
        tables[cohort] = df
    return tables


def load_labels():
    diag = pd.read_csv(PROCESSED_DIR / "diagnosis_labels.csv")
    return diag


def labeled_ids_for(diag_df, cohort):
    sub = diag_df[(diag_df["cohort"] == cohort) & (diag_df["diagnosis"].isin(["AD", "CN"]))]
    return sub.set_index("run_id")["diagnosis"]


def select_genera_training_only(raw_tables, sample_ids_by_cohort):
    """Reuses harmonize.apply_prevalence_filter's exact per-cohort-then-union
    logic (>=20% prevalence in >=1 cohort), restricted to the given sample IDs
    per cohort (the training-side rows of the current split only)."""
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
    clr_arr = clr_mod.clr_transform(counts_replaced)
    return clr_arr


def tune_and_fit(X_train, y_train, model_name):
    inner_cv = StratifiedKFold(n_splits=INNER_FOLDS, shuffle=True, random_state=RANDOM_STATE)
    if model_name == "logreg":
        # NOTE: penalty="elasticnet" is never set -- this is L2, l1_ratio is inert.
        # Matches the historical model construction exactly for comparability.
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


def record_feature_set(experiment, split_id, retained_genera_sorted):
    feature_count_rows.append({
        "experiment": experiment,
        "split_id": split_id,
        "n_retained_genera": len(retained_genera_sorted),
    })
    for g in retained_genera_sorted:
        feature_list_rows.append({"experiment": experiment, "split_id": split_id, "genus": g})


# ---------------------------------------------------------------------------
# A. Within-cohort nested-CV sensitivity
# ---------------------------------------------------------------------------

def run_within_cohort_sensitivity(raw_tables, diag_df):
    print("\n=== A. Within-cohort nested-CV, training-only feature selection ===")
    results = []

    for cohort in SUPERVISED_COHORTS:
        labels = labeled_ids_for(diag_df, cohort)
        ids = np.array(labels.index)
        y = (labels.values == "AD").astype(int)
        n_ad, n_cn = int(y.sum()), int((y == 0).sum())
        print(f"\n  {cohort}  n={len(y)} ({n_ad} AD / {n_cn} CN)")
        if n_ad < 5 or n_cn < 5:
            print("    too few per class, skipping")
            continue

        outer_cv = StratifiedKFold(n_splits=OUTER_FOLDS, shuffle=True, random_state=RANDOM_STATE)

        for model_name in ["logreg", "lgbm"]:
            oof_proba = np.zeros(len(y))
            fold_aucs = []
            for fold_idx, (train_idx, test_idx) in enumerate(outer_cv.split(ids, y)):
                train_ids = ids[train_idx]
                test_ids = ids[test_idx]
                y_train, y_test = y[train_idx], y[test_idx]

                # feature selection uses ONLY the outer-training samples of
                # THIS cohort (the outer-test fold is completely untouched)
                retained = select_genera_training_only(raw_tables, {cohort: train_ids})
                retained_sorted = sorted(retained)
                if model_name == "logreg":
                    record_feature_set(f"within_cohort_{cohort}", f"outer_fold_{fold_idx}", retained_sorted)

                X_train = project_and_clr(raw_tables, cohort, train_ids, retained_sorted)
                X_test = project_and_clr(raw_tables, cohort, test_ids, retained_sorted)

                model, best_params = tune_and_fit(X_train, y_train, model_name)
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    fold_proba = model.predict_proba(X_test)[:, 1]
                oof_proba[test_idx] = fold_proba
                fold_auc = roc_auc_score(y_test, fold_proba)
                fold_aucs.append(fold_auc)
                print(f"    [{model_name}] fold {fold_idx+1:2d}/{OUTER_FOLDS}: "
                      f"AUC={fold_auc:.4f}  n_genera={len(retained_sorted)}  best={best_params}")

            oof_auc = roc_auc_score(y, oof_proba)
            rng = np.random.default_rng(RANDOM_STATE)
            ci_lo, ci_hi = bootstrap_auc_ci(y, oof_proba, rng=rng)
            print(f"    [{model_name}] pooled OOF AUC = {oof_auc:.4f}  95% CI [{ci_lo:.4f}, {ci_hi:.4f}]")

            results.append({
                "cohort": cohort, "model": model_name, "n_samples": len(y),
                "n_AD": n_ad, "n_CN": n_cn,
                "auc_oof_training_only": round(oof_auc, 4),
                "auc_ci_lo_training_only": round(ci_lo, 4),
                "auc_ci_hi_training_only": round(ci_hi, 4),
                "auc_mean_folds_training_only": round(float(np.mean(fold_aucs)), 4),
            })

    return pd.DataFrame(results)


# ---------------------------------------------------------------------------
# B. LOCO sensitivity
# ---------------------------------------------------------------------------

def run_loco_sensitivity(raw_tables, diag_df):
    print("\n=== B. LOCO, training-only feature selection ===")
    results = []
    labels_by_cohort = {c: labeled_ids_for(diag_df, c) for c in SUPERVISED_COHORTS}

    for held_out in SUPERVISED_COHORTS:
        train_cohorts = [c for c in SUPERVISED_COHORTS if c != held_out]
        print(f"\n  LOCO — hold out {held_out}, train on {train_cohorts}")

        train_ids_by_cohort = {c: np.array(labels_by_cohort[c].index) for c in train_cohorts}
        retained = select_genera_training_only(raw_tables, train_ids_by_cohort)
        retained_sorted = sorted(retained)
        record_feature_set("loco", f"holdout_{held_out}", retained_sorted)
        print(f"    n_retained_genera (from {train_cohorts} only) = {len(retained_sorted)}")

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
            model, best_params = tune_and_fit(X_train, y_train, model_name)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                y_score = model.predict_proba(X_test)[:, 1]
            auc = roc_auc_score(y_test, y_score)
            ci_lo, ci_hi = bootstrap_auc_ci(y_test, y_score)
            print(f"    [{model_name}] AUC={auc:.4f}  CI=[{ci_lo:.4f}, {ci_hi:.4f}]  best={best_params}")
            results.append({
                "test_cohort": held_out, "train_cohorts": "+".join(train_cohorts),
                "model": model_name, "n_train": len(y_train), "n_test": len(y_test),
                "n_retained_genera": len(retained_sorted),
                "auc_training_only": round(auc, 4),
                "ci_lo_training_only": round(ci_lo, 4),
                "ci_hi_training_only": round(ci_hi, 4),
            })

    return pd.DataFrame(results)


# ---------------------------------------------------------------------------
# C. Pairwise sensitivity
# ---------------------------------------------------------------------------

def run_pairwise_sensitivity(raw_tables, diag_df):
    print("\n=== C. Pairwise transfer, training-only feature selection ===")
    results = []
    labels_by_cohort = {c: labeled_ids_for(diag_df, c) for c in SUPERVISED_COHORTS}

    for train_cohort in SUPERVISED_COHORTS:
        train_ids = np.array(labels_by_cohort[train_cohort].index)
        retained = select_genera_training_only(raw_tables, {train_cohort: train_ids})
        retained_sorted = sorted(retained)
        record_feature_set("pairwise", f"train_{train_cohort}", retained_sorted)

        y_train = (labels_by_cohort[train_cohort].values == "AD").astype(int)
        X_train = project_and_clr(raw_tables, train_cohort, train_ids, retained_sorted)

        for test_cohort in SUPERVISED_COHORTS:
            if test_cohort == train_cohort:
                continue
            test_ids = np.array(labels_by_cohort[test_cohort].index)
            y_test = (labels_by_cohort[test_cohort].values == "AD").astype(int)
            X_test = project_and_clr(raw_tables, test_cohort, test_ids, retained_sorted)

            print(f"\n  {train_cohort} -> {test_cohort}  (n_genera from {train_cohort} only = {len(retained_sorted)})")
            for model_name in ["logreg", "lgbm"]:
                model, best_params = tune_and_fit(X_train, y_train, model_name)
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    y_score = model.predict_proba(X_test)[:, 1]
                auc = roc_auc_score(y_test, y_score)
                ci_lo, ci_hi = bootstrap_auc_ci(y_test, y_score)
                print(f"    [{model_name}] AUC={auc:.4f}  CI=[{ci_lo:.4f}, {ci_hi:.4f}]")
                results.append({
                    "train_cohort": train_cohort, "test_cohort": test_cohort,
                    "model": model_name, "n_train": len(y_train), "n_test": len(y_test),
                    "n_retained_genera": len(retained_sorted),
                    "auc_training_only": round(auc, 4),
                    "ci_lo_training_only": round(ci_lo, 4),
                    "ci_hi_training_only": round(ci_hi, 4),
                })

    return pd.DataFrame(results)


def build_summary(within_df, loco_df, pairwise_df):
    print("\n=== Summary: original (global-feature) vs strict training-only ===")
    orig_within = pd.read_csv(TABLES_DIR / "within_cohort_auc.csv")
    orig_within = orig_within[~orig_within["cohort"].str.endswith("3class")]
    orig_loco = pd.read_csv(TABLES_DIR / "loco_auc.csv")

    rows = []
    for _, r in within_df.iterrows():
        orig = orig_within[(orig_within["cohort"] == r["cohort"]) & (orig_within["model"] == r["model"])]
        if len(orig) == 0:
            continue
        o = orig.iloc[0]
        rows.append({
            "experiment": "within_cohort",
            "cohort": r["cohort"], "model": r["model"],
            "original_auc": o["auc_oof"], "training_only_auc": r["auc_oof_training_only"],
            "delta": round(r["auc_oof_training_only"] - o["auc_oof"], 4),
            "original_ci": f"[{o['auc_ci_lo']:.4f}, {o['auc_ci_hi']:.4f}]",
            "training_only_ci": f"[{r['auc_ci_lo_training_only']:.4f}, {r['auc_ci_hi_training_only']:.4f}]",
            "n_retained_genera_training_only": None,
        })

    for _, r in loco_df.iterrows():
        orig = orig_loco[(orig_loco["test_cohort"] == r["test_cohort"]) & (orig_loco["model"] == r["model"])]
        if len(orig) == 0:
            continue
        o = orig.iloc[0]
        rows.append({
            "experiment": "loco",
            "cohort": r["test_cohort"], "model": r["model"],
            "original_auc": o["auc"], "training_only_auc": r["auc_training_only"],
            "delta": round(r["auc_training_only"] - o["auc"], 4),
            "original_ci": f"[{o['ci_lo']:.4f}, {o['ci_hi']:.4f}]",
            "training_only_ci": f"[{r['ci_lo_training_only']:.4f}, {r['ci_hi_training_only']:.4f}]",
            "n_retained_genera_training_only": r["n_retained_genera"],
        })

    summary_df = pd.DataFrame(rows)
    return summary_df


def main():
    print("09_training_only_feature_sensitivity.py — Reviewer 3 Round 2, concern #1")

    raw_tables = load_raw_tables()
    diag_df = load_labels()

    within_df = run_within_cohort_sensitivity(raw_tables, diag_df)
    within_df.to_csv(TABLES_DIR / "within_cohort_auc_training_only.csv", index=False)
    print(f"\n  wrote {TABLES_DIR}/within_cohort_auc_training_only.csv")

    loco_df = run_loco_sensitivity(raw_tables, diag_df)
    loco_df.to_csv(TABLES_DIR / "loco_auc_training_only.csv", index=False)
    print(f"  wrote {TABLES_DIR}/loco_auc_training_only.csv")

    pairwise_df = run_pairwise_sensitivity(raw_tables, diag_df)
    pairwise_df.to_csv(TABLES_DIR / "pairwise_auc_training_only.csv", index=False)
    print(f"  wrote {TABLES_DIR}/pairwise_auc_training_only.csv")

    feat_counts_df = pd.DataFrame(feature_count_rows)
    feat_counts_df.to_csv(TABLES_DIR / "training_only_feature_counts.csv", index=False)
    print(f"  wrote {TABLES_DIR}/training_only_feature_counts.csv "
          f"({len(feat_counts_df)} splits)")

    feat_list_df = pd.DataFrame(feature_list_rows)
    feat_list_df.to_csv(TABLES_DIR / "training_only_feature_lists.csv", index=False)
    print(f"  wrote {TABLES_DIR}/training_only_feature_lists.csv "
          f"({len(feat_list_df)} rows, full retained-genus lists per split, for leakage audit)")

    summary_df = build_summary(within_df, loco_df, pairwise_df)
    summary_df.to_csv(TABLES_DIR / "training_only_feature_sensitivity_summary.csv", index=False)
    print(f"  wrote {TABLES_DIR}/training_only_feature_sensitivity_summary.csv")
    print("\n" + summary_df.to_string(index=False))

    print("\nDone.")


if __name__ == "__main__":
    main()
