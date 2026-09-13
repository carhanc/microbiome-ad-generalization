"""Cross-cohort generalization evaluated on label-blind batch-corrected data
(Reviewer 3, Round 2, concern #2).

Mirrors 07_corrected_generalization.py exactly (same model families, grids,
random state, fold counts, bootstrap CI), but consumes the label-blind
corrected matrices from 10_batch_correction_label_blind.R
(clr_matrix_{combatseq,mmuphin}_labelblind.csv) instead of the original
label-informed ones. Writes to separate, clearly-named output files -- the
original label-informed results (loco_auc_corrected.csv, auc_comparison_table.csv)
are left untouched and are retained as an explicitly-labeled exploratory/
oracle-style sensitivity, per the Round 2 instructions.

This is still a TRANSDUCTIVE design (both correction methods estimate batch
effects jointly across all samples in a batch, including the eventual LOCO
test cohort's own feature values) -- it is "label-blind transductive", not
prospective and not inductive. Do not describe it as either of those.
"""
from pathlib import Path
import warnings

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GridSearchCV, StratifiedKFold, StratifiedGroupKFold
from sklearn.metrics import roc_auc_score
import lightgbm as lgb

PROJECT_ROOT = Path(__file__).resolve().parent.parent
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
TABLES_DIR    = PROJECT_ROOT / "results" / "tables"
TABLES_DIR.mkdir(parents=True, exist_ok=True)

SUPERVISED_COHORTS = ["zhuang2018", "ling2020", "shanghai2022", "kazakhstan2022"]
RANDOM_STATE  = 42
INNER_FOLDS   = 5
N_BOOTSTRAP   = 1000

LOGREG_GRID = {"C": [0.001, 0.01, 0.1, 1.0, 10.0]}
LGBM_GRID   = {
    "n_estimators": [50, 100, 200],
    "learning_rate": [0.05, 0.1],
    "max_depth":     [3, 5],
    "num_leaves":    [15, 31],
}

CORRECTION_METHODS = {
    "combatseq_labelblind": PROCESSED_DIR / "clr_matrix_combatseq_labelblind.csv",
    "mmuphin_labelblind":   PROCESSED_DIR / "clr_matrix_mmuphin_labelblind.csv",
}


def bootstrap_auc_ci(y_true, y_score, n_boot=N_BOOTSTRAP, alpha=0.05):
    rng = np.random.default_rng(RANDOM_STATE)
    n = len(y_true)
    boot = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        yt, ys = y_true[idx], y_score[idx]
        if len(np.unique(yt)) < 2:
            continue
        boot.append(roc_auc_score(yt, ys))
    boot = np.array(boot)
    return float(np.percentile(boot, 100*alpha/2)), float(np.percentile(boot, 100*(1-alpha/2)))


def tune_and_fit(X_train, y_train, model_name):
    inner_cv = StratifiedKFold(n_splits=INNER_FOLDS, shuffle=True, random_state=RANDOM_STATE)
    if model_name == "logreg":
        # NOTE: penalty="elasticnet" is never set -- this is L2, l1_ratio is inert.
        # Matches the historical model construction exactly for comparability.
        base  = LogisticRegression(solver="saga", l1_ratio=0.5,
                                   max_iter=10000, random_state=RANDOM_STATE)
        grid  = LOGREG_GRID
    else:
        base  = lgb.LGBMClassifier(objective="binary", n_jobs=4,
                                    random_state=RANDOM_STATE, verbose=-1)
        grid  = LGBM_GRID
    gs = GridSearchCV(base, grid, cv=inner_cv, scoring="roc_auc",
                      n_jobs=-1, refit=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        gs.fit(X_train, y_train)
    return gs.best_estimator_, gs.best_params_


def eval_test(model, X_test, y_test):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        y_score = model.predict_proba(X_test)[:, 1]
    auc = roc_auc_score(y_test, y_score)
    ci_lo, ci_hi = bootstrap_auc_ci(y_test, y_score)
    return {"auc": round(auc, 4), "ci_lo": round(ci_lo, 4),
            "ci_hi": round(ci_hi, 4), "y_score": y_score}


def load_corrected_data(clr_path):
    df = pd.read_csv(clr_path)
    id_col = df.columns[0]
    genus_cols = [c for c in df.columns if c not in [id_col, "cohort", "diagnosis"]]

    by_cohort = {}
    for cohort in SUPERVISED_COHORTS:
        sub = df[(df["cohort"] == cohort) & df["diagnosis"].isin(["AD", "CN"])].copy()
        labels = (sub["diagnosis"] == "AD").astype(int).to_numpy()
        by_cohort[cohort] = {
            "X": sub[genus_cols].to_numpy(),
            "y": labels,
            "ids":  list(sub[id_col]),
            "n":    len(labels),
            "n_AD": int(labels.sum()),
            "n_CN": int((labels == 0).sum()),
        }
    return by_cohort, genus_cols


def run_loco(cohort_data, method_label):
    results = []
    for held_out in SUPERVISED_COHORTS:
        train_cohorts = [c for c in SUPERVISED_COHORTS if c != held_out]
        X_train = np.vstack([cohort_data[c]["X"] for c in train_cohorts])
        y_train = np.concatenate([cohort_data[c]["y"] for c in train_cohorts])
        X_test  = cohort_data[held_out]["X"]
        y_test  = cohort_data[held_out]["y"]

        print(f"  {method_label} LOCO held-out={held_out} "
              f"| train n={len(y_train)}, test n={len(y_test)}")

        for model_name in ["logreg", "lgbm"]:
            model, best_p = tune_and_fit(X_train, y_train, model_name)
            m = eval_test(model, X_test, y_test)
            print(f"    {model_name} AUC={m['auc']:.4f} "
                  f"CI=[{m['ci_lo']:.4f},{m['ci_hi']:.4f}]")
            results.append({
                "correction":    method_label,
                "experiment":    "LOCO",
                "train_cohorts": "+".join(train_cohorts),
                "test_cohort":   held_out,
                "model":         model_name,
                "n_train":       len(y_train),
                "n_test":        len(y_test),
                "auc":           m["auc"],
                "ci_lo":         m["ci_lo"],
                "ci_hi":         m["ci_hi"],
            })
    return results


def run_pairwise(cohort_data, method_label):
    results = []
    for train_c in SUPERVISED_COHORTS:
        for test_c in SUPERVISED_COHORTS:
            if train_c == test_c:
                continue
            X_train = cohort_data[train_c]["X"]
            y_train = cohort_data[train_c]["y"]
            X_test  = cohort_data[test_c]["X"]
            y_test  = cohort_data[test_c]["y"]

            for model_name in ["logreg", "lgbm"]:
                model, _ = tune_and_fit(X_train, y_train, model_name)
                m = eval_test(model, X_test, y_test)
                results.append({
                    "correction":   method_label,
                    "experiment":   "pairwise",
                    "train_cohort": train_c,
                    "test_cohort":  test_c,
                    "model":        model_name,
                    "n_train":      len(y_train),
                    "n_test":       len(y_test),
                    "auc":          m["auc"],
                    "ci_lo":        m["ci_lo"],
                    "ci_hi":        m["ci_hi"],
                })
    return results


def run_within_cohort_check(cohort_data, method_label, diag_df):
    OUTER_FOLDS = 10
    results = []
    for cohort in SUPERVISED_COHORTS:
        cd   = cohort_data[cohort]
        X, y = cd["X"], cd["y"]

        groups = None
        if cohort == "shanghai2022":
            run_id_to_sname = diag_df.set_index("run_id")["sample_name"]
            snames = pd.Series(cd["ids"]).map(run_id_to_sname)
            groups = snames.map(lambda s: "_".join(s.split("_")[1:])).to_numpy()

        outer_cv = (StratifiedGroupKFold(n_splits=OUTER_FOLDS) if groups is not None
                    else StratifiedKFold(n_splits=OUTER_FOLDS, shuffle=True,
                                         random_state=RANDOM_STATE))
        inner_cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=RANDOM_STATE)

        for model_name in ["logreg", "lgbm"]:
            oof = np.zeros(len(y))
            split_iter = (outer_cv.split(X, y, groups=groups) if groups is not None
                          else outer_cv.split(X, y))
            for train_i, test_i in split_iter:
                Xtr, Xte = X[train_i], X[test_i]
                ytr, yte = y[train_i], y[test_i]
                model, _ = tune_and_fit(Xtr, ytr, model_name)
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    oof[test_i] = model.predict_proba(Xte)[:, 1]
            auc = roc_auc_score(y, oof)
            ci_lo, ci_hi = bootstrap_auc_ci(y, oof)
            print(f"  {method_label} within {cohort}/{model_name}: "
                  f"AUC={auc:.4f} CI=[{ci_lo:.4f},{ci_hi:.4f}]")
            results.append({
                "correction": method_label,
                "cohort":     cohort,
                "model":      model_name,
                "auc_within_corrected": round(auc, 4),
                "ci_lo":      round(ci_lo, 4),
                "ci_hi":      round(ci_hi, 4),
            })
    return results


def main():
    print("11_corrected_generalization_label_blind.py")

    diag_df = pd.read_csv(PROCESSED_DIR / "diagnosis_labels.csv")

    uncorrected_loco = pd.read_csv(TABLES_DIR / "loco_auc.csv")
    within_ref       = pd.read_csv(TABLES_DIR / "within_cohort_auc.csv")
    within_ref       = within_ref[~within_ref["cohort"].str.endswith("3class")]
    within_ref_lookup = within_ref.set_index(["cohort", "model"])["auc_oof"]

    all_loco, all_pairwise, all_within = [], [], []

    for method, clr_path in CORRECTION_METHODS.items():
        if not clr_path.exists():
            print(f"\nSkipping {method}: {clr_path} not found. "
                  f"Run 10_batch_correction_label_blind.R first.")
            continue

        print(f"\nMethod: {method}")
        cohort_data, genus_cols = load_corrected_data(clr_path)
        print(f"  Loaded {sum(cd['n'] for cd in cohort_data.values())} samples, "
              f"{len(genus_cols)} genera")

        print("\n  LOCO")
        all_loco += run_loco(cohort_data, method)
        print("\n  Pairwise")
        all_pairwise += run_pairwise(cohort_data, method)
        print("\n  Within-cohort check on label-blind corrected data")
        all_within += run_within_cohort_check(cohort_data, method, diag_df)

    loco_df     = pd.DataFrame(all_loco)
    pairwise_df = pd.DataFrame(all_pairwise)
    within_df   = pd.DataFrame(all_within)

    loco_df.to_csv(TABLES_DIR / "loco_auc_corrected_labelblind.csv", index=False)
    pairwise_df.to_csv(TABLES_DIR / "pairwise_auc_labelblind.csv", index=False)
    within_df.to_csv(TABLES_DIR / "within_cohort_auc_corrected_labelblind.csv", index=False)
    print(f"\n  wrote {TABLES_DIR}/loco_auc_corrected_labelblind.csv")
    print(f"  wrote {TABLES_DIR}/pairwise_auc_labelblind.csv")
    print(f"  wrote {TABLES_DIR}/within_cohort_auc_corrected_labelblind.csv")

    # Comparison table: uncorrected vs label-blind-corrected (existing
    # label-informed results are in the separate auc_comparison_table.csv
    # and are NOT merged into this file)
    uncorr_tagged = uncorrected_loco.copy()
    uncorr_tagged["correction"] = "uncorrected"

    rows = []
    for _, unc_row in uncorr_tagged.iterrows():
        cohort, model = unc_row["test_cohort"], unc_row["model"]
        wc_auc = within_ref_lookup.get((cohort, model), np.nan)
        rows.append({
            "test_cohort": cohort, "model": model, "within_auc": round(wc_auc, 4),
            "correction": "uncorrected", "loco_auc": unc_row["auc"],
            "loco_ci_lo": unc_row["ci_lo"], "loco_ci_hi": unc_row["ci_hi"],
            "loco_drop": round(wc_auc - unc_row["auc"], 4),
        })
    for _, corr_row in loco_df.iterrows():
        cohort, model = corr_row["test_cohort"], corr_row["model"]
        wc_auc = within_ref_lookup.get((cohort, model), np.nan)
        rows.append({
            "test_cohort": cohort, "model": model, "within_auc": round(wc_auc, 4),
            "correction": corr_row["correction"], "loco_auc": corr_row["auc"],
            "loco_ci_lo": corr_row["ci_lo"], "loco_ci_hi": corr_row["ci_hi"],
            "loco_drop": round(wc_auc - corr_row["auc"], 4),
        })

    comparison_df = pd.DataFrame(rows)
    comparison_df.to_csv(TABLES_DIR / "auc_comparison_labelblind.csv", index=False)
    print(f"  wrote {TABLES_DIR}/auc_comparison_labelblind.csv")

    print("\nKey findings (label-blind transductive correction)")
    for model in ["logreg", "lgbm"]:
        sub = comparison_df[comparison_df["model"] == model]
        print(f"\n  Model: {model}")
        for correction in ["uncorrected"] + list(CORRECTION_METHODS.keys()):
            csub = sub[sub["correction"] == correction]
            if csub.empty:
                continue
            print(f"    {correction:24s}  mean LOCO AUC={csub['loco_auc'].mean():.4f}  "
                  f"mean drop={csub['loco_drop'].mean():+.4f}")

    print("\nDone.")


if __name__ == "__main__":
    main()
