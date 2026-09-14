"""Zhu 2022 robustness analysis (Reviewer 4, major statistical request).

The historical single nested-CV partition gave an unusually high within-cohort
AUC for Zhu 2022 fecal AD/CN participants (n=60: 30 AD, 30 CN) -- LogReg
~0.998, LightGBM ~0.979. This script assesses whether that figure is (a)
robust to outer-fold assignment (repeated nested CV) and (b) distinguishable
from chance after label permutation, under the SAME strict training-only
feature-selection pipeline used throughout Round 2 (genus retention
re-derived from only the outer-training participants in every fold; the
outer-test participants never contribute to feature selection or tuning).

Reuses 01_harmonize_taxonomy.py's synonym-map/drop/prevalence-filter and
02_clr_transform.py's zero-replacement/CLR functions, and
09_training_only_feature_sensitivity.py's select_genera_training_only /
project_and_clr / tune_and_fit helpers, unmodified, via direct import --
this is the exact same strict pipeline already used for the cross-cohort
training-only sensitivity, applied here to a single cohort's own internal
nested CV.

Key implementation property, stated explicitly because it is what makes the
permutation test computationally tractable: strict training-only genus
selection and the resulting CLR feature matrices depend ONLY on which
participants are in the training fold, not on their diagnosis labels. So for
a fixed outer-CV partition, the per-fold (X_train, X_test) matrices are
identical across every label permutation and are computed once and reused --
only model tuning/fitting (which does depend on y) is repeated per
permutation. This is a pure implementation optimization; it does not change
what is estimated or relax the leakage constraint in any way.
"""
import sys
import time
import warnings
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GridSearchCV, StratifiedKFold
from sklearn.metrics import roc_auc_score
import lightgbm as lgb

PROJECT_ROOT = Path(__file__).resolve().parent.parent
TABLES_DIR = PROJECT_ROOT / "results" / "tables"
FIGS_DIR = PROJECT_ROOT / "results" / "figures"
TABLES_DIR.mkdir(parents=True, exist_ok=True)

sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

import importlib.util


def _load_module(name, filename):
    spec = importlib.util.spec_from_file_location(name, PROJECT_ROOT / "scripts" / filename)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


harmonize = _load_module("harmonize_12", "01_harmonize_taxonomy.py")
clr_mod = _load_module("clr_transform_12", "02_clr_transform.py")
train_only = _load_module("train_only_12", "09_training_only_feature_sensitivity.py")

COHORT_KEY = "shanghai2022"  # Zhu 2022's internal cohort key throughout this repo
OUTER_FOLDS = 10
INNER_FOLDS = 5
INNER_SEED = 42  # matches RANDOM_STATE used throughout the manuscript pipeline

N_REPEATS = 50
# Pre-specified deterministic seed sequence for the repeated-CV outer partitions.
REPEAT_SEEDS = list(range(1000, 1000 + N_REPEATS))

N_PERM_TARGET = 1000
PERM_SEED = 20260914  # single seed for the permutation RNG stream (label shuffles only)


def load_zhu_data():
    df = harmonize.load_genus_table(COHORT_KEY)
    df = harmonize.apply_synonym_map(df)
    df = harmonize.drop_uninformative_genera(df)

    diag = pd.read_csv(PROJECT_ROOT / "data" / "processed" / "diagnosis_labels.csv")
    labels = diag[(diag["cohort"] == COHORT_KEY) & (diag["diagnosis"].isin(["AD", "CN"]))]
    labels = labels.set_index("run_id")["diagnosis"]
    ids = np.array(labels.index)
    y = (labels.values == "AD").astype(int)
    assert len(ids) == 60 and y.sum() == 30 and (y == 0).sum() == 30, \
        f"expected 60 participants (30 AD/30 CN), got {len(ids)} ({y.sum()} AD)"
    return df, ids, y


def precompute_fold_features(raw_table, ids, outer_seed):
    """For a given outer-CV seed, precompute (train_idx, test_idx, X_train,
    X_test, n_genera) per fold. Depends only on participant membership, not
    on labels -- reused across every label permutation for this outer_seed."""
    # StratifiedKFold needs a y to stratify on for split(); membership itself
    # doesn't depend on which y is used as long as the *real* y is used to
    # define the split (matching the historical/main analysis's own fold
    # structure), so this always uses the TRUE labels to build the partition.
    _, true_ids, true_y = load_zhu_data()
    assert np.array_equal(true_ids, ids)
    outer_cv = StratifiedKFold(n_splits=OUTER_FOLDS, shuffle=True, random_state=outer_seed)
    folds = []
    for train_idx, test_idx in outer_cv.split(ids, true_y):
        train_ids = ids[train_idx]
        test_ids = ids[test_idx]
        retained = train_only.select_genera_training_only(
            {COHORT_KEY: raw_table}, {COHORT_KEY: train_ids}
        )
        retained_sorted = sorted(retained)
        X_train = train_only.project_and_clr({COHORT_KEY: raw_table}, COHORT_KEY, train_ids, retained_sorted)
        X_test = train_only.project_and_clr({COHORT_KEY: raw_table}, COHORT_KEY, test_ids, retained_sorted)
        folds.append({
            "train_idx": train_idx, "test_idx": test_idx,
            "X_train": X_train, "X_test": X_test,
            "n_genera": len(retained_sorted),
        })
    return folds


LOGREG_GRID = {"C": [0.001, 0.01, 0.1, 1.0, 10.0]}
LGBM_GRID = {
    "n_estimators": [50, 100, 200], "learning_rate": [0.05, 0.1],
    "max_depth": [3, 5], "num_leaves": [15, 31],
}


def tune_and_fit_serial(X_train, y_train, model_name):
    """Identical model definitions/grids/inner-CV structure to
    09_training_only_feature_sensitivity.py::tune_and_fit (which itself
    mirrors the historical 03/04 scripts) -- the only difference is
    n_jobs=1 throughout instead of n_jobs=-1/4. This is a pure performance
    choice: for this dataset's tiny per-fit cost, serial execution measured
    as fast as or faster than internal parallelism (see implementation log),
    and forcing n_jobs=1 here is what makes it safe to parallelize across
    repetitions/permutations at the process level without CPU oversubscription.
    GridSearchCV(n_jobs=1) vs (n_jobs=-1) fit identical models and produce
    identical results -- n_jobs never affects which hyperparameters are
    selected or any statistical outcome, only wall-clock time."""
    inner_cv = StratifiedKFold(n_splits=INNER_FOLDS, shuffle=True, random_state=INNER_SEED)
    if model_name == "logreg":
        base = LogisticRegression(solver="saga", l1_ratio=0.5, max_iter=10000, random_state=INNER_SEED)
        grid = LOGREG_GRID
    else:
        base = lgb.LGBMClassifier(objective="binary", n_jobs=1, random_state=INNER_SEED, verbose=-1)
        grid = LGBM_GRID
    gs = GridSearchCV(base, grid, cv=inner_cv, scoring="roc_auc", n_jobs=1, refit=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        gs.fit(X_train, y_train)
    return gs.best_estimator_, gs.best_params_


def run_nested_cv_given_folds(folds, y, model_name):
    n = len(y)
    oof_proba = np.zeros(n)
    for f in folds:
        y_train = y[f["train_idx"]]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            model, _ = tune_and_fit_serial(f["X_train"], y_train, model_name)
            proba = model.predict_proba(f["X_test"])[:, 1]
        oof_proba[f["test_idx"]] = proba
    return roc_auc_score(y, oof_proba)


# ---------------------------------------------------------------------------
# A. Repeated nested CV (fold-assignment stability)
# ---------------------------------------------------------------------------

def _run_one_repeat(args):
    raw_table, ids, y, i, seed = args
    folds = precompute_fold_features(raw_table, ids, seed)
    n_genera_list = [f["n_genera"] for f in folds]
    out = []
    for model_name in ["logreg", "lgbm"]:
        auc = run_nested_cv_given_folds(folds, y, model_name)
        out.append({
            "repeat": i, "outer_seed": seed, "model": model_name,
            "auc": round(auc, 4),
            "mean_n_genera": round(float(np.mean(n_genera_list)), 1),
        })
    return out


def run_repeated_cv(raw_table, ids, y):
    print(f"\n=== A. Repeated nested CV, {N_REPEATS} repetitions ===")
    rows = []
    t0 = time.time()
    args_list = [(raw_table, ids, y, i, seed) for i, seed in enumerate(REPEAT_SEEDS)]
    with ProcessPoolExecutor(max_workers=10) as ex:
        futures = {ex.submit(_run_one_repeat, a): a[3] for a in args_list}
        done = 0
        for fut in as_completed(futures):
            rows.extend(fut.result())
            done += 1
            if done % 10 == 0:
                elapsed = time.time() - t0
                print(f"  repeat {done}/{N_REPEATS}  ({elapsed:.1f}s elapsed)")
    rows.sort(key=lambda r: (r["repeat"], r["model"]))
    return pd.DataFrame(rows)


def summarize_repeated_cv(df, historical_auc, training_only_single_auc):
    rows = []
    for model_name in ["logreg", "lgbm"]:
        sub = df[df["model"] == model_name]["auc"]
        rows.append({
            "model": model_name,
            "historical_single_partition_auc": historical_auc[model_name],
            "training_only_single_partition_auc": training_only_single_auc[model_name],
            "repeated_cv_mean": round(sub.mean(), 4),
            "repeated_cv_sd": round(sub.std(), 4),
            "repeated_cv_median": round(sub.median(), 4),
            "repeated_cv_p2.5": round(sub.quantile(0.025), 4),
            "repeated_cv_p97.5": round(sub.quantile(0.975), 4),
            "repeated_cv_min": round(sub.min(), 4),
            "repeated_cv_max": round(sub.max(), 4),
            "n_repeats": len(sub),
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# B. Label permutation test
# ---------------------------------------------------------------------------

def _run_one_permutation(args):
    folds_light, y_perm, model_name = args
    # folds_light carries only arrays (picklable), rebuilt minimally here
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        n = len(y_perm)
        oof_proba = np.zeros(n)
        for f in folds_light:
            y_train = y_perm[f["train_idx"]]
            model, _ = tune_and_fit_serial(f["X_train"], y_train, model_name)
            oof_proba[f["test_idx"]] = model.predict_proba(f["X_test"])[:, 1]
        return roc_auc_score(y_perm, oof_proba)


def run_permutation_test(raw_table, ids, y, n_perm, observed_auc):
    print(f"\n=== B. Label permutation test, target n_perm={n_perm} ===")
    # Fixed outer-CV partition matching the main/historical analysis's own
    # fold-assignment seed (42) -- only labels are permuted across replicates.
    folds = precompute_fold_features(raw_table, ids, outer_seed=42)
    print(f"  Precomputed {len(folds)} folds once (label-independent); "
          f"mean retained genera = {np.mean([f['n_genera'] for f in folds]):.1f}")

    rng = np.random.default_rng(PERM_SEED)
    perms = [rng.permutation(y) for _ in range(n_perm)]

    results = {"logreg": np.empty(n_perm), "lgbm": np.empty(n_perm)}
    t0 = time.time()

    with ProcessPoolExecutor(max_workers=10) as ex:
        for model_name in ["logreg", "lgbm"]:
            futures = {
                ex.submit(_run_one_permutation, (folds, perms[i], model_name)): i
                for i in range(n_perm)
            }
            done = 0
            for fut in as_completed(futures):
                i = futures[fut]
                results[model_name][i] = fut.result()
                done += 1
                if done % 100 == 0:
                    elapsed = time.time() - t0
                    print(f"  [{model_name}] {done}/{n_perm} permutations "
                          f"({elapsed:.1f}s elapsed)")

    elapsed = time.time() - t0
    print(f"  Total permutation wall time: {elapsed:.1f}s")

    perm_df = pd.DataFrame({
        "permutation": list(range(n_perm)) * 2,
        "model": ["logreg"] * n_perm + ["lgbm"] * n_perm,
        "auc": np.concatenate([results["logreg"], results["lgbm"]]),
    })

    summary_rows = []
    for model_name in ["logreg", "lgbm"]:
        null = results[model_name]
        obs = observed_auc[model_name]
        p = (1 + np.sum(null >= obs)) / (n_perm + 1)
        summary_rows.append({
            "model": model_name,
            "observed_auc_training_only": round(obs, 4),
            "n_perm": n_perm,
            "null_mean": round(float(null.mean()), 4),
            "null_sd": round(float(null.std()), 4),
            "null_p2.5": round(float(np.quantile(null, 0.025)), 4),
            "null_p97.5": round(float(np.quantile(null, 0.975)), 4),
            "empirical_p": round(p, 6),
            "p_resolution": round(1 / (n_perm + 1), 6),
        })
    return perm_df, pd.DataFrame(summary_rows)


def main():
    print("12_zhu_robustness.py -- Reviewer 4 major statistical request")

    raw_table, ids, y = load_zhu_data()
    print(f"Loaded Zhu 2022 (shanghai2022): n={len(ids)} ({y.sum()} AD / {(y==0).sum()} CN)")

    historical_auc = {"logreg": 0.998, "lgbm": 0.979}  # as reported in the manuscript to date

    # Single-partition training-only AUC (outer_seed=42), for direct comparison
    print("\nSingle-partition strict training-only AUC (outer_seed=42), for reference:")
    folds_ref = precompute_fold_features(raw_table, ids, outer_seed=42)
    training_only_single_auc = {}
    for model_name in ["logreg", "lgbm"]:
        auc = run_nested_cv_given_folds(folds_ref, y, model_name)
        training_only_single_auc[model_name] = round(auc, 4)
        print(f"  {model_name}: {auc:.4f}")

    repeated_df = run_repeated_cv(raw_table, ids, y)
    repeated_df.to_csv(TABLES_DIR / "zhu_repeated_nested_cv.csv", index=False)
    print(f"\n  wrote {TABLES_DIR}/zhu_repeated_nested_cv.csv")

    summary_df = summarize_repeated_cv(repeated_df, historical_auc, training_only_single_auc)
    summary_df.to_csv(TABLES_DIR / "zhu_repeated_nested_cv_summary.csv", index=False)
    print(f"  wrote {TABLES_DIR}/zhu_repeated_nested_cv_summary.csv")
    print(summary_df.to_string(index=False))

    perm_df, perm_summary_df = run_permutation_test(
        raw_table, ids, y, N_PERM_TARGET, training_only_single_auc
    )
    perm_df.to_csv(TABLES_DIR / "zhu_label_permutation.csv", index=False)
    perm_summary_df.to_csv(TABLES_DIR / "zhu_label_permutation_summary.csv", index=False)
    print(f"\n  wrote {TABLES_DIR}/zhu_label_permutation.csv")
    print(f"  wrote {TABLES_DIR}/zhu_label_permutation_summary.csv")
    print(perm_summary_df.to_string(index=False))

    # Supplementary figure: repeated-CV distribution + permutation null with observed marked
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(10, 7.5))
    for col, model_name in enumerate(["logreg", "lgbm"]):
        ax = axes[0, col]
        vals = repeated_df[repeated_df["model"] == model_name]["auc"]
        ax.hist(vals, bins=15, color="#4C72B0", alpha=0.85, edgecolor="white")
        ax.axvline(training_only_single_auc[model_name], color="black", ls="--", lw=1.2,
                   label=f"single-partition = {training_only_single_auc[model_name]:.3f}")
        ax.set_title(f"Repeated nested CV ({model_name})\nfold-assignment stability, n={N_REPEATS}",
                    fontsize=9)
        ax.set_xlabel("Pooled OOF AUC", fontsize=8)
        ax.legend(fontsize=7)

        ax2 = axes[1, col]
        null = perm_df[perm_df["model"] == model_name]["auc"]
        obs = training_only_single_auc[model_name]
        ax2.hist(null, bins=30, color="#888888", alpha=0.85, edgecolor="white")
        ax2.axvline(obs, color="#D62728", lw=1.5, label=f"observed = {obs:.3f}")
        p_row = perm_summary_df[perm_summary_df["model"] == model_name].iloc[0]
        ax2.set_title(f"Label permutation null ({model_name})\nn_perm={N_PERM_TARGET}, "
                      f"p={p_row['empirical_p']:.4g}", fontsize=9)
        ax2.set_xlabel("Pooled OOF AUC (permuted labels)", fontsize=8)
        ax2.legend(fontsize=7)

    plt.tight_layout()
    out_fig = FIGS_DIR / "supp_zhu_robustness.png"
    fig.savefig(out_fig, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\n  saved {out_fig}")

    print("\nDone.")


if __name__ == "__main__":
    main()
