"""Formal within-cohort vs. LOCO AUC difference under the strict training-only
pipeline (Reviewer 5, Point 5) -- now the PRIMARY classifier-performance
pipeline. Reuses the exact paired, diagnosis-stratified, participant-level
bootstrap procedure of 13_auc_difference_bootstrap.py (Reviewer 4), applied to
16_training_only_oof_predictions.py's per-sample predictions instead of the
original global-396 predictions. No new statistical test is introduced; only
the input predictions change.

Unlike 13_auc_difference_bootstrap.py, no run_id reconstruction is needed --
16_training_only_oof_predictions.py already saves run_id directly in both the
within-cohort and LOCO prediction files.
"""
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

PROJECT_ROOT = Path(__file__).resolve().parent.parent
MODEL_DIR_WC = PROJECT_ROOT / "results" / "model_outputs" / "within_cohort_cv_training_only"
MODEL_DIR_CC = PROJECT_ROOT / "results" / "model_outputs" / "cross_cohort_training_only"
TABLES_DIR = PROJECT_ROOT / "results" / "tables"

SUPERVISED_COHORTS = ["zhuang2018", "ling2020", "shanghai2022", "kazakhstan2022"]
RANDOM_STATE = 42
N_BOOT = 10_000


def load_paired_predictions():
    paired = {}
    for cohort in SUPERVISED_COHORTS:
        for model in ["logreg", "lgbm"]:
            wc = pd.read_csv(MODEL_DIR_WC / f"{cohort}_{model}_oof_predictions.csv")
            loco = pd.read_csv(MODEL_DIR_CC / f"loco_{cohort}_{model}_predictions.csv")
            merged = wc.merge(loco[["run_id", "y_true", "y_score"]], on="run_id",
                               suffixes=("_within", "_loco"))
            assert len(merged) == len(wc) == len(loco), \
                f"participant count mismatch for {cohort}/{model}"
            assert np.array_equal(merged["y_true_within"].to_numpy(),
                                   merged["y_true_loco"].to_numpy()), \
                f"diagnosis mismatch between within/LOCO files for {cohort}/{model}"
            paired[(cohort, model)] = merged
    return paired


def stratified_bootstrap_indices(y, rng, n_boot):
    ad_idx = np.where(y == 1)[0]
    cn_idx = np.where(y == 0)[0]
    boots = []
    for _ in range(n_boot):
        ad_samp = rng.choice(ad_idx, size=len(ad_idx), replace=True)
        cn_samp = rng.choice(cn_idx, size=len(cn_idx), replace=True)
        boots.append(np.concatenate([ad_samp, cn_samp]))
    return boots


def paired_delta_bootstrap(y, score_within, score_loco, n_boot=N_BOOT, seed=RANDOM_STATE):
    rng = np.random.default_rng(seed)
    boots = stratified_bootstrap_indices(y, rng, n_boot)
    deltas = np.empty(n_boot)
    for i, idx in enumerate(boots):
        yb = y[idx]
        if len(np.unique(yb)) < 2:
            deltas[i] = np.nan
            continue
        a_w = roc_auc_score(yb, score_within[idx])
        a_l = roc_auc_score(yb, score_loco[idx])
        deltas[i] = a_w - a_l
    return deltas


def stratified_bootstrap_ci(y, score, n_boot=N_BOOT, seed=RANDOM_STATE):
    rng = np.random.default_rng(seed)
    boots = stratified_bootstrap_indices(y, rng, n_boot)
    aucs = []
    for idx in boots:
        yb = y[idx]
        if len(np.unique(yb)) < 2:
            continue
        aucs.append(roc_auc_score(yb, score[idx]))
    aucs = np.array(aucs)
    return float(np.percentile(aucs, 2.5)), float(np.percentile(aucs, 97.5))


def main():
    print("17_training_only_auc_difference_bootstrap.py -- Reviewer 5, Point 5 (primary pipeline)")
    paired = load_paired_predictions()

    diff_rows, wc_ci_rows, loco_ci_rows = [], [], []

    for cohort in SUPERVISED_COHORTS:
        for model in ["logreg", "lgbm"]:
            df = paired[(cohort, model)]
            y = df["y_true_within"].to_numpy()
            score_within = df["y_score_within"].to_numpy()
            score_loco = df["y_score_loco"].to_numpy()

            auc_within_obs = roc_auc_score(y, score_within)
            auc_loco_obs = roc_auc_score(y, score_loco)
            observed_delta = auc_within_obs - auc_loco_obs

            deltas = paired_delta_bootstrap(y, score_within, score_loco)
            deltas_valid = deltas[~np.isnan(deltas)]
            ci_lo, ci_hi = np.percentile(deltas_valid, [2.5, 97.5])
            prop_le_zero = float((deltas_valid <= 0).mean())

            print(f"\n{cohort} / {model}")
            print(f"  AUC_within={auc_within_obs:.4f}  AUC_LOCO={auc_loco_obs:.4f}  "
                  f"observed delta={observed_delta:.4f}")
            print(f"  95% CI=[{ci_lo:.4f}, {ci_hi:.4f}]  "
                  f"proportion(delta<=0)={prop_le_zero:.4f}  "
                  f"significant={'YES' if (ci_lo > 0 or ci_hi < 0) else 'no'}")

            diff_rows.append({
                "cohort": cohort, "model": model,
                "n_participants": len(y), "n_AD": int(y.sum()), "n_CN": int((y == 0).sum()),
                "auc_within_observed": round(auc_within_obs, 4),
                "auc_loco_observed": round(auc_loco_obs, 4),
                "observed_delta_auc": round(observed_delta, 4),
                "bootstrap_mean_delta": round(float(deltas_valid.mean()), 4),
                "bootstrap_delta_ci_lo": round(float(ci_lo), 4),
                "bootstrap_delta_ci_hi": round(float(ci_hi), 4),
                "proportion_bootstrap_delta_le_0": round(prop_le_zero, 4),
                "n_boot": N_BOOT, "n_boot_valid": len(deltas_valid),
            })

            wc_lo, wc_hi = stratified_bootstrap_ci(y, score_within)
            wc_ci_rows.append({"cohort": cohort, "model": model,
                                "auc_oof_training_only": round(auc_within_obs, 4),
                                "auc_ci_lo": round(wc_lo, 4), "auc_ci_hi": round(wc_hi, 4)})
            loco_lo, loco_hi = stratified_bootstrap_ci(y, score_loco)
            loco_ci_rows.append({"cohort": cohort, "model": model,
                                  "auc_training_only": round(auc_loco_obs, 4),
                                  "ci_lo": round(loco_lo, 4), "ci_hi": round(loco_hi, 4)})

    diff_df = pd.DataFrame(diff_rows)
    diff_df.to_csv(TABLES_DIR / "within_vs_loco_auc_difference_bootstrap_training_only.csv", index=False)
    print(f"\n  wrote {TABLES_DIR}/within_vs_loco_auc_difference_bootstrap_training_only.csv")

    wc_ci_df = pd.DataFrame(wc_ci_rows)
    wc_ci_df.to_csv(TABLES_DIR / "within_cohort_auc_ci_stratified_bootstrap_training_only.csv", index=False)
    print(f"  wrote {TABLES_DIR}/within_cohort_auc_ci_stratified_bootstrap_training_only.csv")

    loco_ci_df = pd.DataFrame(loco_ci_rows)
    loco_ci_df.to_csv(TABLES_DIR / "loco_auc_ci_stratified_bootstrap_training_only.csv", index=False)
    print(f"  wrote {TABLES_DIR}/loco_auc_ci_stratified_bootstrap_training_only.csv")

    n_sig = sum(1 for r in diff_rows if r["bootstrap_delta_ci_lo"] > 0 or r["bootstrap_delta_ci_hi"] < 0)
    print(f"\n  {n_sig}/8 cohort-model pairs have a 95% CI on the delta excluding zero.")
    print("\nDone.")


if __name__ == "__main__":
    main()
