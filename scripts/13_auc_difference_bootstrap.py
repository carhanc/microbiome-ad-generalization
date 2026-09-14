"""Formal within-cohort vs. LOCO AUC difference (Reviewer 4).

Non-overlapping bootstrap CIs are not a formal statistical comparison of two
AUCs. This script performs a PAIRED, diagnosis-stratified, participant-level
bootstrap: for each of the 4 labeled cohorts x 2 models, the within-cohort OOF
prediction and the LOCO (that cohort held out) prediction refer to the SAME
participants. Each bootstrap replicate resamples participants (with
replacement, separately within AD and within CN) ONCE, and uses that SAME
resampled index set to compute both AUC_within and AUC_LOCO -- so the
resulting delta_auc = AUC_within - AUC_LOCO distribution reflects a genuine
paired comparison, not two independently resampled statistics.

Within-cohort OOF prediction files (results/model_outputs/within_cohort_cv/
<cohort>_<model>_oof_predictions.csv) do not carry a run_id column; this
script reconstructs participant identity by replicating the exact data-
loading/filtering sequence from 03_within_cohort_baseline.py and verifies,
by assertion, that the reconstructed diagnosis vector matches the saved
y_true column exactly before trusting the recovered ordering. LOCO
prediction files already carry run_id directly.

This script also regenerates every within-cohort and LOCO AUC 95% CI
currently reported using the same diagnosis-stratified, participant-level
bootstrap (10,000 replicates), replacing whatever CI procedure originally
produced results/tables/within_cohort_auc.csv / loco_auc.csv, so every CI in
the paper is produced by one consistent, explicitly-specified procedure.
"""
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

PROJECT_ROOT = Path(__file__).resolve().parent.parent
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
MODEL_DIR_WC = PROJECT_ROOT / "results" / "model_outputs" / "within_cohort_cv"
MODEL_DIR_CC = PROJECT_ROOT / "results" / "model_outputs" / "cross_cohort"
TABLES_DIR = PROJECT_ROOT / "results" / "tables"
TABLES_DIR.mkdir(parents=True, exist_ok=True)

SUPERVISED_COHORTS = ["zhuang2018", "ling2020", "shanghai2022", "kazakhstan2022"]
RANDOM_STATE = 42
N_BOOT = 10_000


def reconstruct_within_cohort_run_ids():
    """Replicates 03_within_cohort_baseline.py's exact data-loading/filter
    sequence to recover participant run_id for each cohort's OOF predictions,
    verified against the saved y_true column."""
    clr_df = pd.read_csv(PROCESSED_DIR / "clr_matrix.csv", index_col=0)
    diag_df = pd.read_csv(PROCESSED_DIR / "diagnosis_labels.csv")
    clr_reset = clr_df.reset_index().rename(columns={"sample_id": "run_id"})
    merged = clr_reset.merge(diag_df[["run_id", "cohort", "diagnosis"]],
                              on="run_id", suffixes=("_clr", ""))
    if "cohort_clr" in merged.columns:
        merged = merged.drop(columns=["cohort_clr"])
    merged = merged.set_index("run_id")

    run_ids_by_cohort = {}
    for cohort in SUPERVISED_COHORTS:
        cohort_data = merged[merged["cohort"] == cohort].copy()
        if "MCI" in cohort_data["diagnosis"].unique():
            cohort_binary = cohort_data[cohort_data["diagnosis"].isin(["AD", "CN"])].copy()
        else:
            cohort_binary = cohort_data.copy()
        cohort_binary = cohort_binary[cohort_binary["diagnosis"].notna()]
        y_reconstructed = (cohort_binary["diagnosis"] == "AD").astype(int).to_numpy()
        run_ids = np.array(cohort_binary.index)

        for model in ["logreg", "lgbm"]:
            oof = pd.read_csv(MODEL_DIR_WC / f"{cohort}_{model}_oof_predictions.csv")
            assert np.array_equal(oof["y_true"].to_numpy(), y_reconstructed), \
                f"run_id reconstruction failed for {cohort}/{model}"
        run_ids_by_cohort[cohort] = run_ids
    return run_ids_by_cohort


def load_paired_predictions():
    """Returns {(cohort, model): DataFrame[run_id, y_true, y_score_within, y_score_loco]}"""
    run_ids_by_cohort = reconstruct_within_cohort_run_ids()
    paired = {}
    for cohort in SUPERVISED_COHORTS:
        run_ids = run_ids_by_cohort[cohort]
        for model in ["logreg", "lgbm"]:
            wc = pd.read_csv(MODEL_DIR_WC / f"{cohort}_{model}_oof_predictions.csv")
            wc["run_id"] = run_ids
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
    """Pre-generate n_boot sets of resampled positional indices, stratified
    by class, reused identically for both within and LOCO scores per
    replicate (this is what makes the delta a paired comparison)."""
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
    auc_within_boot = np.empty(n_boot)
    auc_loco_boot = np.empty(n_boot)
    for i, idx in enumerate(boots):
        yb = y[idx]
        if len(np.unique(yb)) < 2:
            deltas[i] = np.nan
            auc_within_boot[i] = np.nan
            auc_loco_boot[i] = np.nan
            continue
        a_w = roc_auc_score(yb, score_within[idx])
        a_l = roc_auc_score(yb, score_loco[idx])
        auc_within_boot[i] = a_w
        auc_loco_boot[i] = a_l
        deltas[i] = a_w - a_l
    return deltas, auc_within_boot, auc_loco_boot


def stratified_bootstrap_ci(y, score, n_boot=N_BOOT, seed=RANDOM_STATE):
    """Single-statistic diagnosis-stratified participant bootstrap CI,
    replacing whatever CI procedure was used historically, for consistency."""
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
    print("13_auc_difference_bootstrap.py -- Reviewer 4")
    paired = load_paired_predictions()

    diff_rows = []
    wc_ci_rows = []
    loco_ci_rows = []

    for cohort in SUPERVISED_COHORTS:
        for model in ["logreg", "lgbm"]:
            df = paired[(cohort, model)]
            y = df["y_true_within"].to_numpy()
            score_within = df["y_score_within"].to_numpy()
            score_loco = df["y_score_loco"].to_numpy()

            auc_within_obs = roc_auc_score(y, score_within)
            auc_loco_obs = roc_auc_score(y, score_loco)
            observed_delta = auc_within_obs - auc_loco_obs

            deltas, auc_w_boot, auc_l_boot = paired_delta_bootstrap(y, score_within, score_loco)
            deltas_valid = deltas[~np.isnan(deltas)]

            ci_lo, ci_hi = np.percentile(deltas_valid, [2.5, 97.5])
            prop_le_zero = float((deltas_valid <= 0).mean())

            print(f"\n{cohort} / {model}")
            print(f"  AUC_within={auc_within_obs:.4f}  AUC_LOCO={auc_loco_obs:.4f}  "
                  f"observed delta={observed_delta:.4f}")
            print(f"  bootstrap mean delta={deltas_valid.mean():.4f}  "
                  f"95% CI=[{ci_lo:.4f}, {ci_hi:.4f}]  "
                  f"proportion(delta<=0)={prop_le_zero:.4f}  "
                  f"(n_valid={len(deltas_valid)}/{N_BOOT})")

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
                "n_boot": N_BOOT,
                "n_boot_valid": len(deltas_valid),
            })

            # Regenerated single-statistic CIs (diagnosis-stratified, participant-level)
            wc_lo, wc_hi = stratified_bootstrap_ci(y, score_within)
            wc_ci_rows.append({
                "cohort": cohort, "model": model, "auc_oof": round(auc_within_obs, 4),
                "auc_ci_lo": round(wc_lo, 4), "auc_ci_hi": round(wc_hi, 4),
            })
            loco_lo, loco_hi = stratified_bootstrap_ci(y, score_loco)
            loco_ci_rows.append({
                "cohort": cohort, "model": model, "auc": round(auc_loco_obs, 4),
                "ci_lo": round(loco_lo, 4), "ci_hi": round(loco_hi, 4),
            })

    diff_df = pd.DataFrame(diff_rows)
    diff_df.to_csv(TABLES_DIR / "within_vs_loco_auc_difference_bootstrap.csv", index=False)
    print(f"\n  wrote {TABLES_DIR}/within_vs_loco_auc_difference_bootstrap.csv")

    wc_ci_df = pd.DataFrame(wc_ci_rows)
    wc_ci_df.to_csv(TABLES_DIR / "within_cohort_auc_ci_stratified_bootstrap.csv", index=False)
    print(f"  wrote {TABLES_DIR}/within_cohort_auc_ci_stratified_bootstrap.csv")

    loco_ci_df = pd.DataFrame(loco_ci_rows)
    loco_ci_df.to_csv(TABLES_DIR / "loco_auc_ci_stratified_bootstrap.csv", index=False)
    print(f"  wrote {TABLES_DIR}/loco_auc_ci_stratified_bootstrap.csv")

    print("\nDone.")


if __name__ == "__main__":
    main()
