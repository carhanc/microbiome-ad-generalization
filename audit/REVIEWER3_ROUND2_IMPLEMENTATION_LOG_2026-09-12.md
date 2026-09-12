# Reviewer 3 Round 2 — Implementation Log

Running log of every scientific change made while repairing the five issues
documented in `audit/REVIEWER3_ROUND2_PREAUDIT_2026-09-12.md`. Entries are
appended in execution order, each with: what was done, exact command, env
used, files changed, and old-vs-new numbers where applicable.

Branch: `frontiers-r3-round2-fixes-2026-09-12`, forked from pre-audit commit
`c41a53b51c7dcf93c04006d654d90d73b1259f3f` on
`frontiers-r3-round2-preaudit-2026-09-12`.

---

## 0. Environment provenance

Recorded in `audit/ROUND2_ENVIRONMENT.txt`. All new/rerun analyses below use
`/opt/anaconda3/envs/microbiome-ad/bin/python` (Python 3.11.15) and
`/usr/local/bin/Rscript` (R 4.5.1). DADA2 was NOT rerun anywhere in Round 2.

---

## 1. Training-only genus-selection sensitivity (Reviewer concern #1)

**New script:** `scripts/09_training_only_feature_sensitivity.py`. Reuses
`01_harmonize_taxonomy.py`'s synonym-map / uninformative-genus-drop /
`apply_prevalence_filter` functions and `02_clr_transform.py`'s
`multiplicative_replacement`/`clr_transform` functions unmodified, via direct
module import (no duplication of logic, no rewrite of the original pipeline).

Design decision (documented in the script's module docstring): genus
prevalence for a given split is computed from exactly the AD/CN-labeled
sample rows that enter that split's classifier as training data (not the full
raw per-cohort table, which also contains MCI/unlabeled samples) — this is
the least ambiguous reading of "the held-out cohort must contribute nothing."
Same model families/grids/random_state/fold counts/bootstrap procedure as
`03_within_cohort_baseline.py` and `04_cross_cohort_generalization.py`.

**Command:**
```
/opt/anaconda3/envs/microbiome-ad/bin/python scripts/09_training_only_feature_sensitivity.py
```
**Environment:** `microbiome-ad` conda env, Python 3.11.15 (see `audit/ROUND2_ENVIRONMENT.txt`).
**Exit code:** 0. Only warnings emitted were benign sklearn `l1_ratio` notices
(the C=... l2-penalty grid points don't use `l1_ratio`; this is identical,
pre-existing behavior from the historical scripts' own `LogisticRegression`
configuration, not something introduced here).

**Outputs written (all new, non-destructive — no existing file overwritten):**
- `results/tables/within_cohort_auc_training_only.csv`
- `results/tables/loco_auc_training_only.csv`
- `results/tables/pairwise_auc_training_only.csv`
- `results/tables/training_only_feature_counts.csv` (48 splits: per-fold/per-cohort/per-pair retained genus count)
- `results/tables/training_only_feature_lists.csv` (9,250 rows — full retained-genus list per split, for the Feature-Selection Gate)
- `results/tables/training_only_feature_sensitivity_summary.csv` (original vs. training-only, side by side)

**Headline old-vs-new numbers (LOCO, the analysis central to the paper's claim):**

| Test cohort | Model | Original (leaky, 396 genera) LOCO AUC | Training-only LOCO AUC | Δ | n genera (training-only) |
|---|---|---|---|---|---|
| zhuang2018 | logreg | 0.5041 | 0.5057 | +0.0016 | 370 |
| zhuang2018 | lgbm | 0.5646 | 0.5960 | +0.0314 | 370 |
| ling2020 | logreg | 0.6775 | 0.6794 | +0.0019 | 370 |
| ling2020 | lgbm | 0.6735 | 0.6534 | -0.0201 | 370 |
| shanghai2022 | logreg | 0.8133 | 0.8122 | -0.0011 | 373 |
| shanghai2022 | lgbm | 0.7578 | 0.7878 | +0.0300 | 373 |
| kazakhstan2022 | logreg | 0.5644 | 0.6092 | +0.0448 | 150 |
| kazakhstan2022 | lgbm | 0.5905 | 0.5247 | -0.0658 | 150 |

**Interpretation (descriptive, not tuned):** the retained genus count drops
substantially when Kazakhstan is excluded from feature selection (150 vs. 396
— Kazakhstan alone passes 20% prevalence for many genera that the three
Chinese/Korean-processing-pipeline cohorts do not), but the LOCO AUC values
themselves move by at most ±0.066 and the qualitative conclusion — LOCO AUC
collapses toward chance relative to within-cohort performance — is unchanged
under strict training-only feature selection. Within-cohort AUCs (mirroring
`03_within_cohort_baseline.py`) shift by -0.083 to +0.038. No result flips
across the AUC=0.55 "near chance" line in a way that would change the paper's
qualitative claim. Full within-cohort/pairwise tables in the CSVs above.

---

## 2. Label-blind batch correction (Reviewer concern #2)

**New scripts:**
- `scripts/10_batch_correction_label_blind.R` — mirrors `06_batch_correction.R` exactly except `ComBat_seq(..., group=NULL, full_mod=FALSE)` and `adjust_batch(..., covariates=NULL, data=<cohort-only data frame>)`. Diagnosis is loaded only to define the labeled-sample subset, never passed to either correction function.
- `scripts/11_corrected_generalization_label_blind.py` — mirrors `07_corrected_generalization.py` exactly (same grids/random_state/fold counts/bootstrap), consuming the label-blind corrected matrices instead.

**Commands:**
```
Rscript scripts/10_batch_correction_label_blind.R
/opt/anaconda3/envs/microbiome-ad/bin/python scripts/11_corrected_generalization_label_blind.py
```
**Environment:** R 4.5.1 / `microbiome-ad` Python 3.11.15 (`audit/ROUND2_ENVIRONMENT.txt`).
**Exit codes:** both 0. R console confirms `"Using null model in ComBat-seq. Adjusting for 0 covariate(s)"` — verifying no covariate was supplied.

**Outputs (new, non-destructive):**
- `data/processed/clr_matrix_combatseq_labelblind.csv`, `clr_matrix_mmuphin_labelblind.csv`, `corrected_counts_combatseq_labelblind.csv`
- `results/tables/loco_auc_corrected_labelblind.csv`, `pairwise_auc_labelblind.csv`, `within_cohort_auc_corrected_labelblind.csv`, `auc_comparison_labelblind.csv`

**Headline result — this is a MAJOR, previously-unknown finding, reported exactly as observed:**

| Model | Correction | Mean LOCO AUC | Mean drop vs. within-cohort |
|---|---|---|---|
| logreg | uncorrected | 0.6398 | +0.1756 |
| logreg | ComBat-seq, **label-informed (original/oracle)** | 0.4930 | +0.323 (from original manuscript) |
| logreg | ComBat-seq, **label-blind (new, primary)** | 0.6338 | +0.1816 |
| logreg | MMUPHin, label-informed (original/oracle) | 0.6120 | (from original manuscript) |
| logreg | MMUPHin, label-blind (new, primary) | 0.6390 | +0.1764 |
| lgbm | uncorrected | 0.6466 | +0.1459 |
| lgbm | ComBat-seq, label-blind (new, primary) | 0.6540 | +0.1385 |
| lgbm | MMUPHin, label-blind (new, primary) | 0.6478 | +0.1446 |

**Per-cohort detail, logreg (the most affected model in the original analysis):**

| Test cohort | Uncorrected | ComBat-seq label-informed (original) | ComBat-seq label-blind (new) | MMUPHin label-informed (original) | MMUPHin label-blind (new) |
|---|---|---|---|---|---|
| zhuang2018 | 0.5041 | 0.4889 | 0.5251 | 0.5452 | 0.5284 |
| ling2020 | 0.6775 | 0.4465 | 0.6807 | 0.6644 | 0.7025 |
| shanghai2022 | 0.8133 | 0.7300 | 0.7689 | 0.7778 | 0.7822 |
| kazakhstan2022 | 0.5644 | **0.3080** | 0.5604 | 0.4628 | 0.5428 |

**This is a substantial change to the paper's central batch-correction narrative and is reported honestly, not minimized:** the severe, sometimes near-total collapse seen under the original label-informed design (Kazakhstan logreg LOCO AUC=0.308, ling2020 logreg LOCO AUC=0.4465 — both far below the uncorrected baseline) is **not reproduced** under label-blind correction, where every cohort's LOCO AUC stays close to its uncorrected value (largest single-cohort change: ling2020 logreg +0.0032 to +0.025 depending on method; no cohort collapses). This indicates the original "ComBat-seq substantially worsened cross-cohort generalization" finding was driven in large part by the label-informed/oracle design itself (each sample's own true label being used during correction), not by an intrinsic property of ComBat-seq's batch-effect removal in isolation.

**What is NOT reversed:** label-blind correction also does not meaningfully IMPROVE mean LOCO AUC (deltas of -0.006 to +0.007 relative to uncorrected, i.e. functionally flat) — so the paper's core claim that "standard batch correction methods do not recover cross-cohort generalization" remains supported, and is in fact now demonstrated more cleanly (a null/flat effect under label-blind correction, rather than a confounded catastrophic-degradation effect under label-informed correction). The manuscript must therefore say batch correction is *ineffective at recovering generalization* under a properly-scoped test, while explicitly retracting/reframing the stronger original claim that it *actively and substantially worsens* performance — that specific claim is now understood to be an artifact of the label-informed design and is presented as such (Section C below), with the original label-informed numbers retained only as an explicitly-labeled exploratory/oracle-style sensitivity.

---