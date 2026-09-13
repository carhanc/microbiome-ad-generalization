# Reviewer 3, Round 2 — Forensic Pre-Audit

**Date:** 2026-09-12
**Branch:** `frontiers-r3-round2-preaudit-2026-09-12`
**Auditor:** Claude Code (read-only forensic pass; no scientific files modified)
**Scope:** Preserve current repository state and document ground truth for five reviewer concerns before any fix is attempted. No analysis code, results, figures, tables, data, or manuscript text was modified as part of this audit.

---

## 1. Git provenance

- **Active branch (created by this audit):** `frontiers-r3-round2-preaudit-2026-09-12`, branched from `main`
- **HEAD SHA before this audit (tip of `main`):** `d412a24511b9d3b4c51e1d3eb80816e864a1fe1e`
  commit subject: "Final pre-submission fixes: Fig4/Fig6 captions, Fig1 label, equation rendering"
- **Remote:** `origin` → `https://github.com/carhanc/microbiome-ad-generalization.git` (fetch + push)
- **Recent commit history (main, at branch point):**
  ```
  d412a24 Final pre-submission fixes: Fig4/Fig6 captions, Fig1 label, equation rendering
  1802e14 Fix Fig3: widen panels, resolve Train Cohort label collisions
  a1d0215 Fix Fig4: widen panels, resolve Panel A label/title overlap
  6547208 Fix figure layout: legend overlaps in Fig2/5, cramping in Fig4, Panel A spacing in Fig6
  1000d7f Final wording: remove genuine, clarify supervised eligibility in S1
  5f6026a Final cleanup: wording, S2 caption, S4 caption, cohort eligibility
  25723d1 Fix text issues: Zhu arithmetic, cohort selection, Harach, SHAP definition, overclaims
  11de9fc Fix Fig4 p-value and Panel A design; fix Fig6A Agathobacter marker
  f447296 Fix formula, Akkermansia SHAP value, overclaiming language, beta-dispersion df
  3b7ade6 Fix Wang et al. 2026 citation author and title
  ```
- **Local uncommitted changes present *before* this audit began:** Yes. Exact list (unchanged by this audit — preserved as-is on the new branch):
  ```
  deleted:    manuscript/manuscript_draft.pdf
  deleted:    results/figures/supp_sensitivity_loco.png
  untracked:  manuscript/manuscript_draft_new.pdf
  untracked:  results/figures/supp_sensitivity_loco (supp_fig_s4).jpg
  untracked:  results/figures/supp_sensitivity_loco (supp_fig_s4).png
  ```
  Verification: `manuscript/manuscript_draft_new.pdf` is byte-identical in size (3,048,542 bytes) to the deleted `manuscript/manuscript_draft.pdf` — this is a local rename, not data loss. Similarly `results/figures/supp_sensitivity_loco (supp_fig_s4).png` (67,077 bytes) matches the deleted `supp_sensitivity_loco.png` exactly, with an additional `.jpg` conversion (161,236 bytes, newer) also present untracked. No content was altered or discarded by this audit; these pre-existing changes are committed as part of the snapshot in Section D below.

---

## 2. Repository inventory relevant to the paper

**Manuscript files** (`manuscript/`):
`draft.md`, `draft.docx`, `supplementary.md`, `supplementary.docx`, `manuscript_draft_new.pdf` (renamed from `manuscript_draft.pdf`), `archive/` (contains `draft_biorxiv.pdf`, `draft_revised.md`, `supplementary.pdf`, stale `.Rhistory`).

**Scripts** (`scripts/`, numbered pipeline + build tooling):
`00_data_acquisition.py`, `01_harmonize_taxonomy.py`, `02_clr_transform.py`, `03_within_cohort_baseline.py`, `04_cross_cohort_generalization.py`, `04b_cohort_classifier.py`, `05_permanova_variance_decomp.R`, `06_batch_correction.R`, `07_corrected_generalization.py`, `08_shap_taxa_comparison.py`, `install_r_packages.R`, `build_pdf.py`, `generate_manuscript_figures.py`, `utils/` (DADA2 R helper).

**Processed input matrices consumed by scripts 02–08** (`data/processed/`, git-tracked):
`unified_genus_matrix.csv` (509×396+cohort), `diagnosis_labels.csv`, `sample_metadata.csv`, `clr_matrix.csv`, `clr_matrix_combatseq.csv`, `clr_matrix_mmuphin.csv`, `corrected_counts_combatseq.csv`, `clr_summary.csv`, `pipeline_summary.csv`, per-cohort `genus_per_cohort/*_genus_relabund.csv`, per-cohort raw `dada2/<cohort>/genus_table.csv` (pre-filter — used in Section 5 below for provenance verification), plus DADA2 logs and QC artifacts.

**results/tables** (21 CSVs, all git-tracked): within/LOCO/pairwise AUC tables, PERMANOVA/betadisper tables, batch-correction comparison tables, SHAP importance/overlap/flip/null tables — full list captured in inventory pass, headline files: `within_cohort_auc.csv`, `loco_auc.csv`, `pairwise_auc.csv`, `auc_drop_summary.csv`, `sensitivity_loco.csv`, `permanova_cohort.csv`, `permanova_cohort_diagnosis.csv`, `betadisper_cohort.csv`, `variance_decomp_summary.csv`, `loco_auc_corrected.csv`, `auc_comparison_table.csv`, `cohort_identity_classifier.csv`, `shap_within_cohort_importance.csv`, `shap_loco_importance.csv`, `shap_taxa_overlap.csv`, `shap_directional_flips.csv`, `jaccard_null.csv`.

**results/figures**: `manuscript_fig1–6.jpg` (submission figures), `supp_fig_s1/s2/s3a/s3b.jpg`, `supp_sensitivity_loco (supp_fig_s4).{png,jpg}`, plus `submission_v1/` and `working/` subdirectories retaining prior-round figure PNGs (pairwise heatmaps, SHAP dot plots, PCoA, betadisper boxplot, etc., generated directly by scripts 04/05/07/08).

**Model outputs** (`results/model_outputs/`): `within_cohort_cv/` (per-cohort OOF prediction CSVs), `cross_cohort/` (LOCO + pairwise prediction CSVs, saved logreg coefficients, saved LGBM feature importances), `corrected/`, `shap/`.

**Large/raw files intentionally absent from Git** (per `.gitignore`, verified via `git check-ignore`):
- `data/raw/` — **55 GB** of per-cohort raw/trimmed FASTQ (5 cohorts) — ignored.
- `data/silva/` — **206 MB** SILVA v138.1 reference training sets — ignored.
- `*.fastq`, `*.fastq.gz`, `data/processed/*.fastq*` — ignored.
- `data/processed/*_ORIG.csv` — ignored.
- `__pycache__/`, `*.pyc`, `.DS_Store`, `.claude/`, `manuscript/~*.docx` — ignored (tooling/OS noise).

67 files under `data/` are git-tracked despite the above — this includes all DADA2 per-cohort outputs (`genus_table.csv`, `asv_table.csv`, `taxonomy_table.csv`, `processing_stats.csv`, quality/error-model PDFs, `.rds` DADA2 objects) and every processed matrix listed above. This is what makes the Section 5 provenance check possible without rerunning DADA2.

---

## 3. Dependency / environment inventory

`environment.yml` declares a conda env `microbiome-ad`: Python 3.11, numpy≥1.26, pandas≥2.0, scipy≥1.11, scikit-learn≥1.3, xgboost≥2.0, lightgbm≥4.0, shap≥0.44, matplotlib≥3.8, seaborn≥0.13, statsmodels≥0.14, pingouin≥0.5, scikit-bio≥0.6, openpyxl≥3.1, h5py≥3.9, jupyter/ipykernel. R packages (installed separately via `scripts/install_r_packages.R`, not conda): vegan, sva (ComBat-seq), MMUPHin, phyloseq, zCompositions.

**Actually-installed versions, checked without modifying the environment** (there are multiple conda envs on this machine; the dedicated `microbiome-ad` env — not `base` — is the one matching the manuscript):

| Package | Manuscript states (§2.6) | Installed in `microbiome-ad` env |
|---|---|---|
| Python | 3.11 | 3.11.15 ✓ |
| scikit-learn | v1.5 | 1.5.2 ✓ |
| lightgbm | v4 | 4.6.0 ✓ |
| shap | v0.51 | 0.51.0 ✓ |
| R | 4.3 | 4.5.1 (newer than stated; not re-verified against original run) |
| vegan | 2.6 | 2.7.1 (newer than stated) |
| sva | not versioned in text | 3.58.0 |
| MMUPHin | v1.15 | 1.24.0 (newer than stated) |

The Python-side versions match the manuscript's stated versions exactly. The R-side installed versions are newer than what the manuscript states (4.3→4.5.1, vegan 2.6→2.7.1, MMUPHin 1.15→1.24.0). This is expected drift from ongoing system updates since the original analysis was run and is **not evidence that the checked-in R results were produced with these newer versions** — no R script was re-run to test this. Flagged for awareness only; not a discrepancy unless Round 2 fixes require rerunning R code, in which case a version pin/reproducibility check would be warranted.

---

## 4. Analysis dependency graph

```
data/raw/<cohort>/fastq(.gz)                                  [55GB, gitignored]
  → cutadapt primer trim → fastq_trimmed/                     [01_harmonize_taxonomy.py Stage A0]
  → DADA2 (scripts/utils/dada2_pipeline.R, per-cohort, independent)
      → data/processed/dada2/<cohort>/genus_table.csv         [git-tracked; pre-filter, per-cohort]
  → apply_synonym_map + drop_uninformative_genera              [01_harmonize_taxonomy.py]
  → apply_prevalence_filter(min_prevalence=0.20, min_cohorts=1) — GLOBAL across all 5 cohorts at once
      → data/processed/unified_genus_matrix.csv (509 × 396)   [git-tracked]
  → multiplicative zero-replacement (δ=0.65, row-wise) + CLR   [02_clr_transform.py]
      → data/processed/clr_matrix.csv (509 × 396)             [git-tracked]
  ├─→ 03_within_cohort_baseline.py (nested nested 10×5 CV, per labeled cohort)
  │     → results/tables/within_cohort_auc.csv
  ├─→ 04_cross_cohort_generalization.py (LOCO + pairwise, 4 labeled cohorts)
  │     → results/tables/loco_auc.csv, pairwise_auc.csv, auc_drop_summary.csv, sensitivity_loco.csv
  │     → results/figures/pairwise_heatmap_*.png, loco_auc_bar.png, supp_sensitivity_loco.png
  ├─→ 04b_cohort_classifier.py (5-class cohort-identity OvR classifier, all 509 samples)
  │     → results/tables/cohort_identity_classifier.csv  (macro-AUC=0.9917, manuscript §3.3/Discussion)
  ├─→ 05_permanova_variance_decomp.R (adonis2 on Aitchison + Bray-Curtis; also consumes
  │     unified_genus_matrix.csv directly for relative abundances)
  │     → results/tables/permanova_cohort.csv, permanova_cohort_diagnosis.csv,
  │       betadisper_cohort.csv, permanova_within_cohort_diagnosis.csv,
  │       variance_decomp_summary.csv  → manuscript Fig 4B/4C, §3.4
  ├─→ 06_batch_correction.R (consumes unified_genus_matrix.csv [raw counts] + clr_matrix.csv;
  │     ComBat-seq fit jointly on all 4 labeled cohorts' true cohort+diagnosis labels;
  │     MMUPHin fit jointly on all 4 labeled cohorts' true cohort+diagnosis labels)
  │     → data/processed/clr_matrix_combatseq.csv, clr_matrix_mmuphin.csv
  │     → 07_corrected_generalization.py (LOCO/pairwise/within-cohort on corrected matrices)
  │           → results/tables/loco_auc_corrected.csv, auc_comparison_table.csv,
  │             within_cohort_auc_corrected.csv → manuscript Fig 5, Table 4, §3.5
  └─→ 08_shap_taxa_comparison.py (consumes clr_matrix.csv + diagnosis_labels.csv only —
        does NOT consume corrected matrices)
        → within-cohort OOF SHAP (10-fold, LinearExplainer/TreeExplainer)
        → LOCO SHAP (train N-1 cohorts, explain held-out)
        → results/tables/shap_within_cohort_importance.csv, shap_loco_importance.csv,
          shap_taxa_overlap.csv, shap_directional_flips.csv, jaccard_null.csv
        → manuscript Fig 6, §3.6

generate_manuscript_figures.py reads the results/tables/*.csv and results/model_outputs/*
above (not raw data) to render manuscript_fig1–6.jpg; build_pdf.py assembles draft.md +
these figures + its own hardcoded FIGURES caption list into the PDF/DOCX.
```

**Every headline manuscript number traces to a specific checked-in file** (spot-verified in Section 10):
- 396 genera / 509 samples → `data/processed/unified_genus_matrix.csv`
- n=401 labeled → `data/processed/diagnosis_labels.csv` (AD+CN rows across 4 labeled cohorts)
- LOCO AUC 0.504–0.813, mean drop 0.176/0.146 → `results/tables/auc_drop_summary.csv`
- Cohort R²=0.172, diagnosis R²=0.014, F=28.03/6.82 → `results/tables/permanova_cohort_diagnosis.csv`
- Cohort R²=0.193 (5-cohort) → `results/tables/variance_decomp_summary.csv`
- Betadisper F=13.93 → `results/tables/betadisper_cohort.csv` (`variance_decomp_summary.csv`: 13.9344)
- Cohort-identity macro-AUC=0.9917 → `results/tables/cohort_identity_classifier.csv`
- Jaccard mean 0.135 (logreg, p=0.0019) / 0.058 (lgbm, p=0.068) → `results/tables/jaccard_null.csv`
- Twelve directional flip taxa (both models) → `results/tables/shap_directional_flips.csv` (re-derived directly, see Section 10)
- ComBat-seq LOCO logreg mean 0.640→0.493, Kazakhstan collapse to 0.308 → `results/tables/auc_comparison_table.csv` / `loco_auc_corrected.csv`

---

## 5. Feature-selection audit — DEFINITIVE FINDING (Reviewer concern #1)

**Files inspected:** `01_harmonize_taxonomy.py`, `02_clr_transform.py`, `03_within_cohort_baseline.py`, `04_cross_cohort_generalization.py`, `04b_cohort_classifier.py`, `08_shap_taxa_comparison.py`.

**Actual implemented rule** (`01_harmonize_taxonomy.py::apply_prevalence_filter`, called from `harmonize_all`):
```python
def apply_prevalence_filter(cohort_tables, min_prevalence=0.20, min_cohorts=1):
    # for each cohort: prevalence = (df > 0).mean(axis=0); passing = prevalence >= min_prevalence
    # retained = genera passing in >= min_cohorts cohorts (default 1)
```
`min_prevalence` defaults to **0.20** (CLI-overridable, but 0.20 is the default used for the checked-in matrix — see reproduction below). This is applied via `harmonize_all(cohorts, min_prevalence)`, where `cohorts` is **all 5 cohorts loaded into one dict before filtering** — i.e., the retained-genus set is computed **globally across every cohort simultaneously**, not per-training-split. There is no train/test distinction anywhere in this code path; it runs once, upstream of every downstream supervised analysis (03, 04, 04b, 08) and both R scripts (05, 06).

**Manuscript's stated rule** (`draft.md` line 66):
> "A genus was retained across the entire study if it was observed at ≥1 count in ≥2 samples from at least one cohort; this threshold was applied per cohort before cross-cohort union."

**These are different rules, and code/manuscript DISAGREE.**

**Which rule actually produced the checked-in 396-genus matrix — verified, not assumed:**
Using the git-tracked pre-filter per-cohort tables (`data/processed/dada2/<cohort>/genus_table.csv`, the exact input `load_genus_table()` reads), both candidate rules were reproduced read-only in a scratch script (no repo files written) with the exact functions imported from `01_harmonize_taxonomy.py`:

| Rule | Genera retained | Exact set-match to checked-in `unified_genus_matrix.csv` columns (396)? |
|---|---|---|
| A — code default: `min_prevalence=0.20`, `min_cohorts=1`, global | **396** | **YES — exact set equality** |
| B — manuscript text: ≥1 count in ≥2 samples, ≥1 cohort | 876 | No |

**Conclusion: the checked-in matrix was generated by the code's actual global 20%-prevalence filter (Rule A), not by the rule stated in the manuscript (Rule B).** This is a confirmed code/manuscript mismatch, not a guess — Rule A reproduces the checked-in 396 genus columns exactly, set-for-set.

**Global vs. training-only:** Confirmed global. The filter runs once in `01_harmonize_taxonomy.py` before any of scripts 03/04/04b/05/06/07/08 execute, using all 5 cohorts' data (including, for every LOCO fold in script 04/04b/08, the cohort that will subsequently be "held out" as the test set). The feature universe itself (which 396 of the ~1,500+ raw genera exist as columns at all) was chosen with the held-out cohort's own prevalence contributing to that decision. Scripts 03/04/07/08 downstream do NOT re-run any feature selection inside their CV loops — they consume the already-fixed 396-column `clr_matrix.csv` as-is. This is exactly the leakage pattern Reviewer 3 describes, and it is real and reproducible, not disputable.

---

## 6. Leakage audit — full preprocessing classification

| Operation | Where | Classification | Notes |
|---|---|---|---|
| Genus retention filter (20% prevalence) | `01_harmonize_taxonomy.py` | **Uses all cohorts (global)** | See Section 5. Feeds every downstream script. |
| Genus synonym merge / uninformative-genus drop | `01_harmonize_taxonomy.py` | Sample-local | Fixed rename/drop rules, no data-driven threshold. |
| Multiplicative zero-replacement (δ=0.65) | `02_clr_transform.py`, mirrored in `06_batch_correction.R` | Sample-local | Row-wise only; confirmed by code inspection (loop over `result.shape[0]`, per-row `total`/`n_zeros`). |
| CLR transform | `02_clr_transform.py` | Sample-local | Row-wise log-ratio to own geometric mean. |
| Model hyperparameter selection (GridSearchCV inner CV) | `03`, `04`, `07`, `08` | Learned from training data only | Inner CV is nested strictly inside the outer train fold / LOCO train-cohorts in all four scripts; confirmed by reading `tune_and_fit`/`run_nested_cv`/`fit_best_*` — `X_test`/held-out cohort is never passed to `GridSearchCV.fit`. |
| LOCO/pairwise model fitting | `04`, `07`, `08` | Learned from training data only | `X_train = vstack(non-held-out cohorts)`; held-out cohort's `X`/`y` only appear in `evaluate_on_test`/SHAP-eval calls. |
| ComBat-seq batch correction | `06_batch_correction.R` | **Uses all cohorts + all diagnosis labels (global, label-informed)** | `ComBat_seq(counts=count_t_lab, batch=batch_vec, group=group_vec, full_mod=TRUE)` is fit once on all 4 labeled cohorts (401 samples) simultaneously, using every sample's own true diagnosis label as the preserved covariate. No train/test split at this stage. |
| MMUPHin batch correction | `06_batch_correction.R` | **Uses all cohorts + all diagnosis labels (global, label-informed)** | `adjust_batch(feature_abd=relab_t_lab, batch="cohort", covariates="diagnosis", data=meta_mmuphin)` — same pattern, fit once on all 401 labeled samples with true labels. |
| Corrected LOCO/pairwise (post-correction) | `07_corrected_generalization.py` | Model fitting: training data only. **Underlying features: label-informed transductive (inherited from `06`)** | The classifier-fitting step itself is properly train/test-separated; the *input features* it fits on were already shaped by a global, label-aware correction step (see Section 7 and Section 9's note on the manuscript's within-cohort-preservation claim). |
| Within-cohort AUC on corrected data | `07_corrected_generalization.py::run_within_cohort_check` | CV structure: training data only. **Feature construction: label-informed (see Section 7)** | Nested CV is correctly implemented (no test-fold leakage into model fitting), but the corrected CLR values being classified were built from a model that used every sample's own true label. |
| SHAP background (within-cohort OOF) | `08_shap_taxa_comparison.py::run_within_cohort_shap` | Learned from training data only | `X_background = X_tr` (same cohort, other 9 folds), `X_eval = X_te`. No leakage, but see Section 8 for the direction-sign caveat this specific choice introduces. |
| SHAP background (LOCO) | `08_shap_taxa_comparison.py::run_loco_shap` | Learned from training data only | `X_background = X_train` (other 3 cohorts), `X_eval = X_test` (held-out cohort). No leakage, but background/eval come from systematically different cohort distributions — see Section 8. |
| Jaccard null baseline | `08_shap_taxa_comparison.py::jaccard_null_distribution` | Sample-local (synthetic) | Draws random genus subsets from `range(n_genera)`; doesn't touch real data at all — its issue is statistical design, not leakage (see Section 8/9). |
| DADA2 ASV calling | `scripts/utils/dada2_pipeline.R`, invoked per-cohort by `01` | Sample-local (per cohort) | Manuscript explicitly states cohorts processed independently "to avoid introducing cross-cohort data leakage during error model estimation" (draft.md line 62) — confirmed by `01_harmonize_taxonomy.py::run_dada2`, which invokes DADA2 once per cohort with no cross-cohort pooling. |

**No supervised step was found to use test-set labels directly during model fitting** (i.e., no classifier ever sees `y_test`). The two label-related leakage findings are (a) the genus feature set itself (Section 5, uses all cohorts' unlabeled prevalence, not diagnosis labels) and (b) batch correction (Section 7, uses all cohorts' diagnosis labels).

---

## 7. Batch-correction audit (Reviewer concern #2)

**Files inspected:** `06_batch_correction.R`, `07_corrected_generalization.py`.

**Where cohort labels enter:** `batch_vec <- meta_lab$cohort` (ComBat-seq) / `batch="cohort"` (MMUPHin), for all 401 labeled samples across all 4 labeled cohorts at once.

**Where diagnosis labels enter:** `group_vec <- ifelse(meta_lab$diagnosis=="AD",1L,2L)` passed to `ComBat_seq(..., group=group_vec, full_mod=TRUE)`; `covariates="diagnosis"` passed to MMUPHin's `adjust_batch(...)`. Both use every sample's **true, real diagnosis label**, for all 4 cohorts simultaneously, in one global model fit.

**Do held-out test samples enter?** Yes. There is no cohort-holdout structure in `06_batch_correction.R` at all — `count_t_lab`/`relab_t_lab` contain all 4 labeled cohorts every time. The "held-out" designation only exists downstream, in `07_corrected_generalization.py`'s LOCO loop, which is applied *after* the correction has already been fit on everyone.

**Do held-out diagnosis labels enter?** Yes, same mechanism — every sample's own true AD/CN label was used as the `group`/`covariates` argument during the correction that produced its own corrected feature vector.

**What ComBat-seq receives:** raw integer genus counts (`count_t_lab`, genera × 401 samples), `batch=cohort`, `group=diagnosis`, `full_mod=TRUE`.
**What MMUPHin receives:** relative abundances (`relab_t_lab`, genera × 401 samples), `batch="cohort"`, `covariates="diagnosis"`.

**Classification of the resulting corrected LOCO AUCs:** **Label-informed transductive**, not prospective and not merely unsupervised-transductive. The correction model had access to every test sample's own cohort identity *and* diagnosis label before the "LOCO" split is even applied in script 07. This is more information than an "unsupervised-transductive" design (which would use test cohort identity/features but not test labels).

**Manuscript disclosure — already present and accurate for the LOCO-transductive point:** `draft.md` explicitly and clearly discloses this exact issue in three places: Section 2.7 (line 129: "both correction methods were fitted jointly on all four labeled cohorts, including the held-out test cohort and its diagnosis labels—a transductive design... corrected LOCO results do not represent a valid prospective external-validation scenario"), Section 3.5 (line 241: "Both methods were applied in a transductive design in which the held-out test cohort and its labels were included during correction parameter estimation"), and a dedicated Limitations bullet (line 335: "Batch correction leakage in LOCO evaluation... should not be extrapolated to prospective deployment"). The code even carries a matching comment (`07_corrected_generalization.py` line 338: `# correction was fit on all cohorts at once — slight leakage in LOCO, noted in paper`). **On this specific point, code and manuscript agree, and the manuscript's caveat is accurate and appropriately prominent** (it appears in Methods, Results, and Limitations).

**Residual, NOT currently disclosed issue — this is the part of Reviewer concern #2 that is still open:** the manuscript's claim at line 253 — *"Critically, within-cohort AUC on corrected data was largely preserved (within ±0.05 of uncorrected)... ruling out the hypothesis that batch correction was simply removing all signal indiscriminately"* — is not adequately caveated. Both `ComBat_seq(..., full_mod=TRUE)` and `adjust_batch(..., covariates="diagnosis")` are, by design, constructed to preserve group-mean structure associated with the covariate they're told to protect (diagnosis), using every sample's own true label in a single global fit. Measuring "within-cohort AUC on the corrected features" afterward is not a clean test of whether *genuine, correction-independent* biological signal survived — part of any preserved separability could be an artifact of the correction algorithm's covariate-preservation mechanism itself having had access to each sample's own label, not evidence that biology survived a label-blind correction. This is distinct from the already-disclosed LOCO-transductive issue (which concerns using test-cohort information in the *cross-cohort* evaluation) and is not addressed by the existing Limitations paragraph, which discusses only the LOCO framing. **This is the more serious, currently-undisclosed half of Reviewer concern #2.**

---

## 8. SHAP audit (Reviewer concern #3 and #4)

**File inspected:** `08_shap_taxa_comparison.py` (in full).

- **Within-cohort OOF procedure:** 10-fold `StratifiedKFold` (or `StratifiedGroupKFold` for shanghai2022, grouped by sample-name prefix), matching `03_within_cohort_baseline.py`'s outer-fold count (`OUTER_FOLDS = 10`, code comment explicitly notes this). For each fold, `fit_best_logreg`/`fit_best_lgbm` refit on the training portion (5-fold inner `GridSearchCV`), then SHAP is computed only for the held-out fold.
- **Logistic-regression SHAP:** `shap.LinearExplainer(model, X_background, feature_perturbation="interventional")` with `X_background = X_train_fold` (i.e., the OOF training partition — explicitly **not** the evaluation set itself). The code contains a correct, self-aware comment explaining *why* background ≠ eval was chosen (using background == eval makes `mean_shap` collapse to ~1e-17 for any linear model under interventional SHAP, since deviations from a set's own mean always sum to zero).
- **LightGBM SHAP:** `shap.TreeExplainer(model).shap_values(X_te, check_additivity=False)`, with a version-compatibility branch for list vs. array return shape.
- **Feature-importance statistic:** `mean_abs_shap = |shap_vals|.mean(axis=0)`, used for ranking/top-N selection.
- **Direction statistic actually used:** `mean_shap = shap_vals.mean(axis=0)` (signed mean, not the model coefficient).
- **Directional-flip criterion:** a taxon is a "flip" if, across cohorts where it appears in that cohort's top-20, `sign(mean_shap)` takes both a positive and a negative value (`directional_analysis` + the flip-set comprehension in `main()`).
- **Are coefficient signs saved per fold?** **No.** `fit_best_logreg` returns `(best_estimator, best_params)`; `model.coef_` is read nowhere in this script. (Contrast with `04_cross_cohort_generalization.py`, which *does* save LOCO logreg coefficients via `_save_logreg_coefs` — but that is a separate script/experiment, not used by the SHAP direction analysis.)
- **Is directional stability formally assessed?** **No.** The only stability-adjacent output is the binary cross-cohort sign-flip tally described above; there is no bootstrap/permutation test of within-cohort sign stability, no confidence measure on the sign itself, and no fold-level variance reported for `mean_shap`.

**Reviewer concern #3, verified as correct:** For `LinearExplainer` with `feature_perturbation="interventional"`, `shap_ij ≈ β_j · (x_ij − mean_j(X_background))`, so `mean_shap_j` for an evaluation set equals `β_j · (mean_j(X_eval) − mean_j(X_background))`. **This is exactly what the code computes.** Its sign therefore depends on both the coefficient sign and the eval-vs-background mean-abundance difference for that genus, not on the coefficient sign alone. For the **within-cohort** OOF analysis, background and eval are different folds of the *same* cohort, so the mean-difference term is close to sampling noise around zero and `sign(mean_shap)` should usually — but not provably — track `sign(β)`. For the **LOCO** analysis (`run_loco_shap`), background = the three *other* cohorts pooled and eval = the entire held-out cohort — here the mean-abundance difference between cohorts is not noise, it is exactly the large, systematic cross-cohort compositional shift that is the paper's central finding (Section 3.4/PERMANOVA). This means LOCO-SHAP "direction" as currently computed is confounded with raw cross-cohort abundance differences, independent of what the coefficient says. **This is a real methodological gap, exactly as the reviewer describes, and it is not currently disclosed or corrected anywhere in code or manuscript.**

Note for the fix-planning stage (not evaluated further under this audit's no-fix mandate): the manuscript's "twelve directional flip taxa for each model" headline claim (Section 3.6) is built from `directional_analysis(within_imp, ...)` — the **within-cohort** OOF table, not the LOCO table — so it is the less-confounded of the two SHAP direction computations, but still uses `mean_shap` rather than a coefficient-based or stability-assessed direction definition.

**Reviewer concern #4, verified as correct:**
- **Observed statistic** (`compute_jaccard_null_baseline`): computed from the full 4×4 pairwise Jaccard matrix among the four real cohorts' top-20 sets, averaged over all off-diagonal entries. Because Jaccard is symmetric, this numerically equals the mean of the **6 unique pairwise overlaps** among the four cohorts (each counted twice in the 12-entry off-diagonal average, which does not change the mean). This matches the manuscript's reported "mean 0.135" / "mean 0.058" statistics exactly (verified against `results/tables/jaccard_null.csv`).
- **Null-generation procedure** (`jaccard_null_distribution`): for each of 10,000 permutation draws, generates **exactly two** independent random top-20 subsets of the 396 genera and computes **one** pairwise Jaccard between them. It does **not** draw four random sets and average six pairwise Jaccards to match the observed statistic's construction.
- **Do they mathematically match?** **No.** The reported null mean (`mean_null_jaccard`, e.g., 0.0266 for logreg) is an unbiased estimate of the correct null's mean by linearity of expectation (any single pairwise Jaccard between two random size-20 sets has the same expectation as the average of six such pairwise Jaccards among four random sets, ≈k/(2n−k) analytically). However, the reported null **95% CI** (`ci_lower_null`/`ci_upper_null`, e.g. [0.0, 0.0811]) reflects the sampling variance of a **single** random pairwise Jaccard, not of a **mean of six** (partially dependent, since the four sets are shared pairwise) random Jaccards — averaging six such quantities reduces variance relative to one draw. The current null is therefore very likely **too wide** (conservative), and the reported `p_value` (comparing the observed mean-of-six statistic against a single-comparison null distribution) does not test the actual sampling distribution of the reported statistic. This does not necessarily overturn the qualitative logreg result (p=0.0019, likely to survive a correctly-calibrated tighter null) but the LightGBM result (p=0.068, "within the random baseline" per the manuscript) is exactly the borderline case most sensitive to the null being mis-specified, and could shift in either direction once refit with a properly matched null. **Confirmed exactly as the reviewer describes: the null generates two random sets per replicate where the observed statistic requires four.**

---

## 9. PERMANOVA audit (Reviewer concern #5)

**Files inspected:** `05_permanova_variance_decomp.R`, `draft.md` (all PERMANOVA-adjacent passages, collected in Section 2/3.4/Discussion).

**What PERMANOVA directly establishes:** `adonis2(dist ~ cohort [+ diagnosis], by="margin")` is a purely descriptive/associational variance-partitioning result: it establishes that cohort-of-origin is associated with a given fraction of total (Aitchison or Bray-Curtis) compositional variance (R²=0.193 across 5 cohorts; marginal R²=0.172 cohort / 0.014 diagnosis in the 4-cohort two-variable model; all p=0.0001 under 9,999 permutations), and that within-cohort dispersion differs significantly by cohort (betadisper F=13.93, p=0.0001). It says nothing, by construction, about *why* a specific downstream algorithm (ComBat-seq/MMUPHin) behaves the way it does.

**What is interpretation/hypothesis in the manuscript, and how it is worded:** The manuscript is largely careful to hedge this. Abstract and Results (draft.md lines 18, 229) use "is consistent with the interpretation that..." rather than a causal claim. The Discussion (line 287) states the mechanistic hypothesis explicitly as conditional reasoning: "A linear batch correction method that targets the cohort dimension of variance **will tend to** also remove disease variance **if** the two sources are not orthogonal — **an assumption that** the empirical worsening of LOCO AUC under correction **suggests** is violated here." This is framed as a plausible mechanism, not a proven causal chain, and is reasonably distinguished from the PERMANOVA measurement itself.

**Manuscript wording that arguably exceeds the analysis:** The same Discussion sentence's lead-in — "The batch correction failure is consistent with the PERMANOVA result... indicates that compositional space is overwhelmingly organized by cohort-of-origin" — blends a direct empirical claim (PERMANOVA does show this) with an implied causal link between "cohort dominates variance" and "therefore correction failed," without a sentence-level flag that the causal direction is inferred, not tested. No experiment in this repository (e.g., an ablation removing the cohort/diagnosis orthogonality assumption, or a simulation under a known ground truth) was run to test the mechanistic hypothesis; it is offered purely as a post-hoc explanation. **Reviewer concern #5 is valid as a request for a sentence-level, unambiguous separation between "PERMANOVA showed X" and "we hypothesize this is part of why Y happened," even though most of the manuscript's language is already reasonably hedged.** This is a wording/framing fix, not a code or results fix — no PERMANOVA number itself is in dispute.

---

## 10. Reproducibility audit

- **All result tables referenced by headline manuscript numbers exist** in `results/tables/` and were checked directly (Section 4 list). None were found missing.
- **Numerical spot-checks performed** (recomputed directly from checked-in CSVs, not from manuscript prose):
  - `variance_decomp_summary.csv`: R²_cohort(5-cohort, Aitchison)=0.193, R²_cohort(marginal,4-cohort)=0.1722, R²_diagnosis(marginal)=0.014, betadisper F=13.9344, p=1e-4 → **matches manuscript exactly** (0.193, 0.172, 0.014, 13.93, p=0.0001).
  - `permanova_cohort_diagnosis.csv`: cohort F=28.03, diagnosis F=6.82, both p=1e-4 → **matches manuscript exactly**.
  - `cohort_identity_classifier.csv`: macro_average auc_ovr=0.9917 → **matches manuscript exactly**.
  - `jaccard_null.csv`: logreg observed_mean=0.135, p=0.0019; lgbm observed_mean=0.0584, p=0.0678 → **matches manuscript's "0.135... p=0.0019" and "0.058... p=0.068" exactly**.
  - `shap_directional_flips.csv`, re-deriving the exact flip-count logic from `08_shap_taxa_comparison.py::main()` (sign-flip across cohorts, ≥2-cohort shared taxa): **logreg = 12 flips, lgbm = 12 flips** → **matches manuscript's "twelve directional flip taxa for each model" exactly**. (An earlier hand-rolled version of this check produced an incorrect 13 for lgbm due to a grouping bug in the ad hoc verification script, not in the underlying data; the figure above uses the actual code's own grouping logic and is the trustworthy value.)
  - `diagnosis_labels.csv`: 509 total rows, 401 AD/CN-labeled rows across the four labeled cohorts → **matches manuscript's n=509 / n=401 exactly**.
  - `auc_drop_summary.csv`: within/LOCO AUC values for all 4 cohorts × 2 models spot-checked against manuscript Table/§3.3 text (e.g., zhuang2018 logreg within=0.6333/LOCO=0.5041) — consistent with the reported 0.504–0.813 LOCO range and 0.176/0.146 mean drops.
- **No manuscript/table/code numerical mismatches were found** among the values checked above.
- **One stale/misleading (but numerically inert) artifact found:** `05_permanova_variance_decomp.R` line 138 prints `cat("PERMANOVA — cohort + diagnosis (n=521)...\n")` — a hardcoded label in a console message that does not match the actual labeled-sample count (401, confirmed via `diagnosis_labels.csv` and the script's own subsequently-computed `nrow(meta_labeled)`). This value is never used in any computation or written to any CSV — it is purely a stale print string, most likely left over from an earlier iteration of the cohort roster. **Severity: cosmetic/low** — does not affect any reported result, but should be corrected for internal consistency during the Round 2 pass.
- **Figures with ambiguous generating source:** none found to be ambiguous — `generate_manuscript_figures.py` reads directly and exclusively from `results/tables/*.csv` and `results/model_outputs/*` (not from raw data), and each `make_figN()` function's source table was traceable.
- **No generated results were found that are not reproducible from committed files** — the one dependency that cannot be regenerated from what's in git is the 55GB raw FASTQ → DADA2 step (by design/gitignore), but its *outputs* (`genus_table.csv` etc.) are committed, which is what made the Section 5 provenance verification possible without rerunning DADA2.
- **No manuscript claims were found unsupported by a checked-in result** among the headline figures audited.

---

## 11. Code/manuscript consistency table

| Manuscript claim/method | Manuscript location | Implemented code location | Agreement? | Severity | Required resolution |
|---|---|---|---|---|---|
| Genus filtering: "≥1 count in ≥2 samples in ≥1 cohort" | draft.md line 66 | `01_harmonize_taxonomy.py::apply_prevalence_filter` (actual: `min_prevalence=0.20`, global) | **NO** | **High** | Either (a) rewrite the manuscript's stated rule to match the code's actual 20%-prevalence global filter, or (b) re-run harmonization with a genuinely training-only / per-fold-safe feature selection and regenerate all downstream results. Reviewer concern #1 requires a *sensitivity analysis*, which is a code+rerun decision outside this audit's scope. |
| 396 genera in unified matrix | draft.md line 66/18 | `data/processed/unified_genus_matrix.csv` | YES (count matches) | — | No action; genus *count* is right, only the *stated rule* producing it is wrong (see row above). |
| Zero replacement δ=0.65, row-wise, 1%-of-total safeguard | draft.md Section 2.3 | `02_clr_transform.py::multiplicative_replacement`, mirrored in `06_batch_correction.R::mult_replace` | YES | — | None. |
| Nested CV (10-fold outer / 5-fold inner) | draft.md Section 2.4 | `03_within_cohort_baseline.py::run_nested_cv` | YES | — | None. |
| LOCO tuning (inner CV on training cohorts only) | draft.md Section 2.5 | `04_cross_cohort_generalization.py::tune_and_fit`/`run_loco` | YES | — | None. |
| AUC bootstrap 95% CIs (1,000 resamples) | draft.md Section 2.4/2.5 | `bootstrap_auc_ci` in `03`/`04`/`07` (`N_BOOTSTRAP=1000`) | YES | — | None. |
| Cohort-exclusion sensitivity analyses (3 configs) | draft.md §3.3/Table, Fig S4 | `04_cross_cohort_generalization.py::SENSITIVITY_CONFIGS`/`run_sensitivity_analysis` | YES | — | None. |
| PERMANOVA n=509 (5-cohort) vs n=401 (4-cohort labeled) | draft.md §3.4 | `05_permanova_variance_decomp.R` (correct sample counts used in each model; only a stale console-print label says "n=521") | YES (computation); print-label bug | Low | Fix the stale `cat("...n=521...")` string in `05_permanova_variance_decomp.R` line 138 to read the actual labeled count. Cosmetic only. |
| Batch correction: transductive design disclosed (LOCO leakage) | draft.md lines 129, 241, 335 | `06_batch_correction.R` (ComBat-seq/MMUPHin fit on all 4 cohorts + true labels jointly) | YES — manuscript accurately discloses this | — | None needed for this specific point. |
| Batch correction: "within-cohort AUC preserved ⇒ signal not indiscriminately removed" | draft.md line 253 | `07_corrected_generalization.py::run_within_cohort_check` on features from `06`'s label-informed global correction | **NO / claim exceeds analysis** | **High** | Manuscript needs to caveat that within-cohort AUC on corrected data cannot cleanly establish "preserved genuine biological signal" because the correction step's covariate-preservation mechanism (`full_mod=TRUE` / `covariates="diagnosis"`) used every sample's own true label globally before this AUC was measured. This is the still-open half of reviewer concern #2 (see Section 7). |
| SHAP direction defined by model coefficient sign | Implied by draft.md §3.6 prose ("association with AD versus cognitive normality") | `08_shap_taxa_comparison.py::compute_importance` uses `mean_shap = shap_vals.mean(axis=0)`, a `LinearExplainer`-background-relative quantity, not `model.coef_` | **NO / definition is not what the prose implies** | **High** | Needs an explicit, coefficient-based (or otherwise background-independent) direction definition plus a documented stability assessment, per reviewer concern #3. Especially urgent for any LOCO-SHAP direction claims, where background/eval cohort mismatch confounds sign with raw abundance differences. |
| Directional flips: "twelve for each model" | draft.md line 305/abstract | `08_shap_taxa_comparison.py` flip-set logic on `shap_within_cohort_importance.csv` → `shap_directional_flips.csv` | YES (numerically, 12/12 verified in Section 10) | — | Numeric claim is correct against current code; however the underlying direction *definition* has the issue in the row above, so the flip *count* could change once direction is redefined. |
| Jaccard null: "significantly above the random baseline of 0.027 [p=0.0019]" (logreg) / "within the random baseline [p=0.068]" (lgbm) | draft.md abstract/§3.6 | `08_shap_taxa_comparison.py::jaccard_null_distribution` draws 2 random sets per replicate; observed statistic is a mean of 6 pairwise Jaccards among 4 real cohorts | **NO / null does not match observed statistic's sampling design** | **High** | Regenerate the null by drawing 4 random top-N sets per replicate and averaging the same 6 pairwise Jaccards the observed statistic uses, per reviewer concern #4. The lgbm p=0.068 borderline result is most likely to change. |
| Cohort identity classifier macro-AUC=0.9917 | draft.md §3.3/Discussion | `04b_cohort_classifier.py` | YES | — | None. |
| PERMANOVA interpretation vs. batch-correction causal explanation | draft.md line 287 (Discussion) | `05_permanova_variance_decomp.R` (descriptive only) | **PARTIAL — wording blurs empirical result and hypothesis** | Medium | Per reviewer concern #5: tighten Discussion wording so the PERMANOVA-established fact ("cohort dominates variance") and the hypothesis ("this is plausibly why correction failed") are unambiguously separated at the sentence level. Most of the manuscript already hedges this reasonably; one Discussion passage needs tightening. |

---

## C. Validation performed (no outputs altered)

- `python -m py_compile` on all 12 `scripts/*.py` files using the project's actual `microbiome-ad` conda env interpreter (Python 3.11.15) — **all compiled cleanly, no syntax errors**.
- `Rscript`-based `parse()` syntax check on all 3 `scripts/*.R` files — **all parsed cleanly, no syntax errors**.
- No table or figure was regenerated. No DADA2, batch-correction, or classifier analysis was rerun. The one computation performed for this audit (Section 5's genus-filter provenance reproduction) was executed from a read-only scratch script outside the repository, importing functions from `01_harmonize_taxonomy.py` without modification and writing no output into the repository; it reads only already-committed, pre-filter per-cohort tables.

---

## D. Snapshot commit

This report, plus the pre-existing (pre-audit) uncommitted working-tree changes listed in Section 1, are committed together on this branch with no other scientific files touched. See commit `chore: snapshot pre-round2 state and add reviewer 3 forensic audit` on `frontiers-r3-round2-preaudit-2026-09-12`.
