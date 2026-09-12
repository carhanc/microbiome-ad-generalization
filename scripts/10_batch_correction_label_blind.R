#!/usr/bin/env Rscript
# Label-blind transductive batch correction (Reviewer 3, Round 2, concern #2).
#
# 06_batch_correction.R fits ComBat-seq with group=diagnosis, full_mod=TRUE
# and MMUPHin with covariates="diagnosis" -- both use every sample's own true
# AD/CN label, for all four labeled cohorts jointly, before any LOCO split is
# applied downstream. That is a LABEL-INFORMED transductive design and cannot
# support claims about prospective external-validation performance, nor can
# post-correction within-cohort AUC be used as evidence that genuine biology
# (as opposed to the correction's own label-preservation mechanism) survived.
#
# This script repeats the same two corrections with NO diagnosis information
# supplied at any point: ComBat-seq is called with group=NULL, full_mod=FALSE;
# MMUPHin's `data` argument carries ONLY cohort membership, never diagnosis.
# Cohort identity and every sample's raw feature values are still used (both
# methods estimate a per-batch correction jointly across all samples in a
# batch), so this remains a TRANSDUCTIVE design -- just a label-blind one. It
# is not prospective and not inductive; downstream language must say so.
#
# Everything else (delta=0.65 multiplicative zero-replacement, CLR, cohort
# roster, output shape) is identical to 06_batch_correction.R so the two are
# directly comparable.

suppressPackageStartupMessages({
  library(sva)
  library(MMUPHin)
  library(readr)
  library(dplyr)
})

set.seed(42)

args        <- commandArgs(trailingOnly = FALSE)
script_flag <- grep("^--file=", args, value = TRUE)
if (length(script_flag)) {
  PROJECT_ROOT <- dirname(dirname(normalizePath(sub("^--file=", "", script_flag))))
} else {
  PROJECT_ROOT <- normalizePath(".")
}
PROCESSED_DIR <- file.path(PROJECT_ROOT, "data", "processed")

cat("10_batch_correction_label_blind.R\n\n")
cat("Label-blind transductive batch correction: ComBat-seq(group=NULL, full_mod=FALSE), ")
cat("MMUPHin(covariates=NULL, diagnosis absent from `data`)\n\n")

count_df <- read_csv(file.path(PROCESSED_DIR, "unified_genus_matrix.csv"),
                     show_col_types = FALSE)
diag_df  <- read_csv(file.path(PROCESSED_DIR, "diagnosis_labels.csv"),
                     show_col_types = FALSE)

sample_ids  <- count_df[[1]]
cohort_vec  <- count_df$cohort
genus_cols  <- setdiff(names(count_df), c(names(count_df)[1], "cohort"))

count_mat_raw <- as.matrix(count_df[, genus_cols])
rownames(count_mat_raw) <- sample_ids
count_t <- t(count_mat_raw)

row_sums   <- rowSums(count_mat_raw)
relab_mat  <- sweep(count_mat_raw, 1, pmax(row_sums, 1), "/")

meta <- data.frame(
  run_id    = sample_ids,
  cohort    = cohort_vec,
  stringsAsFactors = FALSE
) %>%
  left_join(diag_df %>% select(run_id, diagnosis), by = "run_id")

LABELED <- c("zhuang2018", "ling2020", "shanghai2022", "kazakhstan2022")
lab_mask <- meta$cohort %in% LABELED & meta$diagnosis %in% c("AD", "CN")
meta_lab      <- meta[lab_mask, ]
count_t_lab   <- count_t[, lab_mask]
relab_lab     <- relab_mat[lab_mask, ]

cat(sprintf("  Labeled subset: %d samples\n", nrow(meta_lab)))
cat(sprintf("  Cohort sizes: %s\n\n",
            paste(table(meta_lab$cohort), collapse=" | ")))

# Diagnosis is loaded above only to define the labeled-sample SUBSET (which
# samples belong to the four supervised cohorts' AD/CN population) -- exactly
# as 06_batch_correction.R does. It is intentionally NEVER passed into
# ComBat_seq's `group` or into MMUPHin's `data`/`covariates` below.

clr_fast <- function(mat) {
  log_mat   <- log(mat)
  log_gmean <- rowMeans(log_mat)
  sweep(log_mat, 1, log_gmean, "-")
}

mult_replace <- function(mat, delta = 0.65) {
  result <- mat
  for (i in seq_len(nrow(mat))) {
    row <- mat[i, ]
    n_z <- sum(row == 0)
    if (n_z == 0) next
    total <- sum(row)
    if (total == 0) { result[i, ] <- delta / ncol(mat); next }
    rep_total <- n_z * delta
    if (rep_total >= total) {
      adj <- total * 0.01 / n_z
      result[i, row == 0] <- adj
      rep_total <- n_z * adj
    } else {
      result[i, row == 0] <- delta
    }
    scale_f <- (total - rep_total) / sum(row[row != 0])
    result[i, row != 0] <- row[row != 0] * scale_f
  }
  result
}

cat("ComBat-seq (label-blind: group=NULL, full_mod=FALSE)...\n")

batch_vec <- meta_lab$cohort

tryCatch({
  corrected_counts_cs <- ComBat_seq(
    counts    = count_t_lab,
    batch     = batch_vec,
    group     = NULL,
    full_mod  = FALSE
  )
  corrected_counts_cs_t <- t(corrected_counts_cs)
  rownames(corrected_counts_cs_t) <- meta_lab$run_id
  colnames(corrected_counts_cs_t) <- rownames(count_t_lab)

  cat(sprintf("  Done. Output range: [%.1f, %.1f]\n",
              min(corrected_counts_cs_t), max(corrected_counts_cs_t)))
  neg_n <- sum(corrected_counts_cs_t < 0)
  cat(sprintf("  Negative values: %d (will be floored to 0)\n", neg_n))
  corrected_counts_cs_t[corrected_counts_cs_t < 0] <- 0

  cat("  Applying multiplicative replacement + CLR...\n")
  cs_replaced <- mult_replace(corrected_counts_cs_t)
  cs_clr      <- clr_fast(cs_replaced)

  max_rowsum <- max(abs(rowSums(cs_clr)))
  cat(sprintf("  CLR sanity check — max |row_sum| = %.2e\n", max_rowsum))
  stopifnot(max_rowsum < 1e-6)

  cs_out_df <- data.frame(
    sample_id = meta_lab$run_id,
    cohort    = meta_lab$cohort,
    diagnosis = meta_lab$diagnosis,
    as.data.frame(cs_clr)
  )
  write_csv(cs_out_df,
            file.path(PROCESSED_DIR, "clr_matrix_combatseq_labelblind.csv"))
  cs_counts_df <- data.frame(
    sample_id = meta_lab$run_id,
    cohort    = meta_lab$cohort,
    diagnosis = meta_lab$diagnosis,
    as.data.frame(corrected_counts_cs_t)
  )
  write_csv(cs_counts_df,
            file.path(PROCESSED_DIR, "corrected_counts_combatseq_labelblind.csv"))

  cat(sprintf("  Wrote: %s/clr_matrix_combatseq_labelblind.csv  (%d x %d)\n\n",
              PROCESSED_DIR, nrow(cs_clr), ncol(cs_clr)))
  combatseq_ok <- TRUE
}, error = function(e) {
  cat(sprintf("  ComBat_seq (label-blind) failed: %s\n", conditionMessage(e)))
  combatseq_ok <<- FALSE
})

cat("MMUPHin::adjust_batch (label-blind: covariates=NULL, diagnosis absent from `data`)...\n")

relab_t_lab <- t(relab_lab)
rownames(relab_t_lab) <- genus_cols
colnames(relab_t_lab) <- meta_lab$run_id

# `data` carries ONLY cohort -- diagnosis is not present in this data frame at
# all, so it cannot be referenced by adjust_batch even by mistake.
meta_mmuphin_lb <- data.frame(
  cohort    = meta_lab$cohort,
  row.names = meta_lab$run_id
)

tryCatch({
  mmuphin_result <- adjust_batch(
    feature_abd = relab_t_lab,
    batch       = "cohort",
    covariates  = NULL,
    data        = meta_mmuphin_lb,
    control     = list(verbose = FALSE)
  )

  corrected_relab_mm <- t(mmuphin_result$feature_abd_adj)
  rownames(corrected_relab_mm) <- meta_lab$run_id
  colnames(corrected_relab_mm) <- genus_cols

  cat(sprintf("  Done. Output range: [%.4f, %.4f]\n",
              min(corrected_relab_mm), max(corrected_relab_mm)))
  neg_count <- sum(corrected_relab_mm < 0)
  if (neg_count > 0) {
    cat(sprintf("  Negative values: %d — flooring to 0\n", neg_count))
    corrected_relab_mm[corrected_relab_mm < 0] <- 0
  }

  cat("  Applying multiplicative replacement + CLR...\n")
  mm_replaced <- mult_replace(corrected_relab_mm)
  mm_clr      <- clr_fast(mm_replaced)

  max_rowsum_mm <- max(abs(rowSums(mm_clr)))
  cat(sprintf("  CLR sanity check — max |row_sum| = %.2e\n", max_rowsum_mm))
  stopifnot(max_rowsum_mm < 1e-6)

  mm_out_df <- data.frame(
    sample_id = meta_lab$run_id,
    cohort    = meta_lab$cohort,
    diagnosis = meta_lab$diagnosis,
    as.data.frame(mm_clr)
  )
  write_csv(mm_out_df,
            file.path(PROCESSED_DIR, "clr_matrix_mmuphin_labelblind.csv"))

  cat(sprintf("  Wrote: %s/clr_matrix_mmuphin_labelblind.csv  (%d x %d)\n\n",
              PROCESSED_DIR, nrow(mm_clr), ncol(mm_clr)))
  mmuphin_ok <- TRUE
}, error = function(e) {
  cat(sprintf("  MMUPHin (label-blind) failed: %s\n", conditionMessage(e)))
  mmuphin_ok <<- FALSE
})

cat("\nLabel-blind batch correction complete.\n")
cat("Outputs:\n")
if (exists("combatseq_ok") && combatseq_ok) {
  cat(sprintf("  %s/clr_matrix_combatseq_labelblind.csv\n", PROCESSED_DIR))
}
if (exists("mmuphin_ok") && mmuphin_ok) {
  cat(sprintf("  %s/clr_matrix_mmuphin_labelblind.csv\n", PROCESSED_DIR))
}
cat("\nNext step: python scripts/11_corrected_generalization_label_blind.py\n")
