#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(sva)
  library(MMUPHin)
  library(vegan)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
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
FIGURES_DIR   <- file.path(PROJECT_ROOT, "results", "figures")
dir.create(FIGURES_DIR, showWarnings = FALSE, recursive = TRUE)

cat("06_batch_correction.R\n\n")

cat("Loading data...\n")

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

cat(sprintf("  Total samples: %d\n", nrow(meta)))
cat(sprintf("  Labeled (AD/CN): %d\n",
            sum(meta$diagnosis %in% c("AD","CN"), na.rm = TRUE)))
cat(sprintf("  Cohorts: %s\n\n",
            paste(sort(unique(meta$cohort)), collapse = ", ")))

LABELED <- c("zhuang2018", "ling2020", "shanghai2022", "kazakhstan2022")
lab_mask <- meta$cohort %in% LABELED & meta$diagnosis %in% c("AD", "CN")
meta_lab      <- meta[lab_mask, ]
count_t_lab   <- count_t[, lab_mask]
relab_lab     <- relab_mat[lab_mask, ]

cat(sprintf("  Labeled subset: %d samples\n", nrow(meta_lab)))
cat(sprintf("  Cohort sizes: %s\n",
            paste(table(meta_lab$cohort), collapse=" | ")))
cat(sprintf("  AD=%d  CN=%d\n\n",
            sum(meta_lab$diagnosis=="AD"), sum(meta_lab$diagnosis=="CN")))

clr_from_matrix <- function(mat) {
  result <- mat
  for (i in seq_len(nrow(mat))) {
    row <- mat[i, ]
    n_zeros <- sum(row == 0)
    total   <- sum(row)
    if (n_zeros == 0 || total == 0) next
    replace_total <- n_zeros * 0.65
    if (replace_total >= total) {
      adj <- total * 0.01 / n_zeros
      row[row == 0] <- adj
      replace_total <- n_zeros * adj
    } else {
      row[row == 0] <- 0.65
    }
    row[row > 0.65 | (mat[i,] > 0)] <-
      row[row > 0.65 | (mat[i,] > 0)] *
      ((total - replace_total) / sum(mat[i, mat[i,] != 0]))
    result[i, ] <- row
  }
  log_mat    <- log(result)
  log_gmean  <- rowMeans(log_mat)
  sweep(log_mat, 1, log_gmean, "-")
}

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

cat("ComBat-seq...\n")

# NB model on counts — takes a few min
group_vec <- ifelse(meta_lab$diagnosis == "AD", 1L, 2L)
batch_vec <- meta_lab$cohort

cat("  Running ComBat_seq (this may take a few minutes)...\n")
tryCatch({
  corrected_counts_cs <- ComBat_seq(
    counts    = count_t_lab,
    batch     = batch_vec,
    group     = group_vec,
    full_mod  = TRUE
  )
  corrected_counts_cs_t <- t(corrected_counts_cs)
  rownames(corrected_counts_cs_t) <- meta_lab$run_id
  colnames(corrected_counts_cs_t) <- rownames(count_t_lab)

  cat(sprintf("  Done. Output range: [%.1f, %.1f]\n",
              min(corrected_counts_cs_t), max(corrected_counts_cs_t)))
  cat(sprintf("  Negative values: %d (will be floored to 0)\n",
              sum(corrected_counts_cs_t < 0)))

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
            file.path(PROCESSED_DIR, "clr_matrix_combatseq.csv"))
  cs_counts_df <- data.frame(
    sample_id = meta_lab$run_id,
    cohort    = meta_lab$cohort,
    diagnosis = meta_lab$diagnosis,
    as.data.frame(corrected_counts_cs_t)
  )
  write_csv(cs_counts_df,
            file.path(PROCESSED_DIR, "corrected_counts_combatseq.csv"))

  cat(sprintf("  Wrote: %s/clr_matrix_combatseq.csv  (%d × %d)\n\n",
              PROCESSED_DIR, nrow(cs_clr), ncol(cs_clr)))
  combatseq_ok <- TRUE

}, error = function(e) {
  cat(sprintf("  ComBat_seq failed: %s\n", conditionMessage(e)))
  cat("  Falling back: skipping ComBat-seq output.\n\n")
  combatseq_ok <<- FALSE
})

cat("MMUPHin::adjust_batch...\n")

relab_t_lab <- t(relab_lab)
rownames(relab_t_lab) <- genus_cols
colnames(relab_t_lab) <- meta_lab$run_id

meta_mmuphin <- data.frame(
  cohort    = meta_lab$cohort,
  diagnosis = meta_lab$diagnosis,
  row.names = meta_lab$run_id
)

cat("  Running adjust_batch...\n")
tryCatch({
  mmuphin_result <- adjust_batch(
    feature_abd = relab_t_lab,
    batch       = "cohort",
    covariates  = "diagnosis",
    data        = meta_mmuphin,
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
            file.path(PROCESSED_DIR, "clr_matrix_mmuphin.csv"))

  cat(sprintf("  Wrote: %s/clr_matrix_mmuphin.csv  (%d × %d)\n\n",
              PROCESSED_DIR, nrow(mm_clr), ncol(mm_clr)))
  mmuphin_ok <- TRUE

}, error = function(e) {
  cat(sprintf("  MMUPHin failed: %s\n", conditionMessage(e)))
  cat("  Skipping MMUPHin output.\n\n")
  mmuphin_ok <<- FALSE
})

cat("Diagnostic PCoA: before vs after...\n")

pcoa_plot <- function(clr_m, meta_sub, title_suffix, colour_var = "cohort") {
  d    <- dist(clr_m, method = "euclidean")
  pcoa <- cmdscale(d, k = 2, eig = TRUE)
  ve   <- pcoa$eig / sum(pcoa$eig[pcoa$eig > 0])
  df   <- data.frame(
    PC1    = pcoa$points[, 1],
    PC2    = pcoa$points[, 2],
    cohort = meta_sub$cohort,
    diag   = meta_sub$diagnosis
  )
  cohort_pal <- c(
    zhuang2018    = "#E41A1C",
    ling2020      = "#FF7F00",
    shanghai2022  = "#984EA3",
    kazakhstan2022= "#377EB8"
  )
  ggplot(df, aes(x = PC1, y = PC2, colour = cohort)) +
    geom_point(alpha = 0.6, size = 1.5) +
    stat_ellipse(aes(group = cohort), type = "norm", level = 0.95,
                 linewidth = 0.7, linetype = "dashed") +
    scale_colour_manual(values = cohort_pal) +
    labs(
      title = title_suffix,
      x = sprintf("PC1 (%.1f%%)", 100 * ve[1]),
      y = sprintf("PC2 (%.1f%%)", 100 * ve[2])
    ) +
    theme_bw(base_size = 10) +
    theme(legend.position = "right",
          panel.grid.minor = element_blank(),
          plot.title = element_text(size = 10))
}

clr_df <- read_csv(file.path(PROCESSED_DIR, "clr_matrix.csv"),
                   show_col_types = FALSE)
clr_before <- as.matrix(
  clr_df[clr_df[[1]] %in% meta_lab$run_id,
          setdiff(names(clr_df), c(names(clr_df)[1], "cohort"))]
)
rownames(clr_before) <- meta_lab$run_id

p_before <- pcoa_plot(clr_before, meta_lab,
                       "Before correction (uncorrected CLR)")

if (exists("combatseq_ok") && combatseq_ok) {
  p_after_cs <- pcoa_plot(cs_clr, meta_lab,
                           "After ComBat-seq (CLR of corrected counts)")
  library(patchwork)
  combined_cs <- p_before + p_after_cs +
    plot_annotation(
      title    = "Aitchison PCoA — 4 labeled cohorts",
      subtitle = "ComBat-seq batch correction: cohort as batch, diagnosis as biological variable",
      theme    = theme(plot.title    = element_text(size = 12),
                       plot.subtitle = element_text(size = 9, colour = "grey40"))
    )
  ggsave(file.path(FIGURES_DIR, "pcoa_combatseq_before_after.png"),
         combined_cs, width = 12, height = 5, dpi = 150)
  cat(sprintf("  Saved: %s/pcoa_combatseq_before_after.png\n", FIGURES_DIR))
}

if (exists("mmuphin_ok") && mmuphin_ok) {
  p_after_mm <- pcoa_plot(mm_clr, meta_lab,
                           "After MMUPHin (CLR of corrected rel. ab.)")
  library(patchwork)
  combined_mm <- p_before + p_after_mm +
    plot_annotation(
      title    = "Aitchison PCoA — 4 labeled cohorts",
      subtitle = "MMUPHin batch correction: cohort as batch, diagnosis as covariate",
      theme    = theme(plot.title    = element_text(size = 12),
                       plot.subtitle = element_text(size = 9, colour = "grey40"))
    )
  ggsave(file.path(FIGURES_DIR, "pcoa_mmuphin_before_after.png"),
         combined_mm, width = 12, height = 5, dpi = 150)
  cat(sprintf("  Saved: %s/pcoa_mmuphin_before_after.png\n", FIGURES_DIR))
}

cat("\nBatch correction complete.\n")
cat("Outputs:\n")
if (exists("combatseq_ok") && combatseq_ok) {
  cat(sprintf("  %s/clr_matrix_combatseq.csv\n", PROCESSED_DIR))
}
if (exists("mmuphin_ok") && mmuphin_ok) {
  cat(sprintf("  %s/clr_matrix_mmuphin.csv\n", PROCESSED_DIR))
}
cat("\nNext step: python scripts/07_corrected_generalization.py\n")
