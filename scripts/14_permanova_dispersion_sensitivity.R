#!/usr/bin/env Rscript
# PERMANOVA / dispersion robustness (Reviewer 4).
#
# Does NOT reinterpret PERMANOVA as separating technical from biological
# effects -- "cohort" remains a composite study-of-origin variable (technical,
# geographic, recruitment, diagnostic, and dietary differences all bundled
# together), and PERMANOVA/betadisper cannot decompose those sources. This
# script adds two things the original 05_permanova_variance_decomp.R did not:
# (1) multiplicity-adjusted pairwise betadisper comparisons identifying which
# specific cohort pairs drive the overall dispersion difference, and (2) a
# "more technically comparable" (not technically identical) sensitivity
# restricted to the three Chinese paired-end MiSeq cohorts (Zhuang 2018,
# Ling 2020, Zhu 2022), which removes the single-end/NovaSeq (Kazakhstan) and
# cross-country (Kim/KBASE) technical confounds from that specific comparison.

suppressPackageStartupMessages({
  library(vegan)
  library(dplyr)
  library(readr)
})

set.seed(42)
N_PERM <- 9999

args <- commandArgs(trailingOnly = FALSE)
script_flag <- grep("^--file=", args, value = TRUE)
if (length(script_flag)) {
  PROJECT_ROOT <- dirname(dirname(normalizePath(sub("^--file=", "", script_flag))))
} else {
  PROJECT_ROOT <- normalizePath(".")
}
PROCESSED_DIR <- file.path(PROJECT_ROOT, "data", "processed")
TABLES_DIR    <- file.path(PROJECT_ROOT, "results", "tables")
dir.create(TABLES_DIR, showWarnings = FALSE, recursive = TRUE)

cat("14_permanova_dispersion_sensitivity.R\n\n")

clr_df   <- read_csv(file.path(PROCESSED_DIR, "clr_matrix.csv"), show_col_types = FALSE)
diag_df  <- read_csv(file.path(PROCESSED_DIR, "diagnosis_labels.csv"), show_col_types = FALSE)

sample_ids <- clr_df[[1]]
cohort_vec <- clr_df$cohort
genus_cols <- setdiff(names(clr_df), c(names(clr_df)[1], "cohort"))
clr_mat <- as.matrix(clr_df[, genus_cols])
rownames(clr_mat) <- sample_ids

meta <- data.frame(run_id = sample_ids, cohort = cohort_vec, stringsAsFactors = FALSE) %>%
  left_join(diag_df %>% select(run_id, diagnosis), by = "run_id")

# =============================================================================
# 1. Pairwise beta-dispersion post-hoc (all 5 cohorts, Aitchison), with
#    multiplicity-adjusted p-values, identifying which cohort pairs drive the
#    overall dispersion difference reported in the main analysis.
# =============================================================================

cat("1. Pairwise beta-dispersion post-hoc (all 5 cohorts, Aitchison)...\n")

ait_dist <- dist(clr_mat, method = "euclidean")
bd_ait <- betadisper(ait_dist, meta$cohort)
bd_perm_ait <- permutest(bd_ait, permutations = N_PERM, pairwise = TRUE)

cat("  Overall betadisper permutation test:\n")
print(bd_perm_ait)

pw <- bd_perm_ait$pairwise
# permutest$pairwise$permuted holds unadjusted permutation p-values per pair
pw_df <- data.frame(
  comparison   = names(pw$observed),
  F_obs        = unname(pw$observed),
  p_perm       = unname(pw$permuted)
)
pw_df$p_holm <- p.adjust(pw_df$p_perm, method = "holm")
pw_df$p_BH   <- p.adjust(pw_df$p_perm, method = "BH")
pw_df <- pw_df %>% arrange(p_perm)

cat("\n  Pairwise dispersion comparisons (sorted by unadjusted p, with Holm/BH-adjusted):\n")
print(pw_df, row.names = FALSE)

write_csv(pw_df, file.path(TABLES_DIR, "betadisper_pairwise.csv"))
cat(sprintf("\n  Wrote: %s/betadisper_pairwise.csv\n\n", TABLES_DIR))

# =============================================================================
# 2. "More technically comparable" subset: the three Chinese paired-end
#    MiSeq cohorts only (Zhuang 2018, Ling 2020, Zhu 2022). Binary AD/CN only.
#    This removes the single-end/NovaSeq (Kazakhstan) and different-country
#    (Kim/KBASE) technical confounds from this specific comparison -- it is
#    NOT technically identical (library prep, recruitment site, and clinical
#    diagnostic procedure still differ across these three studies).
# =============================================================================

cat("2. Sensitivity: Chinese paired-end MiSeq cohorts only (Zhuang+Ling+Zhu), binary AD/CN...\n")

MISEQ_COHORTS <- c("zhuang2018", "ling2020", "shanghai2022")
mask <- meta$cohort %in% MISEQ_COHORTS & meta$diagnosis %in% c("AD", "CN")
meta_sub <- meta[mask, ]
clr_sub  <- clr_mat[mask, ]

cat(sprintf("  Subset: %d samples (%s)\n", nrow(meta_sub),
            paste(names(table(meta_sub$cohort)), table(meta_sub$cohort), sep="=", collapse=", ")))

ait_sub <- dist(clr_sub, method = "euclidean")

perm_sub <- adonis2(
  ait_sub ~ cohort + diagnosis,
  data         = meta_sub,
  permutations = N_PERM,
  by           = "margin"
)
cat("\n  Marginal PERMANOVA (Aitchison), Chinese-MiSeq subset:\n")
print(perm_sub)

sub_df <- as.data.frame(perm_sub)
sub_df$term <- rownames(sub_df)
sub_df <- sub_df %>%
  rename_with(~ gsub("Pr\\(>F\\)", "p_value", .)) %>%
  rename_with(~ gsub("^F$", "F_stat", .)) %>%
  select(term, Df, SumOfSqs, R2, F_stat, p_value) %>%
  filter(!term %in% c("Total"))
sub_df$subset <- "zhuang2018+ling2020+shanghai2022 (Chinese paired-end MiSeq, binary AD/CN)"

write_csv(sub_df, file.path(TABLES_DIR, "permanova_chinese_miseq_sensitivity.csv"))
cat(sprintf("\n  Wrote: %s/permanova_chinese_miseq_sensitivity.csv\n\n", TABLES_DIR))

bd_sub <- betadisper(ait_sub, meta_sub$cohort)
bd_perm_sub <- permutest(bd_sub, permutations = N_PERM, pairwise = TRUE)
cat("  Beta-dispersion, Chinese-MiSeq subset:\n")
print(bd_perm_sub)

bd_sub_result <- data.frame(
  subset      = "zhuang2018+ling2020+shanghai2022",
  F_stat      = bd_perm_sub$tab$`F`[1],
  df_num      = bd_perm_sub$tab$Df[1],
  df_den      = bd_perm_sub$tab$Df[2],
  p_value     = bd_perm_sub$tab$`Pr(>F)`[1]
)
pw_sub <- bd_perm_sub$pairwise
pw_sub_df <- data.frame(
  subset     = "zhuang2018+ling2020+shanghai2022",
  comparison = names(pw_sub$observed),
  F_obs      = unname(pw_sub$observed),
  p_perm     = unname(pw_sub$permuted)
)
pw_sub_df$p_holm <- p.adjust(pw_sub_df$p_perm, method = "holm")

write_csv(bd_sub_result, file.path(TABLES_DIR, "betadisper_chinese_miseq_sensitivity.csv"))
cat(sprintf("  Wrote: %s/betadisper_chinese_miseq_sensitivity.csv\n", TABLES_DIR))
write_csv(pw_sub_df, file.path(TABLES_DIR, "betadisper_chinese_miseq_pairwise.csv"))
cat(sprintf("  Wrote: %s/betadisper_chinese_miseq_pairwise.csv\n\n", TABLES_DIR))

cat("Summary:\n")
cat(sprintf("  Chinese-MiSeq subset (n=%d): cohort marginal R2=%.4f, diagnosis marginal R2=%.4f\n",
            nrow(meta_sub),
            sub_df$R2[sub_df$term=="cohort"], sub_df$R2[sub_df$term=="diagnosis"]))
cat(sprintf("  Compare to full 4-cohort marginal model: cohort R2=0.172, diagnosis R2=0.014\n"))
cat(sprintf("  Chinese-MiSeq betadisper: F=%.3f, p=%.4f\n", bd_sub_result$F_stat, bd_sub_result$p_value))

cat("\nDone.\n")
