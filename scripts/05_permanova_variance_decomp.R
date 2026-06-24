#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(vegan)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggrepel)
})

set.seed(42)
N_PERM <- 9999

args <- commandArgs(trailingOnly = FALSE)
script_flag <- grep("^--file=", args, value = TRUE)
if (length(script_flag)) {
  SCRIPT_PATH  <- normalizePath(sub("^--file=", "", script_flag))
  PROJECT_ROOT <- dirname(dirname(SCRIPT_PATH))
} else {
  PROJECT_ROOT <- normalizePath(".")
}
PROCESSED_DIR <- file.path(PROJECT_ROOT, "data", "processed")
TABLES_DIR    <- file.path(PROJECT_ROOT, "results", "tables")
FIGURES_DIR   <- file.path(PROJECT_ROOT, "results", "figures")

dir.create(TABLES_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(FIGURES_DIR, showWarnings = FALSE, recursive = TRUE)

cat("05_permanova_variance_decomp.R\n\n")

cat("Loading data...\n")

clr_path   <- file.path(PROCESSED_DIR, "clr_matrix.csv")
count_path <- file.path(PROCESSED_DIR, "unified_genus_matrix.csv")
diag_path  <- file.path(PROCESSED_DIR, "diagnosis_labels.csv")

clr_df   <- read_csv(clr_path,   show_col_types = FALSE)
count_df <- read_csv(count_path, show_col_types = FALSE)
diag_df  <- read_csv(diag_path,  show_col_types = FALSE)

cat(sprintf("  CLR matrix:   %d samples × %d columns\n",
            nrow(clr_df), ncol(clr_df)))

sample_ids <- clr_df[[1]]
cohort_vec <- clr_df$cohort
genus_cols <- setdiff(names(clr_df), c(names(clr_df)[1], "cohort"))

clr_mat <- as.matrix(clr_df[, genus_cols])
rownames(clr_mat) <- sample_ids

cat(sprintf("  CLR genus features: %d\n", ncol(clr_mat)))
cat(sprintf("  Cohorts: %s\n", paste(sort(unique(cohort_vec)), collapse = ", ")))

count_genus_cols <- setdiff(names(count_df), c(names(count_df)[1], "cohort"))
count_mat_full <- as.matrix(count_df[, count_genus_cols])
rownames(count_mat_full) <- count_df[[1]]

count_mat <- count_mat_full[rownames(count_mat_full) %in% sample_ids, , drop = FALSE]
count_mat <- count_mat[sample_ids, , drop = FALSE]

row_sums <- rowSums(count_mat)
relab_mat <- sweep(count_mat, 1, pmax(row_sums, 1), "/")

diag_df <- diag_df %>%
  mutate(
    sample_type = case_when(
      cohort == "shanghai2022" ~ substr(sample_name, 1, 1),
      TRUE ~ NA_character_
    )
  )

meta <- data.frame(
  run_id  = sample_ids,
  cohort  = cohort_vec,
  stringsAsFactors = FALSE
) %>%
  left_join(diag_df %>% select(run_id, diagnosis, sample_type),
            by = "run_id")

cat(sprintf("  Samples with diagnosis labels: %d / %d\n",
            sum(!is.na(meta$diagnosis)), nrow(meta)))
cat("\n")

cat("Computing distance matrices...\n")

cat("  Aitchison (Euclidean on CLR)...\n")
# same distance, different name depending which paper you're reading
ait_dist <- dist(clr_mat, method = "euclidean")

cat("  Bray-Curtis (on relative abundances)...\n")
bc_dist  <- vegdist(relab_mat, method = "bray")

cat(sprintf("  Distance matrix size: %d × %d\n\n", nrow(clr_mat), nrow(clr_mat)))

cat("PERMANOVA — cohort (all 5 cohorts)...\n")

perm_cohort_ait <- adonis2(
  ait_dist ~ cohort,
  data        = meta,
  permutations = N_PERM,
  by          = "margin"
)
cat("  Aitchison distance result:\n")
print(perm_cohort_ait)

perm_cohort_bc <- adonis2(
  bc_dist ~ cohort,
  data        = meta,
  permutations = N_PERM,
  by          = "margin"
)
cat("\n  Bray-Curtis distance result:\n")
print(perm_cohort_bc)

extract_adonis2 <- function(a2, label, distance) {
  df <- as.data.frame(a2)
  df$term     <- rownames(df)
  df$model    <- label
  df$distance <- distance
  df <- df %>%
    rename_with(~ gsub("Pr\\(>F\\)", "p_value", .)) %>%
    rename_with(~ gsub("^F$", "F_stat", .)) %>%
    select(model, distance, term, Df, SumOfSqs, R2, F_stat, p_value) %>%
    filter(!term %in% c("Total"))
  df
}

cohort_results <- bind_rows(
  extract_adonis2(perm_cohort_ait, "cohort_only", "aitchison"),
  extract_adonis2(perm_cohort_bc,  "cohort_only", "bray_curtis")
)

write_csv(cohort_results,
          file.path(TABLES_DIR, "permanova_cohort.csv"))
cat(sprintf("\n  Wrote: %s/permanova_cohort.csv\n\n",
            TABLES_DIR))

cat("PERMANOVA — cohort + diagnosis (n=521)...\n")

labeled_mask <- !is.na(meta$diagnosis) & meta$diagnosis %in% c("AD", "CN")
meta_labeled <- meta[labeled_mask, ]
clr_labeled  <- clr_mat[labeled_mask, ]
relab_labeled <- relab_mat[labeled_mask, ]

cat(sprintf("  Labeled cohorts: %d samples\n", nrow(meta_labeled)))
cat(sprintf("  Cohort sizes: %s\n",
            paste(table(meta_labeled$cohort), collapse = " | ")))

ait_labeled <- dist(clr_labeled, method = "euclidean")
bc_labeled  <- vegdist(relab_labeled, method = "bray")

perm_cd_ait <- adonis2(
  ait_labeled ~ cohort + diagnosis,
  data        = meta_labeled,
  permutations = N_PERM,
  by          = "margin"
)
cat("  Aitchison, cohort + diagnosis:\n")
print(perm_cd_ait)

perm_cd_bc <- adonis2(
  bc_labeled ~ cohort + diagnosis,
  data        = meta_labeled,
  permutations = N_PERM,
  by          = "margin"
)
cat("\n  Bray-Curtis, cohort + diagnosis:\n")
print(perm_cd_bc)

cd_results <- bind_rows(
  extract_adonis2(perm_cd_ait, "cohort+diagnosis", "aitchison"),
  extract_adonis2(perm_cd_bc,  "cohort+diagnosis", "bray_curtis")
)
write_csv(cd_results,
          file.path(TABLES_DIR, "permanova_cohort_diagnosis.csv"))
cat(sprintf("\n  Wrote: %s/permanova_cohort_diagnosis.csv\n\n",
            TABLES_DIR))

cat("Beta-dispersion (betadisper)...\n")

bd_ait <- betadisper(ait_dist, meta$cohort)
bd_perm_ait <- permutest(bd_ait, permutations = N_PERM, pairwise = TRUE)

cat("  Aitchison — betadisper permutation test:\n")
print(bd_perm_ait)

cat("\n  Per-cohort median distance to centroid (Aitchison):\n")
bd_distances <- data.frame(
  cohort   = as.character(bd_ait$group),
  distance = bd_ait$distances
)
bd_summary <- bd_distances %>%
  group_by(cohort) %>%
  summarise(
    median_dist = median(distance),
    mean_dist   = mean(distance),
    sd_dist     = sd(distance),
    n           = n(),
    .groups = "drop"
  )
print(bd_summary)

bd_result_df <- data.frame(
  distance    = "aitchison",
  F_stat      = bd_perm_ait$tab$`F`[1],
  df_num      = bd_perm_ait$tab$Df[1],
  df_den      = bd_perm_ait$tab$Df[2],
  p_value     = bd_perm_ait$tab$`Pr(>F)`[1],
  n_perm      = N_PERM,
  interpretation = ifelse(
    bd_perm_ait$tab$`Pr(>F)`[1] < 0.05,
    "Dispersions differ — PERMANOVA result partially reflects heterogeneity",
    "Dispersions homogeneous — PERMANOVA centroid interpretation is valid"
  )
)

pw_bd <- bd_perm_ait$pairwise
if (!is.null(pw_bd)) {
  pw_df <- data.frame(
    comparison = names(pw_bd$observed),
    F_obs      = unname(pw_bd$observed),
    p_perm     = unname(pw_bd$permuted)
  )
  cat("\n  Pairwise betadisper p-values (Aitchison):\n")
  print(pw_df)
} else {
  pw_df <- NULL
}

write_csv(bd_result_df,
          file.path(TABLES_DIR, "betadisper_cohort.csv"))
cat(sprintf("\n  Wrote: %s/betadisper_cohort.csv\n\n", TABLES_DIR))

cat("Within-cohort PERMANOVA (diagnosis)...\n")

LABELED_COHORTS <- c("zhuang2018", "ling2020", "shanghai2022", "kazakhstan2022")

within_results <- list()

for (cohort_name in LABELED_COHORTS) {
  cohort_mask <- meta_labeled$cohort == cohort_name
  sub_clr     <- clr_labeled[cohort_mask, ]
  sub_meta    <- meta_labeled[cohort_mask, ]

  cat(sprintf("  %s  (n=%d, %d AD / %d CN)\n",
              cohort_name,
              nrow(sub_meta),
              sum(sub_meta$diagnosis == "AD"),
              sum(sub_meta$diagnosis == "CN")))

  sub_dist <- dist(sub_clr, method = "euclidean")

  perm_within <- adonis2(
    sub_dist ~ diagnosis,
    data         = sub_meta,
    permutations = N_PERM,
    by           = "margin"
  )
  print(perm_within)

  df_within <- as.data.frame(perm_within)
  df_within$cohort <- cohort_name
  df_within$term   <- rownames(df_within)
  within_results[[cohort_name]] <- df_within
  cat("\n")
}

within_df <- bind_rows(within_results) %>%
  rename_with(~ gsub("Pr\\(>F\\)", "p_value", .)) %>%
  rename_with(~ gsub("^F$", "F_stat", .)) %>%
  filter(!term %in% c("Total")) %>%
  select(cohort, term, Df, SumOfSqs, R2, F_stat, p_value)

write_csv(within_df,
          file.path(TABLES_DIR, "permanova_within_cohort_diagnosis.csv"))
cat(sprintf("  Wrote: %s/permanova_within_cohort_diagnosis.csv\n\n",
            TABLES_DIR))

cat("Variance decomposition summary...\n")

r2_cohort_ait  <- cohort_results %>%
  filter(distance == "aitchison", term == "cohort") %>%
  pull(R2)
p_cohort_ait   <- cohort_results %>%
  filter(distance == "aitchison", term == "cohort") %>%
  pull(p_value)

r2_cohort_bc   <- cohort_results %>%
  filter(distance == "bray_curtis", term == "cohort") %>%
  pull(R2)

r2_diag_ait    <- cd_results %>%
  filter(distance == "aitchison", term == "diagnosis") %>%
  pull(R2)
r2_cohort_adj_ait <- cd_results %>%
  filter(distance == "aitchison", term == "cohort") %>%
  pull(R2)

within_r2 <- within_df %>%
  filter(term == "diagnosis") %>%
  select(cohort, R2_within_diagnosis = R2, p_value_diagnosis = p_value)

summary_df <- data.frame(
  metric            = c(
    "R2_cohort_aitchison_all5",
    "R2_cohort_bray_curtis_all5",
    "R2_cohort_aitchison_labeled4_marginal",
    "R2_diagnosis_aitchison_labeled4_marginal",
    "betadisper_F",
    "betadisper_p"
  ),
  value             = c(
    round(r2_cohort_ait, 4),
    round(r2_cohort_bc, 4),
    round(r2_cohort_adj_ait, 4),
    round(r2_diag_ait, 4),
    round(bd_result_df$F_stat, 4),
    round(bd_result_df$p_value, 4)
  ),
  interpretation    = c(
    "Fraction of total Aitchison variance explained by cohort (5 cohorts)",
    "Fraction of total Bray-Curtis variance explained by cohort (5 cohorts)",
    "Cohort R2 after adjusting for diagnosis (4 labeled cohorts)",
    "Diagnosis R2 after adjusting for cohort (4 labeled cohorts)",
    "Beta-dispersion F-statistic (heterogeneity of within-group spread)",
    "Beta-dispersion p-value"
  )
)

cat("\n  Key variance decomposition metrics:\n")
print(summary_df, row.names = FALSE)

write_csv(summary_df,
          file.path(TABLES_DIR, "variance_decomp_summary.csv"))
cat(sprintf("\n  Wrote: %s/variance_decomp_summary.csv\n\n",
            TABLES_DIR))

cat("Generating PCoA figures...\n")

pcoa_all <- cmdscale(ait_dist, k = 4, eig = TRUE)

var_exp <- pcoa_all$eig / sum(pcoa_all$eig[pcoa_all$eig > 0])
pct1 <- round(100 * var_exp[1], 1)
pct2 <- round(100 * var_exp[2], 1)

pcoa_df <- data.frame(
  run_id = sample_ids,
  PC1    = pcoa_all$points[, 1],
  PC2    = pcoa_all$points[, 2],
  cohort = meta$cohort,
  diagnosis = meta$diagnosis
)

cohort_labels <- c(
  "zhuang2018"    = "Zhuang 2018 (China)",
  "ling2020"      = "Ling 2020 (China)",
  "shanghai2022"  = "Shanghai 2022 (China)",
  "kbase2022"     = "Kim/KBASE 2022 (South Korea)",
  "kazakhstan2022"= "Kazakhstan 2022 (Kazakhstan)"
)

cohort_colours <- c(
  "zhuang2018"    = "#E41A1C",
  "ling2020"      = "#FF7F00",
  "shanghai2022"  = "#984EA3",
  "kbase2022"     = "#4DAF4A",
  "kazakhstan2022"= "#377EB8"
)

p_cohort <- ggplot(pcoa_df, aes(x = PC1, y = PC2, colour = cohort)) +
  geom_point(alpha = 0.65, size = 1.6) +
  stat_ellipse(aes(group = cohort), type = "norm", level = 0.95,
               linewidth = 0.8, linetype = "dashed") +
  scale_colour_manual(values = cohort_colours,
                      labels = cohort_labels,
                      name   = "Cohort") +
  labs(
    title   = "PCoA of Aitchison Distance — All 5 Cohorts",
    subtitle = sprintf(
      "Cohort explains R²=%.1f%% of total compositional variance (PERMANOVA, p<0.001)\nPC1=%.1f%% | PC2=%.1f%% of total variance",
      100 * r2_cohort_ait, pct1, pct2
    ),
    x = sprintf("PC1 (%.1f%%)", pct1),
    y = sprintf("PC2 (%.1f%%)", pct2)
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "right",
    legend.text     = element_text(size = 8),
    plot.subtitle   = element_text(size = 8, colour = "grey40"),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(FIGURES_DIR, "pcoa_by_cohort.png"),
       p_cohort, width = 8, height = 5.5, dpi = 150)
cat(sprintf("  Saved: %s/pcoa_by_cohort.png\n", FIGURES_DIR))

pcoa_labeled <- pcoa_df %>%
  filter(!is.na(diagnosis), diagnosis %in% c("AD", "CN"))

diag_colours <- c("AD" = "#D62728", "CN" = "#1F77B4")
cohort_shapes <- c(
  "zhuang2018"    = 16,
  "ling2020"      = 17,
  "shanghai2022"  = 15,
  "kazakhstan2022"= 18
)

p_diag <- ggplot(pcoa_labeled,
                 aes(x = PC1, y = PC2,
                     colour = diagnosis, shape = cohort)) +
  geom_point(alpha = 0.7, size = 2.0) +
  scale_colour_manual(values = diag_colours,
                      name   = "Diagnosis") +
  scale_shape_manual(values = cohort_shapes,
                     labels = cohort_labels[names(cohort_shapes)],
                     name   = "Cohort") +
  labs(
    title    = "PCoA of Aitchison Distance — Diagnosis Labels",
    subtitle = sprintf(
      "Diagnosis explains R²=%.1f%% of variance (adjusted for cohort, PERMANOVA)\nCohort explains R²=%.1f%% (adjusted for diagnosis)",
      100 * r2_diag_ait, 100 * r2_cohort_adj_ait
    ),
    x = sprintf("PC1 (%.1f%%)", pct1),
    y = sprintf("PC2 (%.1f%%)", pct2)
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position  = "right",
    legend.text      = element_text(size = 8),
    plot.subtitle    = element_text(size = 8, colour = "grey40"),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(FIGURES_DIR, "pcoa_by_diagnosis.png"),
       p_diag, width = 8, height = 5.5, dpi = 150)
cat(sprintf("  Saved: %s/pcoa_by_diagnosis.png\n", FIGURES_DIR))

bd_plot_df <- bd_distances %>%
  left_join(data.frame(
    cohort = names(cohort_labels),
    cohort_label = unname(cohort_labels)
  ), by = "cohort")

p_bd <- ggplot(bd_plot_df, aes(x = reorder(cohort_label, distance, median),
                                y = distance, fill = cohort)) +
  geom_boxplot(outlier.size = 0.8, alpha = 0.75) +
  scale_fill_manual(values = cohort_colours, guide = "none") +
  coord_flip() +
  labs(
    title    = "Within-Group Beta-Dispersion (Aitchison)",
    subtitle = sprintf(
      "Distance from each sample to its cohort centroid\nbetadisper F=%.2f, p=%.4f",
      bd_result_df$F_stat, bd_result_df$p_value
    ),
    x = NULL,
    y = "Distance to cohort centroid (Aitchison)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.subtitle    = element_text(size = 8, colour = "grey40"),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(FIGURES_DIR, "betadisper_boxplot.png"),
       p_bd, width = 7, height = 4, dpi = 150)
cat(sprintf("  Saved: %s/betadisper_boxplot.png\n\n", FIGURES_DIR))

cat("Key findings:\n\n")

cat(sprintf(
  "  Cohort effect (Aitchison, all 5 cohorts):\n    R²=%.4f (%.1f%%)  p=%.4f\n\n",
  r2_cohort_ait, 100*r2_cohort_ait, p_cohort_ait
))

cat(sprintf(
  "  Cohort effect (Bray-Curtis, all 5 cohorts):\n    R²=%.4f (%.1f%%)\n\n",
  r2_cohort_bc, 100*r2_cohort_bc
))

cat(sprintf(
  "  Cohort + diagnosis model (4 labeled cohorts, Aitchison):\n    Cohort R²=%.4f (%.1f%%) [marginal]\n    Diagnosis R²=%.4f (%.1f%%) [marginal]\n\n",
  r2_cohort_adj_ait, 100*r2_cohort_adj_ait,
  r2_diag_ait, 100*r2_diag_ait
))

cat("  Within-cohort diagnosis R² (Aitchison, per cohort):\n")
print(within_r2, row.names = FALSE)

cat(sprintf(
  "\n  Beta-dispersion: F=%.3f, p=%.4f\n    %s\n",
  bd_result_df$F_stat, bd_result_df$p_value,
  bd_result_df$interpretation
))

cat("\nPhase 4 complete.\n")
cat("Outputs:\n")
cat(sprintf("  %s/permanova_cohort.csv\n", TABLES_DIR))
cat(sprintf("  %s/permanova_cohort_diagnosis.csv\n", TABLES_DIR))
cat(sprintf("  %s/betadisper_cohort.csv\n", TABLES_DIR))
cat(sprintf("  %s/permanova_within_cohort_diagnosis.csv\n", TABLES_DIR))
cat(sprintf("  %s/variance_decomp_summary.csv\n", TABLES_DIR))
cat(sprintf("  %s/pcoa_by_cohort.png\n", FIGURES_DIR))
cat(sprintf("  %s/pcoa_by_diagnosis.png\n", FIGURES_DIR))
cat(sprintf("  %s/betadisper_boxplot.png\n", FIGURES_DIR))
cat("\nNext step: Rscript scripts/06_batch_correction.R\n")