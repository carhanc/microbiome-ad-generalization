#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dada2)
  library(optparse)
})

option_list <- list(
  make_option("--cohort",       type = "character"),
  make_option("--fastq_dir",    type = "character"),
  make_option("--out_dir",      type = "character"),
  make_option("--mode",         type = "character", default = "paired"),
  make_option("--trunc_fwd",    type = "integer",   default = 240L),
  make_option("--trunc_rev",    type = "integer",   default = 200L),
  make_option("--silva_path",   type = "character"),
  make_option("--threads",      type = "integer",   default = 4L),
  make_option("--inspect_only", type = "logical",   default = FALSE)
)

opt <- parse_args(OptionParser(option_list = option_list))

stopifnot(
  !is.null(opt$cohort),
  !is.null(opt$fastq_dir),
  !is.null(opt$out_dir),
  opt$mode %in% c("paired", "single")
)

if (!isTRUE(opt$inspect_only) && is.null(opt$silva_path)) {
  stop("--silva_path is required unless --inspect_only TRUE")
}

dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)
filt_dir <- file.path(opt$out_dir, "filtered")

cat(sprintf("\nCohort: %s | Mode: %s | Threads: %d\n",
            opt$cohort, opt$mode, opt$threads))

fwd_files <- sort(list.files(opt$fastq_dir, pattern = "_1\\.fastq(\\.gz)?$",
                              full.names = TRUE))
if (length(fwd_files) == 0) {
  stop(sprintf("No *_1.fastq files found in %s", opt$fastq_dir))
}
cat(sprintf("Found %d forward read files\n", length(fwd_files)))

if (opt$mode == "paired") {
  rev_files <- sort(list.files(opt$fastq_dir, pattern = "_2\\.fastq(\\.gz)?$",
                                full.names = TRUE))
  if (length(rev_files) != length(fwd_files)) {
    stop(sprintf(
      "Paired mode: found %d fwd files but %d rev files in %s",
      length(fwd_files), length(rev_files), opt$fastq_dir
    ))
  }
  cat(sprintf("Found %d reverse read files\n", length(rev_files)))
}

sample_names <- sub("_1\\.fastq(\\.gz)?$", "", basename(fwd_files))

cat("Generating quality profiles...\n")

n_plot <- min(20, length(fwd_files))
plot_idx <- round(seq(1, length(fwd_files), length.out = n_plot))

pdf(file.path(opt$out_dir, "quality_fwd.pdf"), width = 10, height = 6)
tryCatch(
  print(plotQualityProfile(fwd_files[plot_idx])),
  error = function(e) message("Quality plot (fwd) failed: ", e$message)
)
dev.off()
cat(sprintf("Quality profile written: %s/quality_fwd.pdf\n", opt$out_dir))

if (opt$mode == "paired") {
  pdf(file.path(opt$out_dir, "quality_rev.pdf"), width = 10, height = 6)
  tryCatch(
    print(plotQualityProfile(rev_files[plot_idx])),
    error = function(e) message("Quality plot (rev) failed: ", e$message)
  )
  dev.off()
  cat(sprintf("Quality profile written: %s/quality_rev.pdf\n", opt$out_dir))
}

if (isTRUE(opt$inspect_only)) {
  cat(sprintf(
    "\n--inspect_only TRUE: stopping after quality plots.\n"
  ))
  cat(sprintf(
    "  Review quality_fwd.pdf (and quality_rev.pdf for paired) to choose\n"
  ))
  cat(sprintf(
    "  truncLen values. Forward reads typically retain quality to ~240-260 bp;\n"
  ))
  cat(sprintf(
    "  reverse reads typically degrade faster (truncate at ~160-200 bp).\n"
  ))
  cat(sprintf(
    "  Then re-run without --inspect_only TRUE.\n\n"
  ))
  quit(status = 0)
}

cat(sprintf(
  "Filtering & trimming: truncFwd=%d, truncRev=%d, maxEE=c(2,2)...\n",
  opt$trunc_fwd, opt$trunc_rev
))

filt_fwd <- file.path(filt_dir, paste0(sample_names, "_fwd_filt.fastq.gz"))
names(filt_fwd) <- sample_names

if (opt$mode == "paired") {
  filt_rev <- file.path(filt_dir, paste0(sample_names, "_rev_filt.fastq.gz"))
  names(filt_rev) <- sample_names

  # maxEE=2,2 is the dada2 default starting point for MiSeq
  out_filter <- filterAndTrim(
    fwd = fwd_files,
    filt = filt_fwd,
    rev = rev_files,
    filt.rev = filt_rev,
    truncLen = c(opt$trunc_fwd, opt$trunc_rev),
    maxN = 0,
    maxEE = c(2, 2),
    truncQ = 2,
    rm.phix = TRUE,
    compress = TRUE,
    multithread = opt$threads
  )
} else {
  out_filter <- filterAndTrim(
    fwd = fwd_files,
    filt = filt_fwd,
    truncLen = opt$trunc_fwd,
    maxN = 0,
    maxEE = 2,
    truncQ = 2,
    rm.phix = TRUE,
    compress = TRUE,
    multithread = opt$threads
  )
}

pct_passed <- round(100 * sum(out_filter[, 2]) / sum(out_filter[, 1]), 1)
cat(sprintf("Filtering: %d/%d reads passed (%.1f%%)\n",
            sum(out_filter[, 2]), sum(out_filter[, 1]), pct_passed))

if (pct_passed < 50) {
  warning(sprintf(
    "Only %.1f%% of reads passed filtering for %s. Consider relaxing truncLen or maxEE. Inspect quality_fwd.pdf.",
    pct_passed, opt$cohort
  ))
}

filt_fwd_exist <- filt_fwd[file.exists(filt_fwd)]
if (opt$mode == "paired") {
  filt_rev_exist <- filt_rev[names(filt_fwd_exist)]
}
n_dropped <- length(filt_fwd) - length(filt_fwd_exist)
if (n_dropped > 0) {
  cat(sprintf("Dropped %d samples with zero reads after filtering.\n",
              n_dropped))
}
sample_names_filt <- names(filt_fwd_exist)

cat("Learning error rates (forward)...\n")
err_fwd <- learnErrors(filt_fwd_exist, nbases = 1e8,
                        multithread = opt$threads, randomize = TRUE)

if (opt$mode == "paired") {
  cat("Learning error rates (reverse)...\n")
  err_rev <- learnErrors(filt_rev_exist, nbases = 1e8,
                          multithread = opt$threads, randomize = TRUE)
}

pdf(file.path(opt$out_dir, "error_model_fwd.pdf"), width = 10, height = 8)
tryCatch(print(plotErrors(err_fwd, nominalQ = TRUE)),
         error = function(e) message("Error plot failed: ", e$message))
dev.off()

if (opt$mode == "paired") {
  pdf(file.path(opt$out_dir, "error_model_rev.pdf"), width = 10, height = 8)
  tryCatch(print(plotErrors(err_rev, nominalQ = TRUE)),
           error = function(e) message("Error plot failed: ", e$message))
  dev.off()
}

cat("Running sample inference (pool='pseudo')...\n")
dada_fwd <- dada(filt_fwd_exist, err = err_fwd, pool = "pseudo",
                  multithread = opt$threads)

if (opt$mode == "paired") {
  dada_rev <- dada(filt_rev_exist, err = err_rev, pool = "pseudo",
                    multithread = opt$threads)
}

if (opt$mode == "paired") {
  cat("Merging paired reads...\n")
  mergers <- mergePairs(dada_fwd, filt_fwd_exist,
                         dada_rev, filt_rev_exist,
                         verbose = FALSE)
  cat(sprintf("Merging complete for %d samples.\n", length(mergers)))
}

cat("Building sequence table...\n")
if (opt$mode == "paired") {
  seqtab <- makeSequenceTable(mergers)
} else {
  seqtab <- makeSequenceTable(dada_fwd)
}

seq_lengths <- nchar(getSequences(seqtab))
cat(sprintf("Amplicon length distribution: min=%d, median=%d, max=%d\n",
            min(seq_lengths), as.integer(median(seq_lengths)), max(seq_lengths)))

cat("Removing chimeras...\n")
seqtab_nochim <- removeBimeraDenovo(seqtab, method = "consensus",
                                     multithread = opt$threads, verbose = FALSE)
pct_nochim <- round(100 * sum(seqtab_nochim) / sum(seqtab), 1)
cat(sprintf("After chimera removal: %.1f%% of reads retained.\n",
            pct_nochim))

if (pct_nochim < 70) {
  warning(sprintf(
    "Only %.1f%% of reads survived chimera removal for %s. Verify primer trimming was applied and that truncLen leaves sufficient overlap for mergePairs.",
    pct_nochim, opt$cohort
  ))
}

get_n <- function(x) sum(getUniques(x))

if (opt$mode == "paired") {
  stats <- data.frame(
    sample     = sample_names_filt,
    input      = out_filter[rownames(out_filter) %in%
                              paste0(sample_names_filt, "_1.fastq.gz") |
                              rownames(out_filter) %in%
                              paste0(sample_names_filt, "_1.fastq"), 1],
    filtered   = out_filter[rownames(out_filter) %in%
                              paste0(sample_names_filt, "_1.fastq.gz") |
                              rownames(out_filter) %in%
                              paste0(sample_names_filt, "_1.fastq"), 2],
    denoised_f = sapply(dada_fwd, get_n),
    denoised_r = sapply(dada_rev, get_n),
    merged     = sapply(mergers, function(x) sum(x$accept)),
    nonchim    = rowSums(seqtab_nochim)
  )
} else {
  stats <- data.frame(
    sample     = sample_names_filt,
    input      = out_filter[, 1],
    filtered   = out_filter[, 2],
    denoised_f = sapply(dada_fwd, get_n),
    nonchim    = rowSums(seqtab_nochim)
  )
}

write.csv(stats, file.path(opt$out_dir, "processing_stats.csv"),
          row.names = FALSE)
cat(sprintf("Processing stats written: %s/processing_stats.csv\n",
            opt$out_dir))

cat(sprintf("Assigning taxonomy with SILVA 138.1 (%s)...\n",
            opt$silva_path))

if (!file.exists(opt$silva_path)) {
  stop(sprintf(
    "SILVA training set not found at: %s\nDownload from Zenodo: https://zenodo.org/record/4587955 (silva_nr99_v138.1_train_set.fa.gz for assignTaxonomy, silva_species_assignment_v138.1.fa.gz for addSpecies)",
    opt$silva_path
  ))
}

taxa <- assignTaxonomy(seqtab_nochim, opt$silva_path,
                        minBoot = 80,
                        multithread = opt$threads,
                        verbose = FALSE)

silva_dir    <- dirname(opt$silva_path)
species_path <- file.path(silva_dir, "silva_species_assignment_v138.1.fa.gz")
if (file.exists(species_path)) {
  taxa <- addSpecies(taxa, species_path, verbose = FALSE)
  cat("Species-level assignments added.\n")
} else {
  cat(sprintf(
    "Species database not found at %s; skipping addSpecies.\n",
    species_path
  ))
}

asv_table <- t(seqtab_nochim)
write.csv(as.data.frame(asv_table),
          file.path(opt$out_dir, "asv_table.csv"))

write.csv(as.data.frame(taxa),
          file.path(opt$out_dir, "taxonomy_table.csv"))

genus_vec <- taxa[, "Genus"]
genus_vec[is.na(genus_vec)] <- "Unknown"

genus_table <- t(rowsum(t(seqtab_nochim), genus_vec))
write.csv(as.data.frame(genus_table),
          file.path(opt$out_dir, "genus_table.csv"))

saveRDS(
  list(seqtab_nochim = seqtab_nochim, taxa = taxa, stats = stats,
       err_fwd = err_fwd, cohort = opt$cohort, mode = opt$mode),
  file.path(opt$out_dir, paste0(opt$cohort, "_dada2.rds"))
)

n_asvs  <- ncol(seqtab_nochim)
n_genera <- ncol(genus_table)
cat(sprintf(
  "\nDone: %d ASVs, %d genera, %d samples for cohort '%s'.\n",
  n_asvs, n_genera, nrow(seqtab_nochim), opt$cohort
))
cat(sprintf("Outputs in: %s\n\n", opt$out_dir))
