#!/usr/bin/env Rscript

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

# run once on a fresh machine
BiocManager::install(version = "3.22", ask = FALSE, update = FALSE)

cran_pkgs <- c("optparse", "vegan")
for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    cat(sprintf("Installing %s from CRAN...\n", p))
    install.packages(p, repos = "https://cloud.r-project.org")
  } else {
    cat(sprintf("%-14s already installed (%s)\n", p, packageVersion(p)))
  }
}

bioc_pkgs <- c("dada2", "sva", "MMUPHin", "phyloseq", "zCompositions")
for (p in bioc_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    cat(sprintf("Installing %s from Bioconductor...\n", p))
    BiocManager::install(p, ask = FALSE, update = FALSE)
  } else {
    cat(sprintf("%-14s already installed (%s)\n", p, packageVersion(p)))
  }
}

cat("\nPackage status:\n")
all_pkgs <- c(cran_pkgs, bioc_pkgs)
for (p in all_pkgs) {
  ok <- requireNamespace(p, quietly = TRUE)
  cat(sprintf("%-14s %s\n", p,
              if (ok) paste0("OK (", packageVersion(p), ")") else "FAILED"))
}
