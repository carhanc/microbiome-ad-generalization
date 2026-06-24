# Supplementary Materials

**Title:** Cross-Cohort Generalization Failure in Gut Microbiome-Based Alzheimer's Disease Classifiers: Evidence from Five Independent Cohorts with Batch Correction Analysis

---

## Supplementary Tables

### Table S1. Full Pairwise Cross-Cohort AUC Matrix

Pairwise AUC-ROC for single-cohort train → single-cohort test transfer experiments. All 4×4 combinations excluding diagonal (within-cohort). Values are AUC-ROC with 95% bootstrap confidence intervals (1,000 resamples) on the test cohort predictions.

**Table S1A — Logistic Regression**

| Train \ Test | Zhuang 2018 | Ling 2020 | Zhu 2022 | Kazakhstan |
|---|---|---|---|---|
| Zhuang 2018 | — | 0.583 [0.50–0.67] | 0.696 [0.54–0.82] | 0.566 [0.44–0.69] |
| Ling 2020 | 0.538 [0.41–0.66] | — | 0.879 [0.77–0.96] | 0.593 [0.47–0.72] |
| Zhu 2022 | 0.579 [0.46–0.70] | 0.760 [0.68–0.83] | — | 0.547 [0.43–0.67] |
| Kazakhstan | 0.539 [0.42–0.67] | 0.739 [0.66–0.81] | 0.652 [0.50–0.78] | — |

**Table S1B — LightGBM**

| Train \ Test | Zhuang 2018 | Ling 2020 | Zhu 2022 | Kazakhstan |
|---|---|---|---|---|
| Zhuang 2018 | — | 0.598 [0.51–0.68] | 0.604 [0.46–0.74] | 0.493 [0.37–0.61] |
| Ling 2020 | 0.617 [0.49–0.74] | — | 0.861 [0.77–0.94] | 0.602 [0.47–0.72] |
| Zhu 2022 | 0.528 [0.41–0.65] | 0.591 [0.50–0.68] | — | 0.405 [0.28–0.53] |
| Kazakhstan | 0.522 [0.40–0.64] | 0.580 [0.49–0.67] | 0.603 [0.45–0.75] | — |

Notable observations: The highest pairwise logistic regression AUC is Ling 2020 → Zhu 2022 (AUC=0.88), reflecting transfer between two large Chinese fecal cohorts. Kazakhstan → Ling 2020 also transfers reasonably (logistic regression AUC=0.74). Most pairwise AUCs cluster between 0.54 and 0.70, consistent with the overall pattern of limited cross-cohort transferability documented in Section 3.3.

---

### Table S2. Per-Cohort DADA2 Processing Statistics

Statistics from DADA2 amplicon sequence variant (ASV) inference pipeline. Paired-end cohorts underwent quality filtering, error learning, denoising, merging, and chimera removal; Kazakhstan (single-end) omitted the merging step.

| Cohort | Country | Mode | N Samples | Reads Input | Reads Non-chimeric | % Retention | 16S Region | Platform |
|---|---|---|---|---|---|---|---|---|
| Zhuang 2018 | China | Paired-end | 86 | 3,260,213 | 2,551,870 | 78.3% | V3–V4 | MiSeq 2×300 |
| Ling 2020 | China | Paired-end | 171 | 8,556,209 | 4,100,210 | 47.9% | V3–V4 | MiSeq 2×300 |
| Zhu 2022 | China | Paired-end | 180 | 7,275,079 | 5,485,811 | 75.4% | V3–V4 | MiSeq 2×300 |
| Kim 2022 | South Korea | Paired-end | 78 | 13,514,374 | 6,117,631 | 45.3% | V3–V4 | MiSeq 2×300 |
| Kaiyrlykyzy 2022 | Kazakhstan | Single-end | 84 | 9,467,007 | 6,531,362 | 69.0% | V3–V4 | NovaSeq |
| **Total** | | | **599** | **42,072,882** | **24,786,884** | **58.9%** | | |

Notes: (1) Reads input = reads passing initial quality cutoff prior to DADA2 filtering. (2) Retention rates for Ling 2020 (47.9%) and Kim 2022 (45.3%) are lower than typical for MiSeq 2×300 V3–V4 data; inspection of quality profiles showed elevated error rates in these cohorts' reverse reads at positions >150–180 bp, which may reflect sequencing run-specific quality variation. (3) Kazakhstan single-end processing used truncLen=250 after quality profile inspection; see Supplementary Note S2. (4) ASV count per cohort and per-sample read depth statistics not shown here; final post-harmonization genus matrix contains 396 genera across 509 fecal samples (90 blood microbiome samples from Zhu 2022 excluded prior to analysis).

---

### Table S3. Full PERMANOVA Results

All PERMANOVA models used 9,999 permutations with Type III (marginal) sums of squares (vegan adonis2, `by="margin"`). Aitchison distance = Euclidean distance on CLR-transformed values. Bray-Curtis computed on relative abundances.

**Model 1: Cohort only (all 5 cohorts, Aitchison distance)**

| Term | Df | Sum of Squares | R² | F | p-value |
|---|---|---|---|---|---|
| Cohort | 4 | 67,760.8 | 0.193 | 30.13 | <0.0001 |
| Residual | 504 | 283,347.1 | 0.807 | — | — |

**Model 1B: Cohort only (all 5 cohorts, Bray-Curtis dissimilarity)**

| Term | Df | Sum of Squares | R² | F | p-value |
|---|---|---|---|---|---|
| Cohort | 4 | 20.80 | 0.136 | 19.80 | <0.0001 |
| Residual | 504 | 132.36 | 0.864 | — | — |

**Model 2: Cohort + Diagnosis (4 labeled cohorts, Aitchison, marginal R²)**

| Term | Df | Sum of Squares | R² (marginal) | F | p-value |
|---|---|---|---|---|---|
| Cohort | 3 | 45,476.8 | 0.172 | 28.03 | <0.0001 |
| Diagnosis | 1 | 3,689.7 | 0.014 | 6.82 | <0.0001 |

Note: Marginal R² values for cohort and diagnosis do not sum to total R² because marginal effects are estimated independently (Type III SS), adjusting for the other variable. The sum of marginal R² values underestimates total explained variance.

**Model 3: Within-cohort Diagnosis effects (Aitchison distance)**

| Cohort | Term | Df | R² | F | p-value |
|---|---|---|---|---|---|
| Zhuang 2018 | Diagnosis | 1 | 0.017 | 1.47 | 0.066 |
| Ling 2020 | Diagnosis | 1 | 0.053 | 9.53 | <0.0001 |
| Zhu 2022 | Diagnosis | 1 | 0.094 | 6.01 | <0.0001 |
| Kazakhstan | Diagnosis | 1 | 0.022 | 1.83 | 0.069 |

Zhu 2022 analysis restricted to fecal samples only (n=60 binary: 30 AD + 30 CN). Blood microbiome samples excluded prior to all analyses. Zhuang 2018 (p=0.066) and Kazakhstan (p=0.069) did not reach α=0.05 significance for within-cohort diagnosis effect.

**Beta-Dispersion Test (Aitchison distance, all 5 cohorts)**

| Metric | Value |
|---|---|
| F-statistic | 13.93 |
| Numerator df | 4 |
| Denominator df | 594 |
| p-value (permutation, 9,999 perms) | <0.0001 |
| Interpretation | Within-cohort dispersions differ significantly; PERMANOVA R² reflects both centroid shift and spread heterogeneity |

---

## Supplementary Figures

### Figure S1. Zhu 2022: Fecal-Only Within-Cohort AUC

**Caption:** Within-cohort AUC-ROC for Zhu 2022 fecal samples only (n=60 binary: 30 AD + 30 CN), estimated using standard StratifiedKFold 10-fold nested cross-validation. Blood microbiome (B_*) samples were excluded from all analyses; each participant appears exactly once in the fecal-only dataset, making standard StratifiedKFold appropriate. Logistic regression achieves AUC=0.998 [0.99–1.00] and LightGBM achieves AUC=0.979 [0.95–1.00]. The near-perfect within-cohort AUC should be interpreted cautiously given the small sample size (30 per class); see Section 3.2 and Section 4.3 for discussion. Error bars = 95% bootstrap CI (1,000 resamples).

*Source data:* results/model_outputs/within_cohort_cv/shanghai2022_*.csv

---

### Figure S2. LightGBM Top-15 SHAP Taxa Per Cohort

**Caption:** SHAP feature importance dot plot for LightGBM classifiers, parallel to Figure 6A (logistic regression). Dot size proportional to mean |SHAP| value (feature importance); red = AD-associated (higher mean SHAP in AD samples); blue = CN-associated. Yellow background marks directional flip taxa appearing in top-20 of ≥2 cohorts with opposite sign (★ prefix). Top taxon per cohort: Oscillibacter (Zhuang 2018, CN-associated, |SHAP|=1.025), Akkermansia (Ling 2020, AD-associated, |SHAP|=0.828), Streptococcus (Zhu 2022, AD-associated, |SHAP|=1.287), Castellaniella (Kazakhstan, CN-associated, |SHAP|=0.830). Mean pairwise Jaccard similarity at top-20: 0.026 — lower than logistic regression (0.101), consistent with tree models finding more cohort-specific features. LightGBM directional flips (top-20): Incertae Sedis, Lactobacillus, Oscillibacter.

*Source figure:* `results/figures/shap_top15_per_cohort_lgbm.png`

---

### Figure S3. LOCO SHAP Analysis — Feature Importance in Held-Out Test Cohorts

**Caption:** SHAP values computed for LOCO-trained classifiers (trained on N-1 cohorts) evaluated on each held-out test cohort. Each panel shows the top-15 genera by mean |SHAP| on the held-out test data, with direction (AD or CN associated) as computed on the test samples. (A) Logistic regression. (B) LightGBM.

Top-1 LOCO SHAP taxon per held-out cohort:
- Held-out Zhuang 2018: Akkermansia (LogReg, |SHAP|=0.612), Akkermansia (LGBM)
- Held-out Ling 2020: Eisenbergiella (LogReg), Streptococcus (LGBM)
- Held-out Zhu 2022: Faecalibacterium (LogReg), Akkermansia (LGBM)
- Held-out Kazakhstan: Cutibacterium (LogReg), [Eubacterium] xylanophilum group (LGBM)

No single genus ranks first in more than one held-out condition for logistic regression. The four LOCO top taxa are all different, confirming that cross-cohort models activate different features depending on the test population — there is no universal transferable taxon. This provides the mechanistic explanation for LOCO generalization failure: the model relies on genera that are uninformative or directionally inappropriate for the test population.

*Source figure:* `results/figures/shap_loco_direction_logreg.png`, `results/figures/shap_loco_direction_lgbm.png`

---

### Figure S4. PCoA Before and After ComBat-seq Correction

**Caption:** Principal coordinates analysis (PCoA) on Aitchison distance (CLR-Euclidean) colored by cohort (left column) and diagnosis (right column), before (top row) and after ComBat-seq correction (bottom row). Samples restricted to four labeled cohorts (n=401). After ComBat-seq correction, cohort-level separation on PC1 and PC2 is visually reduced, consistent with the correction algorithm aligning cohort means. However, diagnosis does not become more separable after correction — diagnoses remain interleaved within and across cohorts — providing a visual complement to the LOCO AUC results showing no improvement in cross-cohort discrimination.

*Source figure:* `results/figures/pcoa_combatseq_before_after.png`

---

### Figure S5. PCoA Before and After MMUPHin Correction

**Caption:** Principal coordinates analysis (PCoA) on Aitchison distance colored by cohort (left) and diagnosis (right), before (top) and after MMUPHin correction (bottom). Four labeled cohorts (n=401). MMUPHin correction produces partial cohort alignment — cohorts are more overlapping in PCoA than before correction — but the diagnosis axis remains unresolved after correction, consistent with MMUPHin LOCO AUC failing to improve upon uncorrected performance. The continued interleaving of AD and CN samples in PCoA after correction visualizes the ~12-fold cohort-to-diagnosis variance ratio: the disease signal occupies a small fraction of compositional space that cannot be selectively preserved when removing the dominant cohort signal.

*Source figure:* `results/figures/pcoa_mmuphin_before_after.png`

---

## Supplementary Notes

### Note S1. Kim/KBASE 2022 (PRJEB50447) — IRB Restriction on Per-Sample Diagnosis Labels

**Summary:** The Kim/KBASE 2022 cohort (78 participants; South Korea; 18 amyloid-PET positive preclinical AD, 60 cognitively normal) was excluded from all supervised classifier analyses because per-sample amyloid-PET status and group labels are not publicly deposited in any data archive.

**Documentation of exhaustive metadata search conducted 2026-06-20:**

The following metadata sources were queried for per-sample diagnosis labels:

1. **ENA Portal filereport** (https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB50447&result=read_run&fields=all): All sample fields including `sample_alias`, `sample_title`, `scientific_name`, `sample_description` — no diagnostic labels present. Fields contain sequencing run identifiers only.

2. **ENA BioSample XML**: Retrieved individual BioSample records via ENA Browser. Attributes present: collection date, geographic location, library strategy/layout — no diagnosis, amyloid status, or clinical phenotype fields.

3. **NCBI SRA runinfo** (via `esearch -db sra -query "PRJEB50447" | efetch -format runinfo`): Confirmed SRA runinfo does not contain phenotype data beyond SRA-standard fields. Disease column is empty for all runs.

4. **Sample alias and title fields**: `sample_alias` values are coded identifiers (e.g., "KBASE-CN-001") without diagnostic information inferrable from the coding scheme.

**Published paper finding:** The corresponding publication (Kim et al., 2022; PRJEB50447) reports group-level counts (18 amyloid-positive preclinical AD; 60 cognitively normal) and notes that per-sample data are withheld under institutional IRB protocol. The amyloid-PET classifications are considered clinical data and are protected. The published paper does not provide a supplementary file with per-sample labels.

**Consequence for analysis:** kbase2022 was included only in PERMANOVA variance decomposition (Phase 4), where cohort identity — not per-sample diagnosis — is the grouping variable of interest. All supervised classifier training, evaluation, batch correction, and SHAP analyses were restricted to the four labeled cohorts (Zhuang 2018, Ling 2020, Zhu 2022, Kazakhstan).

**Additional caveat:** Even if per-sample labels were available, direct comparison with the other four cohorts would require caution. Kim/KBASE 2022 uses amyloid-PET positivity as the case definition (preclinical AD — cognitively normal by neuropsychological testing but amyloid-positive), while the other four cohorts use clinical AD diagnosis (cognitive impairment criterion). These phenotype definitions are not equivalent.

**Recommendation for future analyses:** Researchers wishing to include Kim/KBASE 2022 in supervised analyses should contact the study team (KBASE Consortium / corresponding author of Kim et al., 2022) to request a data use agreement for per-sample amyloid-PET status.

---

### Note S2. Kazakhstan Single-End Processing Rationale

**Summary:** The Kazakhstan cohort (PRJNA811324; Kaiyrlykyzy et al., 2022) was confirmed as single-end sequencing and processed using the DADA2 single-end workflow. This document explains the processing decisions made.

**Confirmation of single-end format:** After downloading all 84 samples with `fasterq-dump`, output files were `SRR*_1.fastq` only (no `SRR*_2.fastq` files). This was independently confirmed by inspecting SRA metadata (LibraryLayout=SINGLE). The original publication describes NovaSeq sequencing of V3–V4 amplicons.

**Truncation length decision (250 bp):** Quality profile inspection of Kazakhstan reads showed a characteristic NovaSeq quality pattern: high-quality bases through approximately position 250 bp, with rapid quality decline thereafter. MiSeq 2×300 bp paired-end reads for V3–V4 have effective merged read lengths of approximately 400–450 bp; single-end NovaSeq reads truncated at 250 bp cover a shorter portion of the amplicon. The truncLen parameter was set to 250 bp after visual inspection of per-position quality score profiles across 10 representative samples.

**Implications for taxonomic resolution:** Single-end reads provide lower taxonomic resolution than paired-end merged reads for the same amplicon region, because the additional bases from the reverse read improve amplicon coverage and reduce ambiguous taxonomic assignments. This means genus-level assignments for Kazakhstan may be slightly less accurate than for the four paired-end cohorts, and species-level analysis would be more affected than genus-level. We acknowledge this as a limitation of the mixed-library-design dataset.

**DADA2 workflow difference:** The single-end workflow did not include `mergePairs()`. The full single-end DADA2 pipeline applied: `filterAndTrim()` → `learnErrors()` → `dada()` → `makeSequenceTable()` → `removeBimeraDenovo()`. Error learning was performed independently on Kazakhstan samples; sharing error models across cohorts would be inappropriate because error profiles differ by sequencer and run.

**Interaction with batch correction:** The single-end/paired-end asymmetry is an additional source of between-cohort technical heterogeneity that batch correction methods cannot fully address, because the fundamental difference in read length and sequencing chemistry affects all genera' abundance estimates rather than a subset. This is a study limitation that applies to any analysis including Kazakhstan and should be acknowledged in the Methods section of any manuscript using this dataset.

---

*End of Supplementary Materials*
