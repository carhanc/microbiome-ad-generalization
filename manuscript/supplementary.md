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

All PERMANOVA models used 9,999 permutations with Type III (marginal) sums of squares (vegan adonis2, `by="margin"`). Aitchison distance = Euclidean distance on CLR-transformed values. Bray-Curtis computed on relative abundances. With 9,999 permutations, the minimum attainable permutation p-value is 1/(9,999+1) = 0.0001; p-values reported below as 0.0001 indicate no permuted statistic met or exceeded the observed value (i.e., the true p-value is at or below this resolution limit, not literally less than it).

**Model 1: Cohort only (all 5 cohorts, Aitchison distance)**

| Term | Df | Sum of Squares | R² | F | p-value |
|---|---|---|---|---|---|
| Cohort | 4 | 67,760.8 | 0.193 | 30.13 | 0.0001 |
| Residual | 504 | 283,347.1 | 0.807 | — | — |

**Model 1B: Cohort only (all 5 cohorts, Bray-Curtis dissimilarity)**

| Term | Df | Sum of Squares | R² | F | p-value |
|---|---|---|---|---|---|
| Cohort | 4 | 20.80 | 0.136 | 19.80 | 0.0001 |
| Residual | 504 | 132.36 | 0.864 | — | — |

**Model 2: Cohort + Diagnosis (4 labeled cohorts, Aitchison, marginal R²)**

| Term | Df | Sum of Squares | R² (marginal) | F | p-value |
|---|---|---|---|---|---|
| Cohort | 3 | 45,476.8 | 0.172 | 28.03 | 0.0001 |
| Diagnosis | 1 | 3,689.7 | 0.014 | 6.82 | 0.0001 |

Note: Marginal R² values for cohort and diagnosis do not sum to total R² because marginal effects are estimated independently (Type III SS), adjusting for the other variable. The sum of marginal R² values underestimates total explained variance.

**Model 3: Within-cohort Diagnosis effects (Aitchison distance)**

| Cohort | Term | Df | R² | F | p-value |
|---|---|---|---|---|---|
| Zhuang 2018 | Diagnosis | 1 | 0.017 | 1.47 | 0.066 |
| Ling 2020 | Diagnosis | 1 | 0.053 | 9.53 | 0.0001 |
| Zhu 2022 | Diagnosis | 1 | 0.094 | 6.01 | 0.0001 |
| Kazakhstan | Diagnosis | 1 | 0.022 | 1.83 | 0.069 |

Zhu 2022 analysis restricted to fecal samples only (n=60 binary: 30 AD + 30 CN). Blood microbiome samples excluded prior to all analyses. Zhuang 2018 (p=0.066) and Kazakhstan (p=0.069) did not reach α=0.05 significance for within-cohort diagnosis effect.

**Beta-Dispersion Test (Aitchison distance, all 5 cohorts)**

| Metric | Value |
|---|---|
| F-statistic | 13.93 |
| Numerator df | 4 |
| Denominator df | 504 |
| p-value (permutation, 9,999 perms) | 0.0001 |
| Interpretation | Within-cohort dispersions differ significantly; PERMANOVA R² reflects both centroid shift and spread heterogeneity |

---

### Table S4. Training-Only Feature-Selection Sensitivity (Round 2)

Genus retention (≥20% prevalence) re-derived independently within every split, using only that split's training-side samples (Section 2.5.1), compared against the main analysis's global (all-five-cohort) feature selection.

**LOCO**

| Test Cohort | Model | Original (Global-Feature) LOCO AUC | Training-Only LOCO AUC | Δ | N Retained Genera (Training-Only) |
|---|---|---|---|---|---|
| Zhuang 2018 | LogReg | 0.5041 | 0.5057 | +0.0016 | 370 |
| Zhuang 2018 | LGBM | 0.5646 | 0.5960 | +0.0314 | 370 |
| Ling 2020 | LogReg | 0.6775 | 0.6794 | +0.0019 | 370 |
| Ling 2020 | LGBM | 0.6735 | 0.6534 | −0.0201 | 370 |
| Zhu 2022 | LogReg | 0.8133 | 0.8122 | −0.0011 | 373 |
| Zhu 2022 | LGBM | 0.7578 | 0.7878 | +0.0300 | 373 |
| Kazakhstan | LogReg | 0.5644 | 0.6092 | +0.0448 | 150 |
| Kazakhstan | LGBM | 0.5905 | 0.5247 | −0.0658 | 150 |

**Within-cohort nested CV**

| Cohort | Model | Original (Global-Feature) AUC | Training-Only AUC | Δ |
|---|---|---|---|---|
| Zhuang 2018 | LogReg | 0.6333 | 0.6647 | +0.0314 |
| Zhuang 2018 | LGBM | 0.6349 | 0.6733 | +0.0384 |
| Ling 2020 | LogReg | 0.8648 | 0.8896 | +0.0248 |
| Ling 2020 | LGBM | 0.8618 | 0.8415 | −0.0203 |
| Zhu 2022 | LogReg | 0.9978 | 0.9956 | −0.0022 |
| Zhu 2022 | LGBM | 0.9789 | 0.9700 | −0.0089 |
| Kazakhstan | LogReg | 0.7657 | 0.7368 | −0.0289 |
| Kazakhstan | LGBM | 0.6943 | 0.6109 | −0.0834 |

The retained genus count for LOCO training-only selection is lower when Kazakhstan is excluded from feature selection (150 genera) than when it is included in the training pool (370–373 genera), reflecting that Kazakhstan alone passes the 20%-prevalence threshold for many genera the Chinese/Korean-pipeline cohorts do not. No cohort–model combination changes qualitative classification (near-chance, intermediate, robust) between the global and training-only designs.

*Source data:* `results/tables/within_cohort_auc_training_only.csv`, `results/tables/loco_auc_training_only.csv`, `results/tables/pairwise_auc_training_only.csv`, `results/tables/training_only_feature_sensitivity_summary.csv`, `results/tables/training_only_feature_counts.csv`, `results/tables/training_only_feature_lists.csv` (full per-split retained-genus lists).

---

### Table S5. Within-Cohort AUC on Batch-Corrected Data: Label-Blind vs. Label-Informed (Round 2)

Within-cohort nested-CV AUC (Section 2.4 procedure) computed on batch-corrected features, under both correction designs (Section 2.7). The label-blind comparison is unconfounded by label use during correction; the label-informed comparison is confounded (Section 3.5) and is shown only for completeness/continuity with the original analysis.

| Cohort | Model | Uncorrected | ComBat-seq (label-blind) | MMUPHin (label-blind) | ComBat-seq (label-informed, exploratory) | MMUPHin (label-informed, exploratory) |
|---|---|---|---|---|---|---|
| Zhuang 2018 | LogReg | 0.6333 | 0.6604 | 0.6171 | 0.6647 | 0.6393 |
| Zhuang 2018 | LGBM | 0.6349 | 0.6479 | 0.6247 | 0.6544 | 0.5765 |
| Ling 2020 | LogReg | 0.8648 | 0.8701 | 0.8690 | 0.8754 | 0.8863 |
| Ling 2020 | LGBM | 0.8618 | 0.8255 | 0.8689 | 0.8659 | 0.8844 |
| Zhu 2022 | LogReg | 0.9978 | 0.9878 | 0.9900 | 0.9956 | 0.9933 |
| Zhu 2022 | LGBM | 0.9789 | 0.9533 | 0.9600 | 0.9622 | 0.9650 |
| Kazakhstan | LogReg | 0.7657 | 0.6926 | 0.7663 | 0.7839 | 0.7879 |
| Kazakhstan | LGBM | 0.6943 | 0.6829 | 0.5621 | 0.7351 | 0.6858 |

*Source data:* `results/tables/within_cohort_auc_corrected_labelblind.csv`, `results/tables/within_cohort_auc_corrected.csv`, `results/tables/within_cohort_auc.csv`.

---

### Table S6. Logistic-Regression Coefficient Sign Stability: Screened Candidates and the *Akkermansia* Negative Control (Round 2)

Per (cohort, genus), counts of positive/negative/zero fitted coefficients across the 10 outer OOF folds (Section 2.8). `sign_stability` uses the pre-specified ≥8/10-fold threshold. This is a descriptive screen, not a hypothesis test: the 10 outer folds share overlapping training data and are not independent replications, so no p-value is assigned to any row. Full table (1,584 cohort×genus rows) in `results/tables/logreg_coefficient_stability.csv`; below, *Romboutsia* (the cleanest population-level pattern among the eight genera meeting the screen; Section 3.6) and *Akkermansia* (a negative-control illustration: appears to flip under the old, invalid mean-SHAP statistic but does not meet the coefficient-based screen; Section 3.6, Section 4.4).

| Cohort | Taxon | n_positive/10 | n_negative/10 | Median Coefficient | Sign Stability |
|---|---|---|---|---|---|
| Zhuang 2018 | Romboutsia | 0 | 10 | −0.0107 | stable_negative |
| Ling 2020 | Romboutsia | 0 | 10 | −0.1166 | stable_negative |
| Zhu 2022 | Romboutsia | 0 | 10 | −0.0190 | stable_negative |
| Kazakhstan | Romboutsia | 9 | 1 | 0.0109 | stable_positive |
| Zhuang 2018 | Akkermansia | 10 | 0 | 0.0124 | stable_positive |
| Ling 2020 | Akkermansia | 10 | 0 | 0.1181 | stable_positive |
| Zhu 2022 | Akkermansia | 10 | 0 | 0.0709 | stable_positive |
| Kazakhstan | Akkermansia | 10 | 0 | 0.0682 | stable_positive |

*Romboutsia* is stable-negative in all three Chinese cohorts and stable-positive in Kazakhstan — the cleanest country-level split among the eight genera meeting the screen. *Akkermansia*'s coefficient is stable and positive in all four cohorts — it does not meet the screen at all, despite its pooled mean SHAP value appearing to reverse sign (Section 3.6). See the full CSV for the corresponding rows for each of the eight genera meeting the descriptive screen (*Agathobacter*, *Bifidobacterium*, *Coprococcus*, *Dorea*, *Lactobacillus*, *NK4A214 group*, *Romboutsia*, *Ruminococcus gnavus group*) and for `logreg_directional_flips_stable.csv`, which lists, per candidate taxon, which cohorts carried a stable positive vs. stable negative coefficient.

*Source data:* `results/tables/logreg_coefficient_stability.csv`, `results/tables/logreg_directional_flips_stable.csv`.

---

### Table S7. Matched-Null Jaccard Overlap Across Top-N Thresholds (Round 2)

Observed statistic = mean of the six unique pairwise Jaccard similarities among the four labeled cohorts' top-N SHAP taxa sets. Null = 100,000 replicates, each drawing four independent random top-N gene sets from the 396-genus pool and averaging the same six pairwise Jaccards among them (Section 2.8). Empirical one-sided p = (1 + #{null ≥ observed}) / (n_perm + 1).

| Model | Top-N | Observed Mean (6 pairs) | Null Mean | Null 95% CI | Empirical p |
|---|---|---|---|---|---|
| LogReg | 10 | 0.0850 | 0.0134 | [0.0000, 0.0361] | <0.0001 |
| LogReg | **20 (primary)** | **0.1350** | 0.0265 | [0.0085, 0.0491] | **<0.0001** |
| LogReg | 50 | 0.2296 | 0.0680 | [0.0491, 0.0895] | <0.0001 |
| LightGBM | 10 | 0.0175 | 0.0135 | [0.0000, 0.0361] | 0.4494 (not significant) |
| LightGBM | **20 (primary)** | **0.0584** | 0.0265 | [0.0085, 0.0491] | **0.0048** |
| LightGBM | 50 | 0.1430 | 0.0680 | [0.0492, 0.0894] | <0.0001 |

LightGBM's overlap is significant at N=20 and N=50 but not at N=10, indicating a real but small-magnitude and N-sensitive signal — smaller and less consistent than logistic regression's overlap at every threshold tested. This supersedes an earlier, incorrectly-calibrated null (Section 2.8) that had reported LightGBM's N=20 overlap as "within the random baseline" (p=0.068); that earlier null drew only two random sets per replicate rather than the four-set, six-pair statistic actually being tested.

*Source data:* `results/tables/jaccard_null_corrected.csv`. The superseded `results/tables/jaccard_null.csv` is retained on disk for audit-trail provenance only and is not used by this manuscript.

---

## Supplementary Figures

### Figure S1. Zhu 2022: Fecal-Only Within-Cohort AUC

**Caption:** Within-cohort AUC-ROC for Zhu 2022 fecal samples only (n=60 binary: 30 AD + 30 CN), estimated using standard StratifiedKFold 10-fold nested cross-validation. Blood microbiome (B_*) samples were excluded from all analyses; each participant appears exactly once in the fecal-only dataset, making standard StratifiedKFold appropriate. Logistic regression achieves AUC=0.998 [0.99–1.00] and LightGBM achieves AUC=0.979 [0.95–1.00]. The near-perfect within-cohort AUC should be interpreted cautiously given the small sample size (30 per class); see Section 3.2 and Section 4.3 for discussion. Error bars = 95% bootstrap CI (1,000 resamples).

*Source data:* results/model_outputs/within_cohort_cv/shanghai2022_*.csv

---

### Figure S2. Top-15 SHAP Taxa Per Cohort (LightGBM, Out-of-Fold)

**Caption (revised, Round 2):** Horizontal bar charts showing mean |SHAP| value (CLR units) for the top-15 genera by importance in each of the four labeled cohorts under out-of-fold LightGBM evaluation. Bars are a single neutral color and show feature-importance magnitude only — **no AD/CN direction is assigned**, because a tree ensemble has no single coefficient-like global direction and LightGBM's SHAP-feature relationships may be nonlinear or non-monotonic (Section 2.8). Top taxon per cohort (by |SHAP| magnitude only): Lachnoclostridium (Zhuang 2018, |SHAP|=0.560), Akkermansia (Ling 2020, |SHAP|=0.947), Bacteroides (Zhu 2022, |SHAP|=1.243), Castellaniella (Kazakhstan, |SHAP|=0.614). Mean pairwise Jaccard similarity at top-20: 0.058, significantly above a null distribution matched to this statistic's construction (p=0.0048 at N=20; not significant at N=10, p=0.45; Supplementary Table S7). **No LightGBM directional-flip claim is made in this paper**; an earlier draft's report of "12 LightGBM directional flips" (based on an invalid mean-signed-SHAP statistic) is superseded and should not be cited.

*Source figure:* `results/figures/supp_fig_s2.jpg` (built from `results/figures/working/shap_top15_per_cohort_lgbm.png`, regenerated Round 2 with neutral-color bars)

---

### Figure S3a/S3b. LOCO SHAP Direction — Removed in Round 2

**These two figures (originally: out-of-fold SHAP values for LOCO-trained classifiers, plotted with a signed x-position labeled "AD-associated"/"CN-associated") have been removed.** Reviewer 3, Round 2, concern #3 established that mean signed SHAP is not a valid direction statistic in general, and this is especially severe for LOCO SHAP specifically: the SHAP background (the three pooled training cohorts) and the evaluation set (the held-out cohort) are drawn from systematically different populations by construction, so the sign of a LOCO mean SHAP value is confounded with the raw cross-cohort compositional difference quantified by PERMANOVA (Section 3.4), independent of anything the model learned. No coefficient-based replacement direction is offered for LOCO SHAP either, because LOCO models are refit on pooled multi-cohort training data each time and a single stable-fold coefficient analysis analogous to Section 2.8's within-cohort procedure was not part of this design.

LOCO SHAP **feature importance** (mean |SHAP|, magnitude only, no direction) remains available and is unaffected by this issue: see `results/tables/shap_loco_importance.csv`. For reference, the previously reported top-importance genus per held-out cohort was: Zhuang 2018 — *Akkermansia* (logistic regression |SHAP|=0.470); Ling 2020 — *Subdoligranulum* (|SHAP|=0.269); Zhu 2022 — *Akkermansia* (|SHAP|=0.372); Kazakhstan — *Christensenellaceae R-7 group* (|SHAP|=0.368); no genus ranked first across all four held-out conditions for either model (Section 3.6).

---

### Figure S4. Sensitivity Analysis: LOCO AUC Under Cohort Exclusion

**Caption:** Sensitivity analysis: LOCO AUC under three cohort exclusion configurations. Each panel shows the 3-cohort LOCO experiment with one cohort excluded from the full analysis. Left: excluding Kazakhstan (Zhuang+Ling+Zhu trained, each held out in turn). Center: excluding Zhuang 2018. Right: excluding Ling 2020. Error bars = 95% bootstrap CI. Generalization degradation persisted across all three prespecified cohort-exclusion configurations.

*Source data:* `results/tables/sensitivity_loco.csv`
*Source figure:* `results/figures/supp_sensitivity_loco (supp_fig_s4).jpg`

---

### Figure S5. Batch Correction, Label-Informed Design (Exploratory Sensitivity, Round 2)

**Caption:** LOCO AUC-ROC under the label-informed transductive batch-correction design (ComBat-seq/MMUPHin fit using every sample's own true diagnosis label; Section 2.7), shown for (**A**) logistic regression and (**B**) LightGBM. This is the same comparison previously presented as the paper's primary batch-correction result; it is retained here as an **explicitly-labeled exploratory/oracle-style sensitivity only**, because a correction step that has already seen each sample's true label cannot estimate prospective external-validation performance, and because the resulting within-cohort AUC on corrected data cannot be used as evidence that genuine, correction-independent biological signal survived (Section 3.5). The main-text Figure 5 shows the label-blind design (no diagnosis information used during correction), which we now treat as the primary batch-correction result. Note the severe, cohort-specific degradation visible here (e.g., Kazakhstan logistic regression LOCO AUC=0.308) that is **not** reproduced under the label-blind design (Kazakhstan label-blind LOCO AUC=0.560, closely matching the uncorrected value of 0.564) — indicating this degradation was driven substantially by the label-informed design itself.

*Source data:* `results/tables/auc_comparison_table.csv`, `results/tables/loco_auc_corrected.csv`
*Source figure:* `results/figures/supp_fig_s5.jpg`

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

**Cohort search strategy and full inclusion/exclusion record:** The cohort identification search (Section 2.1) queried NCBI SRA and EBI-ENA for studies deposited before January 2025, using the following search string in the SRA search interface: `("Alzheimer's disease" OR "Alzheimer disease" OR "mild cognitive impairment" OR "MCI" OR "dementia") AND ("gut microbiome" OR "gut microbiota" OR "intestinal microbiome" OR "16S rRNA") AND ("Homo sapiens")`. For supervised analyses, inclusion required: (1) per-sample diagnosis labels (AD, MCI, or CN) accessible in public metadata; (2) raw FASTQ files or quality-filtered reads deposited in a public archive; (3) a minimum of 40 participants per cohort. Studies were excluded from supervised analyses if per-sample labels were IRB-restricted, if data had not been publicly released at the time of analysis, or if samples represented a tissue compartment other than fecal microbiome as the primary analyte; shotgun metagenomic cohorts meeting all other criteria were identified but excluded to avoid mixing amplicon and whole-genome sequencing platforms within the same analysis. Beyond the Kim/KBASE 2022 IRB restriction documented above, one further identified cohort was excluded on public-availability grounds: **ALBION 2025** (Maraki et al., *Journal of Alzheimer's Disease*, 2025; PRJNA1297934), a Greek cohort of 99 participants (50 MCI, 49 CN; V3–V4) — the only European 16S cohort identified in the search — was confirmed as the correct accession from the published paper, but NCBI BioProject reported the deposit as not yet publicly released as of the search date. ALBION was therefore excluded from the active cohort set; it remains the most direct candidate for a future sixth cohort should the data be released.

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
