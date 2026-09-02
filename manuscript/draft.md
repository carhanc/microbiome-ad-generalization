# Cross-Cohort Generalization Failure in Gut Microbiome-Based Alzheimer's Disease Classifiers: Evidence from Five Independent Cohorts with Batch Correction Analysis

**Running title:** Microbiome AD classifier generalization failure

**Authors:** Arhan Chakravarthy¹

**Affiliations:**
¹ Independent Researcher, Cupertino, CA

**Correspondence:** arhanc360@gmail.com

**Keywords:** gut microbiome; Alzheimer's disease; machine learning; cross-cohort generalization; batch correction; SHAP; PERMANOVA; 16S rRNA sequencing

---

## Abstract

Gut microbiome dysbiosis has been proposed as a biomarker for Alzheimer's disease (AD) and mild cognitive impairment (MCI), yet individual studies consistently report contradictory taxonomic signatures across populations, and cross-cohort classifier performance has received limited systematic evaluation. We assembled five publicly available 16S rRNA gene sequencing cohorts (509 samples; China, Kazakhstan, South Korea; Zhuang 2018, Ling 2020, Zhu 2022, Kaiyrlykyzy 2022, Kim 2022) and evaluated whether machine learning classifiers trained within one cohort could generalize to independent cohorts, whether standard batch correction methods could recover transferability, and which taxonomic features drove cohort-specific predictions. DADA2 amplicon sequence variant (ASV) processing and genus-level harmonization produced a unified matrix of 396 genera across all cohorts, transformed to centered log-ratio (CLR) coordinates for compositional analysis. The Zhu 2022 cohort deposited paired blood microbiome and fecal samples; only fecal samples (n=90) were retained for all analyses. Within four labeled cohorts (n=401 for binary AD/CN analyses), logistic regression and LightGBM classifiers achieved within-cohort AUC-ROC values ranging from 0.633 to 0.998 under nested cross-validation, broadly consistent with originally reported performance. In leave-one-cohort-out (LOCO) transfer experiments, LOCO AUC-ROC fell to 0.504–0.813, a mean drop of 0.176 AUC points for logistic regression and 0.146 for LightGBM, with one of eight model–cohort combinations collapsing to near chance (AUC = 0.504). PERMANOVA on Aitchison distance showed that cohort identity explained 19.3% of total Aitchison variance across all five cohorts (F=30.13, p=0.0001); in the four-cohort marginal model, cohort explained R²=0.172 versus only 0.014 for disease label—an approximately 12-fold disparity that is consistent with the interpretation that batch correction cannot recover generalization: both correction methods reduced mean LOCO AUC, with ComBat-seq consistently worsening and MMUPHin failing to improve transferability. SHAP analysis using out-of-fold values revealed near-zero pairwise Jaccard similarity between cohorts' top-20 predictive taxa (mean 0.135 for logistic regression, significantly above the random baseline of 0.027 [p=0.0019]; 0.058 for LightGBM, within the random baseline [p=0.068]) and twelve directional flip taxa for each model—genera whose association with AD versus cognitive normality reverses sign across cohorts, including *Akkermansia* and *Lactobacillus*. These results demonstrate that gut microbiome-based AD classifiers trained in one population do not transfer reliably to other populations, that the cohort effect is structurally much larger than the disease effect, and that linear batch correction methods are insufficient to bridge this gap. Multi-cohort harmonized training data and prospective validation across geographically and demographically diverse populations are prerequisites before microbiome-based AD classifiers can be considered clinically actionable.

---

## 1. Introduction

Alzheimer's disease (AD) is the leading cause of dementia worldwide, affecting an estimated 55 million people, with projections exceeding 150 million by 2050 (World Health Organization, 2023). The gut–brain axis has emerged as a mechanistic pathway linking intestinal microbial communities to neurodegeneration: altered gut microbiome composition has been associated with AD pathology through neuroinflammatory, metabolic, and neuroendocrine routes (Cryan et al., 2019; Kowalski and Mulak, 2019). Since the landmark observation that germ-free APPPS1 transgenic mice showed reduced amyloid-β pathology in the absence of gut microbiota compared with conventionally colonized APPPS1 controls (Harach et al., 2017), dozens of 16S rRNA and shotgun metagenomic studies have reported significant differences in gut microbial composition between AD or MCI patients and cognitively normal (CN) controls. MCI is included as a diagnostic category in this literature given its role as a preclinical or prodromal stage of AD, though the primary supervised analyses in this study focus on binary AD versus CN comparisons (see Section 2.4).

The clinical appeal of a gut microbiome-based diagnostic signature is considerable: stool is non-invasive, inexpensive to collect, and microbiome profiling is increasingly accessible. Multiple independent research groups have trained machine learning classifiers on gut microbiome abundance profiles and reported promising within-cohort discrimination performance, with area under the receiver operating characteristic curve (AUC-ROC) values of 0.70–0.90 in several published studies (Zhuang et al., 2018; Ling et al., 2021; Zhu et al., 2022). However, a critical property of any diagnostic biomarker—transferability across independent populations—has rarely been tested in this literature.

The limited cross-cohort comparisons that do exist reveal a striking inconsistency in the direction of reported associations. A systematic review and meta-analysis synthesizing 17 published studies (679 AD/MCI patients, 632 controls) documented that *Bacteroides* abundance was significantly elevated in North American cohorts yet significantly reduced in Chinese cohorts—a directional reversal across geographically distinct populations (Jemimah et al., 2023). Similar contradictions have been observed for *Akkermansia*, *Faecalibacterium*, and *Lactobacillus* across individual studies (Vogt et al., 2017; Liu et al., 2019; Zhuang et al., 2018). These population-level contradictions suggest that published within-cohort classifiers may be detecting cohort-specific compositional signatures—driven by diet, geography, host genetics, sequencing protocols, or a combination thereof—rather than a universal microbial correlate of AD pathology.

A recent multi-cohort study integrating 15 public AD microbiome projects identified candidate gut microbiota indicators across multiple geographic populations using integrated data (Wang et al., 2026), confirming that multi-population integration remains a methodological frontier. The present study is complementary, providing a systematic quantification of cross-cohort generalization failure and an empirical analysis of why batch correction methods fail to bridge the gap—questions that remain open even where multi-cohort classifiers have been attempted.

Despite this recognized heterogeneity, few studies have jointly quantified: (1) the magnitude of cross-cohort AUC collapse when trained classifiers are applied to independent cohorts; (2) the relative contribution of technical batch effects versus biological heterogeneity to this collapse, as estimated by variance partitioning; and (3) whether state-of-the-art microbiome batch correction methods can recover cross-cohort generalization. Answering these three questions jointly in one rigorous analysis is the primary contribution of this paper.

We assembled five publicly available 16S rRNA gene sequencing cohorts spanning three countries (China, Kazakhstan, South Korea) and conducted a comprehensive cross-cohort generalization evaluation. We report that: classifiers trained within one cohort consistently fail to generalize to other cohorts; cohort identity explains approximately 12 times more compositional variance than disease status; standard batch correction methods (ComBat-seq and MMUPHin) reduced mean cross-cohort AUC rather than improving transferability; and the specific taxa driving predictions are almost entirely cohort-specific and frequently contradict in direction across populations. Together, these findings constitute a systematic quantification of the generalization failure in the gut microbiome-based AD classifier literature, with direct implications for the field's trajectory toward clinical translation.

---

## 2. Materials and Methods

### 2.1 Cohort Identification and Data Acquisition

We performed a systematic search of NCBI SRA and ENA for publicly available 16S rRNA amplicon sequencing studies of gut microbiota in AD, MCI, or cognitively normal populations, deposited before January 2025. The search used the following terms in the SRA search interface: ("Alzheimer's disease" OR "Alzheimer disease" OR "mild cognitive impairment" OR "MCI" OR "dementia") AND ("gut microbiome" OR "gut microbiota" OR "intestinal microbiome" OR "16S rRNA") AND ("Homo sapiens"). Studies were included if they provided: (1) per-sample diagnosis labels (AD, MCI, or CN) accessible in public metadata; (2) raw FASTQ files or quality-filtered reads deposited in a public archive; (3) minimum sample sizes of 40 participants per cohort. We restricted this analysis to 16S rRNA studies; shotgun metagenomic cohorts were identified but excluded to avoid the additional platform confound introduced by mixing amplicon and whole-genome sequencing data. Cohorts were excluded if: per-sample labels were IRB-restricted (for supervised analyses); data had not been publicly released at the time of analysis; or samples represented a tissue compartment other than fecal microbiome as the primary analyte. A complete record of the search procedure, identified cohorts, and exclusion rationale is provided in Supplementary Note S1.

Four cohorts met all criteria for supervised analysis. A fifth cohort, Kim/KBASE 2022 (PRJEB50447), lacked publicly available per-sample diagnosis labels and was retained only for unsupervised cohort-level analyses. The ALBION 2025 Greek cohort (PRJNA1297934) was identified but excluded because data had not yet been publicly released at the time of analysis. Group-level counts from the published paper (18 amyloid-positive preclinical AD; 60 cognitively normal) were available but insufficient for supervised classifier training or evaluation.

The four cohorts with per-sample diagnosis labels used in supervised analyses were:

**Zhuang 2018 (China; PRJNA554111):** 86 participants (43 AD, 43 CN) from the Third Military Medical University, Chongqing. Paired-end MiSeq 2×300 bp, V3–V4 region (Zhuang et al., 2018).

**Ling 2020 (China; PRJNA633959):** 171 participants (100 AD, 71 CN). Paired-end MiSeq 2×300 bp, V3–V4 region (Ling et al., 2021).

**Zhu 2022 (China; PRJNA489760):** 180 samples deposited in SRA from a larger 302-participant study, representing 90 participants (30 AD, 30 MCI, 30 CN) each contributing 2 sample types = 180 samples. Each participant contributed two samples: a fecal (F_*) sample and a blood microbiome (B_*) sample. Only the 90 fecal samples (30 AD, 30 MCI, 30 CN) were used in all analyses; the 90 blood microbiome samples were excluded as they represent a distinct biological compartment and would confound within-cohort classifiers. Binary analyses used 30 AD and 30 CN fecal samples (n=60). Paired-end MiSeq 2×300 bp, V3–V4 (Zhu et al., 2022). Diagnostic criteria across cohorts relied on clinical AD diagnosis (Zhuang, Ling, Zhu) or amyloid-PET positivity (Kim); these differences in diagnostic ascertainment are acknowledged as a source of between-cohort biological heterogeneity.

**Kazakhstan 2022 (PRJNA811324):** 84 participants (41 AD, 43 CN). Single-end NovaSeq sequencing, V3–V4 region (Kaiyrlykyzy et al., 2022). This cohort is the only single-end dataset in the study; implications for read length and taxonomic resolution are discussed in Section 2.2 and the Limitations section.

The Kim/KBASE 2022 cohort (Jung et al., 2022; PRJEB50447; 78 participants; 18 preclinical AD by amyloid-PET, 60 CN) was included in PERMANOVA analyses (Phase 4) only and is explicitly excluded from all supervised classifier training and evaluation.

All raw FASTQ data were downloaded from NCBI SRA using fasterq-dump (SRA Toolkit v3.0.10). Cohort-level sequencing statistics are reported in Supplementary Table S2.

### 2.2 DADA2 Amplicon Sequence Variant Processing

Raw paired-end reads (four cohorts) and single-end reads (Kazakhstan) were processed independently per cohort using DADA2 v1.28 in R v4.3 (Callahan et al., 2016), to avoid introducing cross-cohort data leakage during error model estimation.

For paired-end cohorts, reads were trimmed to truncation lengths selected per cohort based on per-base quality score profiles (generally 240 bp forward / 200 bp reverse for V3–V4 MiSeq 2×300 data), filtered to remove reads with >2 expected errors, chimeras were removed using the consensus method, and forward and reverse reads were merged. For Kazakhstan (single-end), the single-end DADA2 workflow was applied without merging; truncation length was adjusted to 250 bp after inspection of quality profiles (see Supplementary Note S2 for rationale). The mixed single-end/paired-end design is an additional technical source of between-cohort heterogeneity acknowledged as a study limitation.

Taxonomic assignment used the SILVA reference database v138.1. ASVs were collapsed to genus level by summing raw counts per sample. A genus was retained across the entire study if it was observed at ≥1 count in ≥2 samples from at least one cohort; this threshold was applied per cohort before cross-cohort union. The union of retained genera across all five cohorts produced a unified feature set of 396 genera. Genera absent in a given sample received a count of zero.

### 2.3 Compositional Data Transformation

Raw genus-level count tables were converted to centered log-ratio (CLR) coordinates for all downstream analyses requiring distributional assumptions on unconstrained data (machine learning classifiers, PERMANOVA on Euclidean distance). The CLR transformation is defined as:

$$\text{clr}(x_i) = \log\left(\frac{x_i}{g(\mathbf{x})}\right)$$

where $g(\mathbf{x})$ is the geometric mean of the composition $\mathbf{x}$.

Zero counts were replaced prior to transformation using the multiplicative replacement strategy (Martín-Fernández et al., 2003) with $\delta = 0.65$:

$$x_i^* = \begin{cases} \delta & \text{if } x_i = 0 \\ \left(1 - n_0 \cdot \delta / T\right) \cdot x_i & \text{if } x_i > 0 \end{cases}$$

where $T = \sum_j x_j$ is the total read count of the sample and $n_0$ is the count of zeros in the sample. The unified genus matrix contained 75.4% zero entries (151,926 of 201,564 cells across the 509 fecal samples and 396 genera). The value $\delta = 0.65$ was chosen such that the replacement mass remained small relative to typical per-sample sequencing depth (median 29,521 reads): in the worst-observed sample, the total mass allocated to zero-filled entries was 2.53% of that sample's total reads. A built-in safeguard rescales $\delta$ downward whenever $n_0 \cdot \delta \geq T$, capping replacement mass at 1% of sample total in such cases. This strategy preserves the ratios among non-zero components while introducing a minimal perturbation proportional to the relative frequency of zeros, and is preferable to the commonly used additive pseudocount which distorts low-abundance taxa more severely.

CLR-transformed values were used as features for all machine learning models and for Aitchison distance computation (Euclidean distance on CLR values). Raw genus-level count matrices (not CLR-transformed) were provided as input to ComBat-seq, which requires non-negative integer counts. Relative abundances (proportions) were provided to MMUPHin. Bray-Curtis dissimilarity was computed on relative abundances. CLR values were used for all other analyses.

### 2.4 Within-Cohort Classifier Baselines

We trained binary classifiers (AD vs. CN) independently for each of the four labeled cohorts using two model families:

**Logistic Regression with Elastic Net regularization** (scikit-learn v1.5; Pedregosa et al., 2011; LogisticRegression, solver="saga", L1-ratio=0.5). Regularization strength C was tuned over {0.01, 0.1, 1.0, 10.0} in inner cross-validation.

**LightGBM gradient-boosted trees** (lightgbm v4; Ke et al., 2017). Hyperparameters tuned: number of estimators ∈ {100, 200, 500}, learning rate ∈ {0.05, 0.1, 0.2}, max depth ∈ {3, 5, -1}.

For Zhu 2022, AD was defined as the positive class and MCI samples were excluded from binary analyses (resulting in n=60: 30 AD, 30 CN fecal samples); a separate three-class analysis (AD / MCI / CN) used all 90 fecal samples and is reported as a secondary finding. Because only fecal samples were retained, each patient appears exactly once; standard StratifiedKFold was used without any paired-sample grouping structure.

Performance was estimated using nested cross-validation: a 10-fold outer loop for performance estimation and a 5-fold inner loop for hyperparameter tuning. The primary metric was AUC-ROC computed on out-of-fold predictions pooled across the 10 outer folds, providing a single unbiased AUC estimate per cohort–model combination. Bootstrap 95% confidence intervals (1,000 resamples) on the pooled out-of-fold predictions are reported throughout. We also report sensitivity and specificity at the Youden's J optimal threshold. For comparison with originally reported performance, the published AUC values from the original cohort papers are noted in Table 2 where available.

### 2.5 Cross-Cohort Generalization Experiment

We evaluated two cross-cohort transfer designs:

**Leave-one-cohort-out (LOCO):** For each test cohort, a classifier was trained on the pooled samples from the remaining three labeled cohorts (no cross-validation within this training step; hyperparameters were selected via 5-fold inner CV on the training pool). The trained model was then evaluated on the held-out cohort. This produced four LOCO AUC estimates per model type. We report the LOCO AUC drop as the difference between within-cohort (nested CV) AUC and LOCO AUC; bootstrap 95% CIs are reported for each LOCO AUC estimate, and the non-overlap of within-cohort and LOCO CIs serves as an indicator of a statistically reliable drop.

**Pairwise transfer:** A classifier was trained on one cohort and tested on each of the other three cohorts individually, yielding a 4×4 matrix (diagonal excluded) of cross-cohort AUC values per model type. This design is more stringent than LOCO and isolates specific directional failures.

Bootstrap 95% CIs (1,000 resamples) on the test set predictions are reported for all LOCO and pairwise estimates.

**Sensitivity analyses:** To assess whether generalization failure was attributable to any single cohort, we repeated the LOCO experiment under three additional configurations: (A) excluding Kazakhstan entirely, with LOCO evaluated on the remaining three cohorts; (B) excluding Zhuang 2018; (C) excluding Ling 2020. Results are reported in Section 3.3 and Supplementary Figure S4.

### 2.6 PERMANOVA Variance Decomposition

To quantify the relative contribution of cohort identity and disease status to total compositional variance, we conducted PERMANOVA (permutational multivariate analysis of variance) using the adonis2 function in the R package vegan v2.6 (Oksanen et al., 2024). Two distance matrices were used:

- **Aitchison distance** (Euclidean distance on CLR values): the compositionally appropriate distance that does not distort due to the simplex constraint.
- **Bray-Curtis dissimilarity** (on relative abundances): included for comparability with the broader microbiome literature.

All PERMANOVA models used 9,999 permutations with Type III (marginal) sums of squares (by="margin") to allow independent assessment of cohort and diagnosis effects. The minimum attainable permutation p-value with 9,999 permutations is reported as p=0.0001 throughout. Three models were fitted: (1) cohort only, all five cohorts (n=509); (2) cohort + diagnosis, restricted to the four labeled cohorts (n=401: AD and CN samples only, Zhu 2022 MCI excluded), with marginal R² for each term; (3) within each labeled cohort, diagnosis as the sole predictor (within-cohort disease effect size). Note that Figure 4B displays the PCoA from Model 1 (n=509, all five cohorts including MCI samples), while the marginal R² values reported throughout the text derive from Model 2 (n=401).

Homogeneity of multivariate dispersions was tested using betadisper and permutest (9,999 permutations) to assess whether observed cohort differences reflected centroid separation, within-group spread heterogeneity, or both. A significant beta-dispersion test indicates that PERMANOVA R² values partially reflect variation in within-group spread, not only centroid shifts.

To further quantify the degree to which cohort-specific information is retained in the preprocessed feature space, we trained a logistic regression classifier (one-vs-rest, stratified 5-fold CV) to predict cohort identity from the CLR-transformed features across all 509 samples. The macro-AUC for cohort prediction provides a complementary measure of cohort separability.

### 2.7 Batch Correction

We evaluated one count-based batch-correction method (ComBat-seq) and one microbiome-specific meta-analysis framework (MMUPHin), applying each with cohort identity as the batch variable and disease status as the biological covariate of interest:

**ComBat-seq** (Zhang et al., 2020): A negative binomial regression-based batch correction method adapted from the original ComBat algorithm for use with RNA-seq count data. We provided raw genus-level count matrices (genera × samples, four labeled cohorts only) and specified group = disease_status to preserve the biological variable. Corrected count tables were CLR-transformed using the same multiplicative zero-replacement procedure described in Section 2.3.

**MMUPHin** v1.15 (Ma et al., 2022): A meta-analysis framework for microbiome data that fits linear mixed-effects models to remove batch effects from relative abundance profiles while adjusting for covariates. Relative abundance profiles from the four labeled cohorts were provided, with batch = cohort and covariates = diagnosis.

After correction, the LOCO generalization experiment (Section 2.5) was repeated on the batch-corrected CLR matrices. We also re-computed within-cohort AUC on corrected data to assess whether correction degraded within-cohort biological signal. Importantly, both correction methods were fitted jointly on all four labeled cohorts, including the held-out test cohort and its diagnosis labels—a transductive design in which test-set information is used during correction. This means the corrected LOCO results do not represent a valid prospective external-validation scenario: providing test-cohort labels may improve or worsen performance relative to a strictly prospective design in unpredictable ways. The corrected LOCO values should be interpreted as the performance achievable under a transductive batch correction scenario, not as an estimate of prospective deployability.

### 2.8 SHAP Feature Importance Analysis

To identify which genera drove predictions in each cohort and whether these features were consistent across cohorts, we computed SHapley Additive exPlanations (SHAP) values (Lundberg and Lee, 2017) using an out-of-fold (OOF) procedure to avoid the in-sample overfitting bias that arises when SHAP values are computed on the same data used to fit the model.

For each labeled cohort, SHAP values were computed using the same 10-fold outer CV structure used for AUC estimation (Section 2.4). In each fold, the classifier was trained on the 9 training folds; SHAP values were then computed for the held-out fold using the training fold as the SHAP background distribution. For logistic regression, shap.LinearExplainer was used with the training-fold data as background (feature_perturbation="interventional"), and the held-out fold as the evaluation set; this ensures the mean SHAP value over the evaluation set is not algebraically forced to zero (as would occur if background and evaluation sets were identical). For LightGBM, shap.TreeExplainer was used (check_additivity=False) with the same OOF structure. SHAP values were aggregated across all 10 folds: the mean absolute SHAP value per genus (|SHAP|) was used as the importance metric, and the mean signed SHAP value per genus across all samples indicated direction of association (positive = AD-associated; negative = CN-associated).

For LOCO analyses, SHAP values were computed on the held-out test cohort from each LOCO-trained model, with the three-cohort training pool as the SHAP background.

Pairwise taxonomic overlap between cohorts was quantified using Jaccard similarity on the top-N SHAP genera (N ∈ {10, 20, 50}); N=20 is the primary threshold, chosen as a balance between capturing the most informative features and avoiding the long tail of near-zero-importance genera. A null distribution for Jaccard similarity was estimated by randomly sampling 20 genera (without replacement) from the full 396-genus set for each cohort pair and computing Jaccard similarity; this was repeated 10,000 times per model. The analytical expectation for two random top-20 sets from 396 genera is E[J] = 20/(2×396−20) ≈ 0.026. The empirical p-value was computed as the proportion of permutations with Jaccard ≥ observed mean. A genus was classified as a "directional flip" if it appeared in the top-20 SHAP features of at least two cohorts and had opposite sign (positive mean signed SHAP across all samples in one cohort; negative in another), constituting an empirical cross-population directional contradiction. SHAP directions should be interpreted as model-derived feature contributions under the fitted model, not as direct estimates of biological association or its statistical significance; directions may be influenced by regularization, collinearity among CLR features, and nonlinear interactions (for LightGBM).

### 2.9 AI Tool Disclosure

Initial manuscript drafts and computational scripts were prepared with minimal assistance from Claude (Anthropic). All analytical decisions, parameter choices, experimental design, statistical interpretations, and scientific conclusions were made by the author. All code and text were reviewed, validated, and substantially revised under the author's direction.

### 2.10 Statistical Analysis Environment

All Python analyses used scikit-learn v1.5, lightgbm v4, shap v0.51, pandas v2.x, and numpy v1.x in a conda environment (Python 3.11). R analyses used DADA2 v1.28, vegan v2.6, sva v3.50 (ComBat-seq), and MMUPHin v1.15 in R v4.3. Figures were generated with matplotlib v3.8 at 300 DPI.

---

## 3. Results

### 3.1 Cohort Characteristics and Processing

Five cohorts comprising 599 deposited samples from three countries were assembled (Figure 1A; Table 1). After DADA2 processing and chimera removal, read retention rates varied substantially across cohorts (47.9%–78.3%; Supplementary Table S2), reflecting differences in sequencing platform, library preparation, and amplicon length. The Kazakhstan cohort (single-end, NovaSeq) had notably different quality dynamics compared to the paired-end MiSeq cohorts. SILVA v138.1 genus-level classification and cross-cohort union produced a unified feature matrix of 396 genera. The 90 blood microbiome (B_*) samples from Zhu 2022 were excluded prior to all analyses (Section 2.4), leaving a working dataset of 509 samples across five cohorts. CLR transformation was applied to all 509 samples jointly after union, with multiplicative zero-replacement at δ=0.65.

For supervised classifier analyses (Phases 2–3 and 5–6), the Kim/KBASE 2022 cohort was excluded because per-sample diagnosis labels are IRB-restricted (Supplementary Note S1), leaving four labeled cohorts with a total of 401 samples for binary AD/CN analyses (86 Zhuang + 171 Ling + 60 Zhu fecal binary subset + 84 Kazakhstan; Zhu 2022 MCI samples excluded from binary analyses).

### 3.2 Within-Cohort Classification Performance

Within-cohort nested cross-validation yielded AUC-ROC values broadly consistent with those reported in the original publications (Figure 2; Table 2). Logistic regression achieved AUC-ROC of 0.633 [95% CI: 0.51–0.75] for Zhuang 2018, 0.865 [0.81–0.92] for Ling 2020, 0.998 [0.99–1.00] for Zhu 2022, and 0.766 [0.66–0.86] for Kazakhstan. LightGBM yielded similar values: 0.635 [0.52–0.76], 0.862 [0.81–0.91], 0.979 [0.95–1.00], and 0.694 [0.57–0.80] respectively.

The near-perfect within-cohort AUC for Zhu 2022 fecal samples (logistic regression 0.998; LightGBM 0.979; n=60 binary) is notable given the small sample size (30 per class). With 10-fold CV on 60 samples, each test fold contains approximately six samples; fold-level AUC estimates can be unstable at this scale. PERMANOVA analysis (Section 3.4) reveals that diagnosis explains R²=0.094 of within-cohort compositional variance (F=6.01, p=0.0001), suggesting statistically detectable diagnosis-associated compositional differences in this cohort. However, the high AUC should be interpreted cautiously and the LOCO experiment (Section 3.3) provides a more informative assessment of how transferable this signal is.

The three-class analysis (AD / MCI / CN) on all 90 Zhu 2022 fecal samples using logistic regression (one-vs-rest macro-AUC under nested cross-validation) yielded macro-AUC=0.761, confirming the multi-class discriminability of this cohort's fecal compositional data.

For Zhuang 2018, within-cohort AUC values (0.633–0.635) were modest, consistent with the small sample size and weak within-cohort diagnosis effect observed in PERMANOVA (R²=0.017, p=0.066; Section 3.4).

**Table 1.** Cohort characteristics.

| Cohort | Country | Accession | N (analyzed) | AD | CN | MCI | Sequencing | Library | Supervised Analysis |
|---|---|---|---|---|---|---|---|---|---|
| Zhuang 2018 | China | PRJNA554111 | 86 | 43 | 43 | — | MiSeq 2×300 | Paired-end | Phases 2–6 |
| Ling 2020 | China | PRJNA633959 | 171 | 100 | 71 | — | MiSeq 2×300 | Paired-end | Phases 2–6 |
| Zhu 2022† | China | PRJNA489760 | 90 (fecal) | 30 | 30 | 30 | MiSeq 2×300 | Paired-end | Phases 2–6 |
| Kaiyrlykyzy 2022 | Kazakhstan | PRJNA811324 | 84 | 41 | 43 | — | NovaSeq | Single-end | Phases 2–6 |
| Kim 2022 | South Korea | PRJEB50447 | 78 | 18* | 60 | — | MiSeq 2×300 | Paired-end | Phase 4 only |
| **Total** | | | **509** | **232** | **247** | **30** | | | |

*Amyloid-PET positive preclinical AD; per-sample labels IRB-restricted and not publicly available.
†Zhu 2022 deposited 180 samples (90 fecal F_* + 90 blood B_*); only fecal samples retained for analysis.

**Table 2.** Within-cohort nested cross-validation performance (10-fold outer, 5-fold inner CV).

| Cohort | Model | N | AUC-ROC | 95% CI | Sensitivity | Specificity | Originally Reported AUC |
|---|---|---|---|---|---|---|---|
| Zhuang 2018 | LogReg | 86 | 0.633 | [0.51–0.75] | 0.814 | 0.512 | ~0.73 (Zhuang et al., 2018) |
| Zhuang 2018 | LGBM | 86 | 0.635 | [0.52–0.76] | 0.884 | 0.442 | — |
| Ling 2020 | LogReg | 171 | 0.865 | [0.81–0.92] | 0.780 | 0.859 | ~0.87 (Ling et al., 2021) |
| Ling 2020 | LGBM | 171 | 0.862 | [0.81–0.91] | 0.740 | 0.887 | — |
| Zhu 2022 | LogReg | 60 | 0.998 | [0.99–1.00] | 0.967 | 1.000 | ~0.94 (Zhu et al., 2022) |
| Zhu 2022 | LGBM | 60 | 0.979 | [0.95–1.00] | 0.900 | 1.000 | — |
| Kaiyrlykyzy 2022 | LogReg | 84 | 0.766 | [0.66–0.86] | 0.829 | 0.605 | ~0.78 (Kaiyrlykyzy et al., 2022) |
| Kaiyrlykyzy 2022 | LGBM | 84 | 0.694 | [0.57–0.80] | 0.610 | 0.767 | — |

N=60 for Zhu 2022 binary analysis (30 AD + 30 CN fecal samples; MCI and blood samples excluded). Originally reported AUC values are approximate, derived from the primary publications using the best-performing classifier reported; direct numerical comparison is limited by differences in cohort subset, model type, and validation design.

### 3.3 Cross-Cohort Generalization Failure

LOCO cross-cohort AUC-ROC values revealed substantial and variable degradation in cross-cohort performance across all cohorts and model types (Figure 3A; Table 3). Mean LOCO AUC across the four cohorts was 0.640 for logistic regression and 0.647 for LightGBM—corresponding to mean drops of 0.176 and 0.146 AUC points respectively from within-cohort performance. One of eight model–cohort combinations produced LOCO AUC near chance: Zhuang 2018 with logistic regression (AUC=0.504 [0.37–0.63]).

The largest absolute logistic regression AUC drop was observed for Kazakhstan (within-cohort 0.766 → LOCO 0.564, Δ=0.202). Zhu 2022 showed the largest LightGBM drop: within-cohort 0.979 → LOCO 0.758 (Δ=0.221); its logistic regression drop was also substantial (0.998 → 0.813, Δ=0.185). These results are consistent with LOCO models relying on predictive patterns learned from training cohorts that do not transfer fully to the held-out cohort, where cohort-specific compositional structure may contribute to within-cohort performance.

Ling 2020 also experienced substantial LOCO collapse: logistic regression 0.865 → 0.678 (Δ=0.187); LightGBM 0.862 → 0.674 (Δ=0.188). LightGBM confidence intervals for the Ling 2020 LOCO estimate [0.59–0.75] do not overlap with the within-cohort estimate [0.81–0.91], suggesting a performance drop unlikely to be explained by sampling variability alone.

LOCO performance was heterogeneous across cohorts: Zhu 2022 achieved the highest LOCO AUC (logistic regression 0.813), while Zhuang 2018 collapsed to near chance. It should be noted that only one of the eight main model–cohort combinations fell to near chance, and five of eight exceeded 0.55; the results therefore demonstrate substantial and variable degradation rather than uniform failure.

**Sensitivity analyses** confirmed that generalization failure was not attributable to Kazakhstan, Zhuang 2018, or Ling 2020 individually (Supplementary Figure S4). When Kazakhstan was excluded and LOCO was performed on the three remaining Chinese cohorts, Zhuang 2018 remained near chance (LogReg 0.512, LightGBM 0.567), demonstrating that the generalization failure for Zhuang 2018 is not a Kazakhstan artifact. When Zhuang 2018 was excluded, LOCO on the remaining three cohorts still showed AUC values that were substantially below within-cohort baselines for Kazakhstan (LogReg 0.599, LightGBM 0.495). When Ling 2020 was excluded, Kazakhstan LOCO LogReg was 0.604. These sensitivity results confirm the robustness of the main finding across cohort compositions.

**Table 3.** LOCO cross-cohort AUC-ROC versus within-cohort nested CV baseline.

| Test Cohort | Model | Within-Cohort AUC | LOCO AUC | LOCO 95% CI | AUC Drop |
|---|---|---|---|---|---|
| Zhuang 2018 | LogReg | 0.633 | 0.504 | [0.37–0.63] | 0.129 |
| Zhuang 2018 | LGBM | 0.635 | 0.565 | [0.44–0.69] | 0.070 |
| Ling 2020 | LogReg | 0.865 | 0.678 | [0.59–0.76] | 0.187 |
| Ling 2020 | LGBM | 0.862 | 0.674 | [0.59–0.75] | 0.188 |
| Zhu 2022 | LogReg | 0.998 | 0.813 | [0.69–0.92] | 0.185 |
| Zhu 2022 | LGBM | 0.979 | 0.758 | [0.63–0.87] | 0.221 |
| Kazakhstan | LogReg | 0.766 | 0.564 | [0.44–0.69] | 0.202 |
| Kazakhstan | LGBM | 0.694 | 0.591 | [0.47–0.71] | 0.103 |
| **Mean** | **LogReg** | **0.816** | **0.640** | — | **0.176** |
| **Mean** | **LGBM** | **0.793** | **0.647** | — | **0.146** |

Pairwise single-cohort transfer AUC values (Figure 3B–C; Supplementary Table S1) ranged from 0.54 to 0.88 for logistic regression. The highest pairwise AUC was Ling 2020 → Zhu 2022 (logistic regression AUC=0.88), reflecting transfer between two large Chinese fecal cohorts with overlapping compositional structure. Kazakhstan → Ling 2020 also transferred reasonably (logistic regression AUC=0.74). Most other pairwise AUCs clustered between 0.54 and 0.70, consistent with the overall pattern of limited cross-cohort transferability.

### 3.4 PERMANOVA Variance Decomposition

PERMANOVA on Aitchison distance across all five cohorts showed that cohort identity explained 19.3% of total compositional variance (R²=0.193, F=30.13, p=0.0001, 9,999 permutations; Figure 4A). On Bray-Curtis dissimilarity, the cohort R² was 0.136 (F=19.80, p=0.0001). The cohort effect was substantially larger than disease signal: in a marginal two-variable model restricted to the four labeled cohorts (n=401), cohort (marginal) explained R²=0.172 (F=28.03, p=0.0001) versus diagnosis (marginal) R²=0.014 (F=6.82, p=0.0001)—an approximately 12-fold disparity within the same model. This disparity is consistent with the interpretation that batch correction methods struggle to recover cross-cohort generalization, as discussed in Section 3.5.

A logistic regression classifier trained to predict cohort identity from the CLR-transformed features across all 509 samples achieved macro-AUC=0.9917 (one-vs-rest, stratified 5-fold CV; per-cohort AUC: Kazakhstan 1.000, Kim 2022 0.997, Ling 2020 0.993, Zhuang 2018 0.989, Zhu 2022 0.980). Cohort membership is thus near-perfectly recoverable from the same features used for disease classification, confirming that cohort-specific information dominates the compositional feature space.

Principal coordinates analysis (PCoA) on Aitchison distance (Figure 4B; n=509) visually confirmed strong cohort-level clustering, with partial overlap between Chinese cohorts (Zhuang 2018, Ling 2020) but clear separation of Kazakhstan and Zhu 2022. The Kim/KBASE 2022 cohort was clearly separated from all others along PC1, consistent with its preclinical amyloid-based staging, South Korean origin, and potentially distinct dietary patterns.

Beta-dispersion testing revealed significant heterogeneity in within-group spread across cohorts (F=13.93, p=0.0001, 9,999 permutations; Figure 4C), meaning that the PERMANOVA R² of 0.193 reflects both compositional centroid shifts and differences in within-cohort spread. Kazakhstan showed the greatest within-cohort variance, while the Chinese cohorts (Zhuang, Ling) were more tightly dispersed. This dispersion heterogeneity is relevant to interpreting PERMANOVA effect sizes and is acknowledged in the Methods section.

Within-cohort PERMANOVA revealed highly variable disease-associated compositional effects: Ling 2020 (R²=0.053, F=9.53, p=0.0001) and Zhu 2022 (R²=0.094, F=6.01, p=0.0001) showed significant diagnosis effects, while Zhuang 2018 (R²=0.017, p=0.066) and Kazakhstan (R²=0.022, p=0.069) did not reach statistical significance at α=0.05. These within-cohort disease effect sizes broadly rank-order with within-cohort AUC performance, providing a triangulating validation of the classifier results. The Zhu 2022 fecal-only R²=0.094 represents the highest within-cohort disease effect among the four labeled cohorts, consistent with a detectable diagnosis-associated compositional difference in this dataset.

### 3.5 Batch Correction Does Not Recover Cross-Cohort Generalization

Batch correction was applied using ComBat-seq (on raw counts) and MMUPHin (on relative abundances) to the four labeled cohorts with cohort as the batch variable (Figure 5; Table 4). Both methods were applied in a transductive design in which the held-out test cohort and its labels were included during correction parameter estimation. Both methods failed to recover cross-cohort generalization, and ComBat-seq substantially worsened it.

**Table 4.** Mean LOCO AUC across four cohorts, by correction method and model.

| Method | LogReg Mean LOCO AUC | LogReg Mean Drop from Within-Cohort | LGBM Mean LOCO AUC | LGBM Mean Drop |
|---|---|---|---|---|
| Uncorrected | 0.640 | 0.176 | 0.647 | 0.146 |
| ComBat-seq | 0.493 | 0.322 | 0.555 | 0.238 |
| MMUPHin | 0.612 | 0.203 | 0.582 | 0.210 |

ComBat-seq reduced mean logistic regression LOCO AUC from 0.640 to 0.493 and degraded all four individual cohort LOCO AUCs, with the Kazakhstan cohort collapsing to AUC=0.308 [0.20–0.42]. LightGBM was somewhat more robust: ComBat-seq reduced LGBM mean LOCO from 0.647 to 0.555. MMUPHin was less destructive but still failed to improve upon uncorrected LOCO performance (LogReg mean: 0.640 → 0.612; LGBM: 0.647 → 0.582).

Critically, within-cohort AUC on corrected data was largely preserved (within ±0.05 of uncorrected) for most cohort–method combinations, ruling out the hypothesis that batch correction was simply removing all signal indiscriminately. The exception was Kazakhstan, where both correction methods moderately improved within-cohort AUC (logistic regression: 0.766 → 0.784 with ComBat-seq, 0.766 → 0.788 with MMUPHin), while simultaneously collapsing LOCO AUC to below chance. This pattern—improved within-cohort AUC coupled with worsened cross-cohort AUC—is consistent with the correction algorithm overfitting the batch-correction transformation to Kazakhstan's within-cohort structure, making its representation less compatible with the other cohorts, not more.

The pattern is consistent with the interpretation that when cohort explains approximately 12 times more compositional variance than disease status (marginal R²: 0.172 vs. 0.014), linear correction methods that target cohort-associated variance will tend to also remove disease-associated variance, because those sources of variation are not cleanly separable in this composition space. This interpretation should be understood as a plausible mechanistic account, not as a formal proof that cohort and disease axes overlap; a large cohort R² could, in principle, coexist with an orthogonal disease axis. However, the empirical worsening of LOCO AUC under correction—including when correction was applied with the test cohort's labels available—is consistent with the non-orthogonality account.

### 3.6 SHAP Analysis Reveals Cohort-Specific and Contradictory Taxonomic Signatures

To characterize the feature-level patterns underlying generalization failure, we computed out-of-fold SHAP values for within-cohort classifiers and LOCO-trained classifiers (Figure 6). Using an OOF design avoids the in-sample overfitting bias identified in prior SHAP analyses of microbiome classifiers, where training-set SHAP values can reflect model memorization rather than the model's out-of-sample feature importance.

**Minimal overlap in top predictive taxa across cohorts.** Pairwise Jaccard similarity of the top-20 SHAP genera between cohort pairs averaged 0.135 for logistic regression and 0.058 for LightGBM (Figure 6C). Against a permutation null distribution (10,000 draws of two random top-20 sets from 396 genera; null mean 0.027 [95% CI: 0.000–0.081]), logistic regression overlap was significantly above chance (p=0.0019), whereas LightGBM overlap was within the null 95% CI (p=0.068). These results confirm that within-cohort logistic regression classifiers share a small but statistically detectable set of top predictive genera across cohorts, whereas LightGBM concentrates importance in sets that are essentially non-overlapping—consistent with the greater sparsity of tree model feature selection.

**Cohort-specific top predictors.** The most discriminative genus for logistic regression was *Subdoligranulum* (CN-associated) in Zhuang 2018, while *Akkermansia* was the top predictor in the remaining three cohorts—CN-associated in Ling 2020 (|SHAP|=0.270) but AD-associated in Zhu 2022 (|SHAP|=0.228) and Kazakhstan (|SHAP|=0.136)—illustrating that even the single most important taxon can reverse direction across populations. For LightGBM, the top taxon per cohort was entirely distinct: *Lachnoclostridium* (Zhuang 2018), *Akkermansia* (Ling 2020), *Bacteroides* (Zhu 2022), *Castellaniella* (Kazakhstan) (Figure 6A; Supplementary Figure S2).

**Directional flip taxa.** Among taxa appearing in the top-20 for two or more cohorts, twelve showed directional reversal for logistic regression and twelve for LightGBM (Figure 6B), confirming that directional inconsistency is pervasive across both model types under the OOF evaluation.

*Logistic regression flips (12 taxa):* Akkermansia (CN-associated in Zhuang/Ling; AD-associated in Zhu 2022/Kazakhstan), Bifidobacterium (AD-associated in Zhuang/Zhu 2022/Kazakhstan; CN-associated in Ling), Christensenellaceae R-7 group (CN-associated in Zhuang/Ling/Zhu 2022; AD-associated in Kazakhstan), Coprococcus (CN-associated in Zhuang/Ling/Kazakhstan; AD-associated in Zhu 2022), Eisenbergiella, Enterococcus, Lactobacillus (CN-associated in Chinese cohorts; AD-associated in Kazakhstan), NK4A214 group (CN-associated in Zhuang/Ling; AD-associated in Zhu 2022/Kazakhstan), Parasutterella (CN-associated in Zhuang/Kazakhstan; AD-associated in Ling/Zhu 2022), Romboutsia (AD-associated in Zhuang/Kazakhstan; CN-associated in Ling/Zhu 2022), Ruminococcus gnavus group, and Subdoligranulum.

*LightGBM flips (12 taxa):* Bacteroides, Akkermansia, Incertae Sedis, Faecalibacterium, Romboutsia, Oscillibacter, Lachnospira, Lactobacillus (CN-associated in Ling/Zhu 2022; AD-associated in Kazakhstan), [Eubacterium] siraeum group, Roseburia, Agathobacter, and Colidextribacter.

*Akkermansia* shows a complex pattern: for logistic regression it is CN-associated in Zhuang 2018 and Ling 2020 but AD-associated in Zhu 2022 and Kazakhstan, constituting a directional flip across cohorts. For LightGBM, *Akkermansia* is the dominant predictor in Ling 2020 (|SHAP|=0.947, CN-associated) and also appears in other cohorts with mixed directionality. *Akkermansia muciniphila* abundance is strongly influenced by host diet, particularly fiber and polyphenol intake (Plovier et al., 2017), and by geographic-population-specific dietary patterns that likely modulate its abundance independently of any AD-specific effect. *Lactobacillus* shows a consistent pattern: CN-associated in Chinese cohorts but AD-associated in Kazakhstan across both model types, plausibly reflecting the prominence of fermented dairy (kumis, kefir) in Central Asian diets. These dietary explanations are speculative in the absence of dietary co-variate data and are offered as mechanistic hypotheses rather than established conclusions.

**LOCO SHAP analysis: no single taxon transfers universally.** When SHAP values were computed for LOCO-trained classifiers applied to held-out cohorts, the top predictor was *Akkermansia* (CN-associated) for two of four held-out cohorts—Zhuang 2018 (|SHAP|=0.470) and Zhu 2022 (|SHAP|=0.372)—but *Subdoligranulum* (AD-associated, |SHAP|=0.269) for held-out Ling 2020 and *Christensenellaceae R-7 group* (AD-associated, |SHAP|=0.368) for held-out Kazakhstan. No genus ranked first across all four held-out conditions, and the recurrence of Akkermansia in two conditions does not indicate a stable transferable signal: it is the top predictor for two cohorts while playing a directionally inconsistent role in the underlying within-cohort analyses (Section 3.6), consistent with the classifier activating whichever genera most closely resemble its training-cohort signature rather than a disease marker that generalizes across populations.

---

## 4. Discussion

### 4.1 Systematic Quantification of Generalization Failure

This study provides a systematic, quantitative analysis of cross-cohort generalization failure in gut microbiome-based AD classifiers. Across four labeled cohorts, two model types, two cross-cohort evaluation designs, and two batch correction methods, we consistently observe the same result: classifiers trained within one population fail to generalize reliably to another. The mean LOCO AUC drop of 0.176 points for logistic regression and 0.146 for LightGBM, with one of eight model–cohort combinations near chance, demonstrates that previously reported within-cohort performance figures cannot be taken as evidence of a generalizable diagnostic signal.

This finding is consistent with the broader literature. The *Bacteroides* directionality reversal documented by Jemimah et al. (2023) across North American and Chinese cohorts, the contradictory reports for *Akkermansia* (Liu et al., 2020), and the highly cohort-specific patterns reported in individual papers all suggest population-level heterogeneity. A recent multi-cohort integration study (Wang et al., 2026, mSystems) similarly identified cross-geographic heterogeneity as a major obstacle; the present analysis complements that work by providing a systematic quantification of how severely within-cohort classifiers fail when applied to external populations, and by demonstrating that on average, neither method recovered cross-cohort transferability.

### 4.2 Why Batch Correction Failed: The 12:1 Disparity

The batch correction failure is consistent with the PERMANOVA result. Cohort identity explains 19.3% of Aitchison variance (R²=0.193 across all five cohorts; marginal R²=0.172 within the four-cohort two-variable model); disease status explains a marginal 1.4% (R²=0.014) in the same model. This approximately 12-fold disparity, combined with near-perfect cohort predictability from CLR features (macro-AUC=0.9917), indicates that compositional space is overwhelmingly organized by cohort-of-origin. A linear batch correction method that targets the cohort dimension of variance will tend to also remove disease variance if the two sources are not orthogonal—an assumption that the empirical worsening of LOCO AUC under correction suggests is violated here.

ComBat-seq's aggressive worsening of LOCO AUC (mean LogReg: 0.640 → 0.493) is consistent with this interpretation. ComBat-seq models each feature's batch effect separately using a negative binomial regression framework well-suited to RNA-seq data but potentially poorly calibrated for the extreme sparsity (75.4% zero entries) and compositional structure of 16S microbiome data. The within-cohort AUC preservation (and even improvement for Kazakhstan) after ComBat-seq correction shows that the correction did not remove all biological signal—it reshaped the compositional space in a way that made within-cohort structure more accessible to linear classifiers while making cross-cohort alignment worse.

MMUPHin's milder but still negative effect is similarly interpretable: linear mixed-effects models on relative abundances may not sufficiently align the full multivariate compositional space when the batch effect is substantially larger than the biological effect.

These results have direct implications for the microbiome harmonization literature. When the batch effect dominates by a factor of ~12, the two correction methods evaluated here could not reliably separate the two sources of variance. The field may need fundamentally different approaches, such as multi-cohort training with cohort as a stratification variable, domain adaptation methods, or meta-analytic frameworks that explicitly model between-cohort heterogeneity rather than attempting to remove it.

### 4.3 Implications of the Zhu 2022 Paired-Sample Design

The Zhu 2022 cohort deposited both blood microbiome (B_*) and fecal (F_*) samples from the same 90 participants (PRJNA489760). The blood samples were excluded from this analysis because blood and fecal microbiome communities represent fundamentally distinct compartments; including both would introduce a sample-type confound identical in structure to a batch effect, inflating within-cohort AUC through tissue-compartment discrimination rather than disease discrimination.

With fecal-only samples (n=60 binary), the near-perfect within-cohort AUC (logistic regression 0.998, LightGBM 0.979) is more difficult to attribute to a single confound. The small sample size (30 per class) is relevant: with 10-fold CV on 60 samples, each test fold contains approximately 6 samples, and single-fold AUC estimates can be highly variable. Zhu 2022 showed the largest detectable diagnosis-associated compositional effect (R²=0.094, F=6.01, p=0.0001) of the four labeled cohorts, but an AUC of 0.998 likely still reflects optimistic estimation at this sample size.

The LOCO result for Zhu 2022 (logistic regression AUC=0.813, the highest LOCO value among all cohorts) is informative: classifiers trained on the other three cohorts transfer to Zhu 2022 fecal samples better than to any other cohort. This suggests that the fecal disease signal in this cohort is somewhat consistent with signals in the Chinese and Kazakh training cohorts, unlike the spurious compartment-discrimination signal that would have characterized the original mixed dataset. We recommend that future AD microbiome studies report sample types explicitly and analyze each compartment independently before pooling.

### 4.4 Directional Flip Taxa: Empirical Confirmation of Cross-Population Contradictions

The directional flip taxa identified under OOF SHAP—twelve for logistic regression and twelve for LightGBM—confirm empirically what the meta-analytic literature has suggested qualitatively. The prevalence of directional flips across both model types is a stronger finding than previously appreciated: the original in-sample SHAP analysis had identified ten flips for logistic regression and three for LightGBM, but the OOF analysis reveals that directional inconsistency is comparably severe across both model families.

*Akkermansia* shows a complex cross-cohort pattern. For logistic regression, it is CN-associated in Zhuang 2018 and Ling 2020 but AD-associated in Zhu 2022 and Kazakhstan. For LightGBM, it is the dominant predictor in Ling 2020 (CN-associated) but shows mixed directionality across other cohorts. *A. muciniphila* abundance is strongly influenced by host diet, particularly fiber and polyphenol intake (Plovier et al., 2017), and geographic-population-specific dietary patterns likely modulate its abundance independently of any AD-specific effect. These dietary explanations are speculative in the absence of co-variate dietary data.

*Lactobacillus* shows a more consistent directional pattern: CN-associated in Chinese cohorts but AD-associated in Kazakhstan across both model types. Fermented dairy consumption (kumis, kefir) is culturally prominent in Central Asian populations and would be expected to substantially elevate *Lactobacillus* abundance in Kazakh controls. This speculative dietary hypothesis, if correct, would imply that the same relative abundance of *Lactobacillus* carries different disease-predictive information in Chinese versus Kazakh populations—a biological generalization problem not correctable by standard batch correction.

The critical implication of directional flips is that they may reflect biological population heterogeneity, technical study effects, or interactions between the two; the present design cannot disentangle these sources. Standard batch correction methods, designed to remove technical batch effects, would not be expected to resolve directional disagreements to the extent that these differences are driven by population-level biological differences rather than technical artifacts alone.

### 4.5 Implications for Clinical Translation

The results of this study suggest that gut microbiome-based AD classifiers are not ready for cross-population clinical application in their current form. A classifier trained on a Chinese cohort and applied to a Kazakh patient population would perform at or below chance for logistic regression, and only modestly above chance even for LightGBM. These performance levels are clinically inadequate for any diagnostic application.

We propose three prerequisites for meaningful clinical translation of microbiome-based AD diagnostics: (1) **Multi-cohort harmonized training data**: classifiers should be trained on samples from multiple geographic and demographic populations simultaneously, with cohort or population as a stratification variable, rather than training within a single cohort and expecting transfer. (2) **Prospective cross-cohort validation as a publication standard**: any future study reporting a gut microbiome AD classifier should be required to validate performance on at least one independent cohort from a different geographic population before claiming clinical relevance. The within-cohort AUC, no matter how high, is an insufficient evidential standard. (3) **Covariate-adjusted models**: diet, geographic origin, antibiotic history, and other known microbiome-modulating covariates should be incorporated into classification models either as features or confounders, to reduce the influence of population-specific compositional baselines.

The multi-cohort harmonized collection approach has succeeded in other biomarker domains (e.g., the UK Biobank for genomic predictors, the ADNI consortium for neuroimaging). A parallel initiative—perhaps a consortium of AD microbiome studies committed to harmonized collection protocols and prospective cross-cohort sharing—would be the most direct route to a clinically actionable microbiome diagnostic.

### 4.6 Limitations

Several limitations of this study should be acknowledged explicitly.

**Cohort search and selection:** The search strategy documented in Section 2.1 and Supplementary Note S1 identified five publicly available cohorts meeting inclusion criteria. The cohort set is not exhaustive of all existing AD microbiome studies, and the inclusion of additional cohorts from non-Asian populations would strengthen the conclusions. The comparability of diagnostic criteria across cohorts (clinical AD diagnosis vs. amyloid-PET positivity) adds biological heterogeneity that cannot be fully disentangled from technical batch effects.

**Mixed single-end and paired-end sequencing:** Kazakhstan used single-end NovaSeq reads while all other cohorts used paired-end MiSeq 2×300 bp. Single-end reads yield lower taxonomic resolution at equivalent read length due to the absence of overlap-based error correction in DADA2, and different error models were fitted per cohort. This technical asymmetry contributes to between-cohort heterogeneity and cannot be disentangled from biological population differences in this dataset.

**Kim/KBASE 2022 excluded from supervised analyses:** Per-sample amyloid-PET diagnosis labels are IRB-restricted and not publicly deposited. The Kim 2022 cohort (78 participants; 18 amyloid-positive preclinical AD, 60 CN) therefore contributed only to PERMANOVA variance decomposition.

**China overrepresentation:** Three of the four labeled cohorts are Chinese. The cross-cohort experiment is therefore primarily a China-internal transfer (Zhuang ↔ Ling ↔ Zhu 2022) plus transfer to/from a single Central Asian cohort. A more geographically balanced design—incorporating North American, European, and African populations—would be needed to generalize the conclusions to the full range of human diversity.

**Partial SRA deposit for Zhu 2022:** The public archive contains 180 samples representing 90 participants, a subset of the 302 participants analyzed in the published study (PRJNA489760). The deposited subset may not be a random sample of the full study population, potentially introducing selection bias of unknown direction.

**Batch correction leakage in LOCO evaluation:** As noted in Section 2.7, the batch correction analysis is transductive: test-cohort data and labels were included during correction parameter estimation. These corrected LOCO results reflect a transductive design and should not be extrapolated to prospective deployment.

**Cross-sectional design and observational data:** All cohorts used cross-sectional sampling. Causal inference about microbiome–AD relationships is not supported by this design. All language throughout this paper refers to associations and predictive performance, not causal effects.

**SHAP interpretability limitations:** SHAP directions reflect model-derived feature contributions, not direct estimates of biological association. For logistic regression, directions are conditional on regularization and collinearity among CLR features; for LightGBM, SHAP direction may reflect nonlinear interactions. The directional flip classifications and dietary hypotheses should be treated as empirical patterns worthy of targeted follow-up, not as established mechanistic findings.

---

## 5. Conclusion

We conducted a systematic cross-cohort generalization evaluation of gut microbiome-based Alzheimer's disease classifiers, spanning five independent cohorts from three countries. Within-cohort classifiers achieved AUC-ROC values of 0.633–0.998, consistent with published reports; however, LOCO cross-cohort AUC values fell to 0.504–0.813, a mean drop of 0.176 AUC points for logistic regression. Standard batch correction methods (ComBat-seq and MMUPHin) reduced mean LOCO performance and did not recover cross-cohort transferability. PERMANOVA revealed that cohort identity explained approximately 12 times more compositional variance than disease status (marginal R²: 0.172 vs. 0.014), and a cohort identity classifier achieved macro-AUC=0.9917, confirming that cohort structure dominates the feature space. Sensitivity analyses confirmed that generalization failure persists across all three prespecified cohort-exclusion configurations. Out-of-fold SHAP analysis identified near-zero taxon overlap between cohort-specific predictive signatures (mean Jaccard 0.135 for logistic regression, above the random baseline; 0.058 for LightGBM, within the random baseline) and twelve directional flip taxa for each model type—genera whose AD or CN association reverses across populations. These findings demonstrate that the gut microbiome-based AD classifier literature has a fundamental generalization problem that was not recovered by the two post-hoc batch-correction methods evaluated here. Multi-cohort harmonized prospective studies and mandatory cross-cohort validation are needed before microbiome-based AD biomarkers can achieve clinical utility.

---

## Data Availability Statement

All raw sequencing data analyzed in this study are publicly available from NCBI SRA or ENA under accession numbers PRJNA554111 (Zhuang 2018), PRJNA633959 (Ling 2020), PRJNA489760 (Zhu 2022), PRJNA811324 (Kaiyrlykyzy 2022), and PRJEB50447 (Kim 2022). Analysis code, processed data tables, and figure generation scripts are available at https://github.com/carhanc/microbiome-ad-generalization.

---

## Author Contributions

A.C. conceived the study, developed and executed the complete analysis pipeline, interpreted the results, and wrote the manuscript.

---

## Funding

The author declares that no external funding was received for the conduct of this research.

---

## Acknowledgements

The author thanks the original study teams of all five cohorts for depositing their data in publicly accessible archives, without which this analysis would not have been possible.

---

## Conflict of Interest

The author declares no conflict of interest.

---

## References

Callahan, B.J., McMurdie, P.J., Rosen, M.J., Han, A.W., Johnson, A.J.A., Holmes, S.P. (2016). DADA2: High-resolution sample inference from Illumina amplicon data. *Nature Methods*, 13(7), 581–583. https://doi.org/10.1038/nmeth.3869

Cryan, J.F., O'Riordan, K.J., Cowan, C.S.M., Sandhu, K.V., Bastiaanssen, T.F.S., Boehme, M., et al. (2019). The microbiota–gut–brain axis. *Physiological Reviews*, 99(4), 1877–2013. https://doi.org/10.1152/physrev.00018.2018

Harach, T., Marungruang, N., Duthilleul, N., Cheatham, V., Mc Coy, K.D., Frisoni, G., et al. (2017). Reduction of Abeta amyloid pathology in APPPS1 transgenic mice in the absence of gut microbiota. *Scientific Reports*, 7(1), 41802. https://doi.org/10.1038/srep41802

Jemimah, S., Chabib, C.M.M., Hadjileontiadis, L., AlShehhi, A. (2023). Gut microbiome dysbiosis in Alzheimer's disease and mild cognitive impairment: a systematic review and meta-analysis. *PLOS ONE*, 18(5), e0285346. https://doi.org/10.1371/journal.pone.0285346

Jung, J.H., Kim, G., Byun, M.S., Lee, J.H., Yi, D., Park, H., Lee, D.Y., KBASE Research Group. (2022). Gut microbiome alterations in preclinical Alzheimer's disease. *PLOS ONE*, 17(12), e0278276. https://doi.org/10.1371/journal.pone.0278276

Kaiyrlykyzy, A., Kozhakhmetov, S., Babenko, D., Zholdasbekova, G., Alzhanova, D., Olzhayev, F., Baibulatova, A., Kushugulova, A.R., Askarova, S. (2022). Study of gut microbiota alterations in Alzheimer's dementia patients from Kazakhstan. *Scientific Reports*, 12, 15115. https://doi.org/10.1038/s41598-022-19393-0

Ke, G., Meng, Q., Finley, T., Wang, T., Chen, W., Ma, W., et al. (2017). LightGBM: A highly efficient gradient boosting decision tree. In *Advances in Neural Information Processing Systems* 30, pp. 3146–3154.

Kowalski, K., Mulak, A. (2019). Brain-gut-microbiota axis in Alzheimer's disease. *Journal of Neurogastroenterology and Motility*, 25(1), 48–60. https://doi.org/10.5056/jnm18087

Ling, Z., Zhu, M., Liu, X., Shao, L., Cheng, Y., Yan, X., et al. (2021). Structural and functional dysbiosis of fecal microbiota in Chinese patients with Alzheimer's disease. *Frontiers in Cell and Developmental Biology*, 8, 634069. https://doi.org/10.3389/fcell.2020.634069

Liu, P., Wu, L., Peng, G., Han, Y., Tang, R., Ge, J., et al. (2019). Altered microbiomes distinguish Alzheimer's disease from amnestic mild cognitive impairment and health in a Chinese cohort. *Brain, Behavior, and Immunity*, 80, 633–643. https://doi.org/10.1016/j.bbi.2019.05.008

Liu, S., Gao, J., Zhu, M., Liu, K., Zhang, H.L. (2020). Gut microbiota and dysbiosis in Alzheimer's disease: implications for pathogenesis and treatment. *Molecular Neurobiology*, 57(12), 5026–5043. https://doi.org/10.1007/s12035-020-02073-3

Lundberg, S.M., Lee, S.I. (2017). A unified approach to interpreting model predictions. In *Advances in Neural Information Processing Systems* 30, pp. 4765–4774.

Ma, S., Shungin, D., Mallick, H., Schirmer, M., Nguyen, L.H., Kolde, R., Franzosa, E., Vlamakis, H., Xavier, R., Huttenhower, C. (2022). Population structure discovery in meta-analyzed microbial communities and inflammatory bowel disease using MMUPHin. *Genome Biology*, 23(1), 208. https://doi.org/10.1186/s13059-022-02753-4

Martín-Fernández, J.A., Barceló-Vidal, C., Pawlowsky-Glahn, V. (2003). Dealing with zeros and missing values in compositional data sets using nonparametric imputation. *Mathematical Geology*, 35(3), 253–278. https://doi.org/10.1023/A:1023866030544

Oksanen, J., Simpson, G.L., Blanchet, F.G., Kindt, R., Legendre, P., Minchin, P.R., et al. (2024). *vegan: Community Ecology Package*. R package version 2.6-8. https://CRAN.R-project.org/package=vegan

Pedregosa, F., Varoquaux, G., Gramfort, A., Michel, V., Thirion, B., Grisel, O., et al. (2011). Scikit-learn: Machine learning in Python. *Journal of Machine Learning Research*, 12, 2825–2830.

Plovier, H., Everard, A., Druart, C., Depommier, C., Van Hul, M., Geurts, L., et al. (2017). A purified membrane protein from *Akkermansia muciniphila* or the pasteurized bacterium improves metabolism in obese and diabetic mice. *Nature Medicine*, 23(1), 107–113. https://doi.org/10.1038/nm.4236

Vogt, N.M., Kerby, R.L., Dill-McFarland, K.A., Harding, S.J., Merluzzi, A.P., Johnson, S.C., et al. (2017). Gut microbiome alterations in Alzheimer's disease. *Scientific Reports*, 7(1), 13537. https://doi.org/10.1038/s41598-017-13601-y

Wang, J., Liu, H., Lai, H., Bao, Y., Shen, M., Li, C., Ma, L., Wu, T., Yang, S., Du, X., O'Brien, T.J., Zhang, J., Zhang, L. (2026). Identifying candidate gut microbiota indicators for Alzheimer's disease through integrated data. *mSystems*, 11(4). https://doi.org/10.1128/msystems.01754-25

World Health Organization. (2023). *Dementia*. WHO Fact Sheet. https://www.who.int/news-room/fact-sheets/detail/dementia

Zhang, Y., Parmigiani, G., Johnson, W.E. (2020). ComBat-seq: Batch effect adjustment for RNA-seq count data. *NAR Genomics and Bioinformatics*, 2(3), lqaa078. https://doi.org/10.1093/nargab/lqaa078

Zhu, Z., Ma, X., Wu, J., Xiao, Z., Wu, W., Ding, S., Zheng, L., Liang, X., Luo, J., Ding, D., Zhao, Q. (2022). Altered gut microbiota and its clinical relevance in mild cognitive impairment and Alzheimer's disease: Shanghai Aging Study and Shanghai Memory Study. *Nutrients*, 14(19), 3959. https://doi.org/10.3390/nu14193959

Zhuang, Z.Q., Shen, L.L., Li, W.W., Fu, X., Zeng, F., Gui, L., et al. (2018). Gut microbiota is altered in patients with Alzheimer's disease. *Journal of Alzheimer's Disease*, 63(4), 1337–1346. https://doi.org/10.3233/JAD-180176

---
