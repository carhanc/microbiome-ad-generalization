# Response to Reviewers

**Manuscript 1918251**
*Cross-Cohort Generalization Failure in Gut Microbiome-Based Alzheimer's Disease Classifiers: Evidence from Five Independent Cohorts with Batch Correction Analysis*
*Frontiers in Aging Neuroscience*

---

We thank the editor and both reviewers for their careful, constructive engagement with this manuscript across two rounds of review. Every substantive methodological and statistical concern raised — by Reviewer 3 in the first round and by Reviewer 4 in the independent report that followed — has been addressed with a concrete analysis, code change, or manuscript revision, described point by point below with the exact new numerical result and its location in the current manuscript. Two of these checks changed a previously reported conclusion, and we report both changes directly rather than preserving the earlier result: the LightGBM cross-cohort taxonomic-overlap finding is no longer significant under a more realistic null model, and the ~12-fold cohort-versus-diagnosis PERMANOVA disparity is now shown to be highly sensitive to cohort composition (falling to ~2.6-fold in a more technically comparable subset). Where a point raised by Reviewer 4 had already been addressed during the revision prepared in response to Reviewer 3, we say so explicitly and describe what has been added specifically in response to Reviewer 4's independent report.

---

## Section 1: Reviewer 3 — Round 2 Follow-Up Comments

Reviewer 3 raised five methodological concerns in their second-round report. We note that Reviewer 3 subsequently withdrew from the review process; nevertheless, all five comments identified genuine gaps, and we addressed every one of them in full, as described below.

### 1.1 Global feature selection used held-out-cohort information

**Comment:** The genus/feature set used in supervised cross-cohort analyses was defined using all five cohorts, including the eventual held-out cohort, rather than being re-derived from training data alone.

**Response:** Agreed. We built a new training-only feature-selection sensitivity analysis (`scripts/09_training_only_feature_sensitivity.py`) that re-derives the ≥20%-prevalence genus filter independently within every outer fold (within-cohort), every LOCO hold-out, and every pairwise train cohort, using only that split's training-side samples. LOCO AUC under this strict design changed by at most 0.066 in either direction (largest: Kazakhstan LightGBM Δ=−0.066); no cohort–model combination changed its qualitative classification. **Location:** Methods 2.2, 2.5.1; Results 3.3; Limitations; Supplementary Table S4.

### 1.2 Batch correction used held-out diagnosis labels

**Comment:** The batch-correction analysis (ComBat-seq/MMUPHin) used held-out diagnosis labels during correction, which cannot estimate prospective external-validation performance, and corrected within-cohort AUC cannot establish that biological signal was preserved.

**Response:** Agreed with both halves. We built a label-blind transductive correction design (`scripts/10_batch_correction_label_blind.R`, `scripts/11_corrected_generalization_label_blind.py`) that uses cohort membership and feature values but never a sample's diagnosis label, and we now present this as the manuscript's primary batch-correction result. Under label-blind correction, mean LOCO AUC changed by no more than 0.007 relative to uncorrected for either method or model; no cohort collapsed. The originally reported severe degradation (e.g., Kazakhstan ComBat-seq LOCO AUC=0.308) occurred only under the label-informed design and is retained solely as an explicitly labeled exploratory sensitivity analysis, not as evidence about batch correction's effect on biological signal in general. **Location:** Methods 2.7; Results 3.5 (retitled "Label-Blind Transductive Batch Adjustment Does Not Improve Cross-Cohort Transferability"); Discussion 4.2 (retitled "Batch-Adjustment Results and Limits of Mechanistic Interpretation"); Limitations; Supplementary Table S5; Supplementary Figure S5.

### 1.3 Mean signed SHAP is not a valid direction statistic

**Comment:** SHAP direction was defined from mean signed SHAP, which for a linear explainer under interventional feature perturbation depends on both the fitted coefficient and the background-vs-evaluation abundance difference, not coefficient sign alone.

**Response:** Agreed. Direction for logistic regression is now defined from the fitted model's own coefficient, using a pre-specified descriptive stability screen (stable sign in ≥8 of 10 outer folds), not from mean signed SHAP anywhere in this paper. This screen identifies eight genera with opposite, stable coefficient directions across cohorts (down from the invalid statistic's count of twelve), and while implementing this fix we found a concrete illustration of why it mattered: *Akkermansia*'s fitted coefficient is stable and positive in all four cohorts even though its pooled mean SHAP value appears to reverse sign — it is not a flip taxon under the correct definition. No comparable direction statistic is computed for LightGBM, which has no single coefficient-like global direction; the previously reported LightGBM directional-flip count is not used anywhere in the current manuscript. **Location:** Methods 2.8; Results 3.6; Discussion 4.4 (retitled "Coefficient-Sign Stability Screen: Descriptive Cross-Cohort Direction Heterogeneity"); Limitations; Figure 6; Supplementary Figure S2; Supplementary Table S6.

### 1.4 Jaccard null did not match the observed statistic's construction

**Comment:** The reported overlap statistic is the mean of six pairwise Jaccard values among four cohorts' top-N sets, but the null distribution it was compared against drew only two random sets and computed a single pairwise Jaccard per replicate.

**Response:** Agreed; confirmed exactly as described. The null now draws four independent random top-N sets per replicate (100,000 replicates) and averages the same six pairwise Jaccards among them, matching the observed statistic's construction exactly. Under this corrected null, logistic regression's overlap remains significant at all tested N, and LightGBM's overlap — previously reported as "within the random baseline" (p=0.068 under the old, mismatched null) — is significant at N=20 (p=0.0048) and N=50, though not at N=10. We report this as a real but N-dependent and comparatively small signal rather than either the earlier "no shared features" framing or an overstated "confirmed overlap" framing. **Location:** Methods 2.8; Results 3.6; Abstract; Conclusion; Supplementary Table S7.

### 1.5 PERMANOVA does not by itself explain the batch-correction result

**Comment:** PERMANOVA demonstrates cohort-associated compositional heterogeneity but does not by itself establish that this is why batch correction changed performance in a particular direction; empirical result and causal hypothesis must be kept separate.

**Response:** Agreed. We restructured the relevant Discussion passage into an explicit "Empirical result" / "Hypothesis, clearly flagged as such" structure, and every place the cohort-versus-diagnosis PERMANOVA disparity is invoked in connection with batch correction now states explicitly that it provides context but does not establish that the two sources of variance overlap or explain the specific pattern observed. **Location:** Results 3.4; Discussion 4.2.

---

## Section 2: Reviewer 4 — Independent Report

Reviewer 4's central assessment — that the main empirical finding of reduced cross-cohort performance is supported and robust — is one we share. Reviewer 4 evaluated the version of the manuscript available at the time of their review; where a point below overlaps with a change already made during the intervening revision prepared in response to Reviewer 3, we say so and describe what has been added specifically in response to Reviewer 4's report.

### 2.1 Batch correction: primary analysis and framing

**Request:** Clarify that the primary batch-correction analysis is label-blind and transductive, that the label-informed analysis is exploratory only, and that neither is described as prospective.

**Response:** During the intervening revision, we had already implemented the label-blind design as the primary analysis and the label-informed design as an explicitly labeled exploratory sensitivity (Section 1.2 above). In response to Reviewer 4, we additionally: retitled Results 3.5 and Discussion 4.2 to state this scope directly in the section headings; removed language suggesting a specific causal mechanism for the label-informed collapse ("most likely because...") in favor of hedged hypothesis language; and removed the term "oracle-style" from the published manuscript, replacing it with "label-informed transductive sensitivity analysis" throughout. **Location:** Methods 2.7; Results 3.5; Discussion 4.2; Limitations; Supplementary Figure S5 caption.

### 2.2 Cohort as a composite variable; no technical/biological decomposition claim

**Request:** State explicitly that "cohort" bundles technical and biological differences that PERMANOVA cannot decompose, and remove language implying PERMANOVA separates technical batch effects from biological heterogeneity.

**Response:** During the intervening revision this framing was still present in the Introduction ("the relative contribution of technical batch effects versus biological heterogeneity... as estimated by variance partitioning"). In response to Reviewer 4, we rewrote this passage to state that cohort is a composite, study-of-origin variable that PERMANOVA cannot decompose into technical and biological sources, and added the same explicit statement to Methods 2.6 and to every place the ~12-fold cohort-versus-diagnosis disparity is reported (Abstract, Results 3.4, Discussion 4.2, Limitations, Conclusion). **Location:** Introduction; Methods 2.6; Results 3.4; Discussion 4.2; Limitations; Abstract; Conclusion.

### 2.3 PERMANOVA/dispersion sensitivity

**Request:** Check which cohort pairs drive the significant beta-dispersion difference, and test the cohort-versus-diagnosis variance disparity on a more technically comparable subset of cohorts.

**Response:** New (`scripts/14_permanova_dispersion_sensitivity.R`). Pairwise post-hoc beta-dispersion comparisons (Holm/BH-adjusted) show the overall significant dispersion difference is driven predominantly by Ling 2020, which differs significantly from every other cohort (all Holm-adjusted p≤0.0063); Zhu 2022 vs. Zhuang 2018 and Kazakhstan vs. Kim/KBASE are not significantly different from one another. Restricted to the three more technically comparable Chinese paired-end MiSeq cohorts (Zhuang 2018, Ling 2020, Zhu 2022 — still not technically identical), marginal PERMANOVA gives cohort R²=0.0588 (F=10.05, p=0.0001) versus diagnosis R²=0.0229 (F=7.83, p=0.0001) — **an approximately 2.6-fold disparity, not the ~12-fold figure from the full four-cohort model.** We report this directly as evidence that the ~12-fold figure is not a fixed, generalizable property of "cohort effects" and depends substantially on which cohorts are compared, and we no longer display it as a standalone annotation on Figure 4. **Location:** Methods 2.6; Results 3.4; Discussion 4.2; Limitations; Abstract; Conclusion; Figure 4 (annotation removed); Supplementary Table S11, S12.

### 2.4 Zhu 2022 robustness

**Request:** Assess whether the near-perfect Zhu 2022 within-cohort AUC (n=60) is robust to fold assignment and formally distinguishable from a chance-label null.

**Response:** New (`scripts/12_zhu_robustness.py`), both run under the strict training-only pipeline. **Repeated nested cross-validation** (50 repetitions, different pre-specified seed each time): mean pooled-OOF AUC 0.9955 (SD 0.0027) for logistic regression and 0.9760 (SD 0.0103) for LightGBM — a fold-assignment stability check, not independent replication, since the same 60 participants are used in every repetition. **Label permutation test** (exactly 1,000 permutations, same strict pipeline, fixed partition): null AUC centered near chance for both models (mean ≈0.47); no permutation among the 1,000 reached the observed AUC for either model, giving empirical p = 1/1001 ≈ 0.000999 for both models — the resolution floor at this permutation count. **Location:** Methods 2.4.1; Results 3.2; Discussion 4.3; Abstract; Conclusion; Supplementary Table S8, S9; Supplementary Figure S6.

### 2.5 Preprocessing leakage: explicit classification

**Request:** Explicitly classify the leakage status of every preprocessing step rather than leaving it to be inferred.

**Response:** New text added to Methods 2.2 explicitly classifying every preprocessing operation: DADA2 (cohort-specific, label-free, no leakage pathway); global genus selection (uses all cohorts, transparently identified as not leakage-free, bounded by the training-only sensitivity analysis in 2.5.1); zero-replacement/CLR (deterministic, sample-local, no leakage pathway — and, per a related wording fix, applied independently to each sample, not "jointly" across samples); and model tuning/fitting (training-side samples only). No new preprocessing was run; the existing training-only sensitivity analysis already quantifies the one step that is not leakage-free. **Location:** Methods 2.2.

### 2.6 Formal within-cohort-vs-LOCO AUC comparison

**Request:** Replace the non-overlapping-CI heuristic with a formal statistical comparison of within-cohort and LOCO AUC.

**Response:** New (`scripts/13_auc_difference_bootstrap.py`). A paired, diagnosis-stratified, participant-level bootstrap (10,000 replicates) uses identical resampled participant indices to compute both the within-cohort and LOCO AUC in each replicate, so the resulting delta distribution reflects a genuine paired comparison. Five of eight cohort–model pairs (Ling 2020 both models, Zhu 2022 both models, Kazakhstan logistic regression) have a 95% CI on the delta that excludes zero; the remaining three (Zhuang 2018 both models, Kazakhstan LightGBM) show the same-signed delta but a CI that includes zero, and are not formally significant at this sample size. We report this distinction directly rather than treating all eight LOCO drops as equally well-supported, and we removed the "non-overlap of CIs indicates a statistically reliable drop" language from Methods, replacing it with this formal comparison. Every AUC 95% CI in the paper (Table 2, Table 3, Table S1) was regenerated using this same bootstrap procedure for consistency. **Location:** Methods 2.5.2; Results 3.3; Table 3; Abstract; Discussion 4.1; Conclusion; Limitations; Supplementary Table S10.

### 2.7 SHAP direction terminology

**Request:** Use "positive/negative fitted coefficient" language rather than "AD-associated"/"CN-associated"/"biological association" framing for a fitted-model coefficient pattern that has not been independently validated or covariate-adjusted.

**Response:** During the intervening revision, the move to coefficient-based direction (Section 1.3 above) had already eliminated mean-signed-SHAP-based direction claims. In response to Reviewer 4, we did a full-text sweep for "AD-associated," "CN-associated," "association with AD status," live "biological association" claims, and "directionally consistent," and rewrote every remaining instance to use only positive/negative fitted-coefficient language. We also substantially shortened the Akkermansia/Lactobacillus dietary-speculation passages, retaining only the fitted-coefficient pattern, a brief statement that diet is a plausible unmeasured confounder, and an explicit statement that no dietary covariate data exist and no mechanistic attribution is made. **Location:** Results 3.6; Discussion 4.4 (retitled).

### 2.8 Jaccard null model realism

**Request:** Check whether the Jaccard-overlap null is realistic given that genus prevalence varies substantially by cohort, and do not assume previously reported significance survives a more realistic null.

**Response:** New (`scripts/15_jaccard_prevalence_null.py`), distinct from the four-set/six-pair statistic-matching fix made in response to Reviewer 3 (Section 1.4 above). Each cohort's top-N ranking and random draws are now additionally restricted to genera reaching ≥20% prevalence within that cohort's own AD/CN samples (a per-cohort eligible universe), rather than the shared 396-genus pool. Logistic regression remains significant at all three tested N under this more realistic null. **LightGBM's overlap, significant at N=20/50 under the shared-pool null, is not significant at any of N=10/20/50 under the prevalence-restricted null** (p=0.80, 0.24, 0.061). We report this as a genuine reversal: LightGBM's apparent cross-cohort top-predictor overlap does not survive this more rigorous check. **Location:** Methods 2.8; Results 3.6; Abstract; Limitations; Figure 6 caption; Supplementary Table S13.

### 2.9 Geographic and population scope

**Request:** Scope claims to the cohorts/datasets actually evaluated rather than implying a universal population-level or cross-population biological effect.

**Response:** We replaced "classifiers trained within one cohort consistently fail to generalize to other cohorts" with "classifiers showed substantial and variable degradation when transferred to independent cohorts" in the Introduction; removed "the cohort effect is structurally much larger than the disease effect" from the Abstract; removed the unrestricted claim that classifiers are "not ready for cross-population clinical application," replacing it with "not ready for clinical deployment across heterogeneous external cohorts without prospective validation"; and retitled Discussion 4.4 to use "cross-cohort" rather than "cross-population" language for this study's own findings (prior-literature discussion of geographically distinct populations is retained where it refers to the cited literature). The Conclusion is scoped to "substantial cross-cohort generalization limitations among the evaluated publicly available 16S datasets," avoiding "fundamental generalization problem," "universal failure," or "all populations" language, and uses "prospective cross-cohort validation" rather than "mandatory cross-cohort validation." **Location:** Abstract; Introduction; Discussion 4.4, 4.5; Conclusion.

### 2.10 Speculative dietary/mechanistic language

**Request:** Shorten and explicitly caveat dietary/mechanistic speculation about individual taxa (Akkermansia, Lactobacillus).

**Response:** Addressed together with 2.7 above: the Akkermansia and Lactobacillus discussions in both Results 3.6 and Discussion 4.4 were substantially shortened, keeping only the fitted-coefficient pattern, a brief unmeasured-confounder statement, and an explicit no-diet-data/no-mechanistic-attribution statement. **Location:** Results 3.6; Discussion 4.4.

### 2.11 Main conclusion wording

**Request:** State the main conclusion using language scoped to the cohorts and datasets actually evaluated.

**Response:** The Conclusion states: "These findings demonstrate substantial cross-cohort generalization limitations among the evaluated publicly available 16S datasets." This wording is used consistently across the Abstract, Introduction, and Conclusion. **Location:** Conclusion; Abstract; Introduction.

---

## Summary of findings that changed as a result of this process

We want two changes flagged explicitly rather than left for a reader to notice on their own:

1. **LightGBM's cross-cohort top-predictor Jaccard overlap is no longer reported as statistically significant** under a more realistic, prevalence-restricted null model (Section 2.8). Logistic regression's overlap finding is unchanged and remains significant under both null designs.
2. **The ~12-fold cohort-vs-diagnosis PERMANOVA variance disparity is now explicitly shown to be sensitive to cohort composition** (falling to ~2.6-fold in a more technically comparable three-cohort subset) and is no longer presented as a fixed, generalizable figure or displayed as a standalone figure annotation (Section 2.3).

Additionally, **three of eight within-cohort-vs-LOCO AUC comparisons do not reach formal statistical significance** under the new paired bootstrap test, though all eight show the same directional pattern (Section 2.6). None of these changes alters the paper's central conclusion that cross-cohort generalization is substantially and measurably limited among the cohorts evaluated; the manuscript is now more conservative and more precisely scoped than the version originally reviewed.

We thank the editor and both reviewers again for the opportunity to substantially strengthen this manuscript's statistical rigor and precision.
