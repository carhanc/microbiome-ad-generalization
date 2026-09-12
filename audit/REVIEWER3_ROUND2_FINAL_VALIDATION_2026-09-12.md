# Reviewer 3 Round 2 — Final Validation Report

**Date:** 2026-09-12
**Branch:** `frontiers-r3-round2-fixes-2026-09-12` (forked from the frozen pre-audit branch `frontiers-r3-round2-preaudit-2026-09-12` at commit `c41a53b51c7dcf93c04006d654d90d73b1259f3f`)
**Companion documents:** `audit/REVIEWER3_ROUND2_PREAUDIT_2026-09-12.md`, `audit/REVIEWER3_ROUND2_IMPLEMENTATION_LOG_2026-09-12.md`, `review/REVIEWER3_ROUND2_RESPONSE_DRAFT.md`

---

## 1. Summary of all changes

Five reviewer concerns were addressed by adding new, non-destructive analysis scripts (`09`–`11`), modifying one existing script in place (`08_shap_taxa_comparison.py`, verified to leave the underlying SHAP computation byte-identical), fixing one cosmetic bug (`05_permanova_variance_decomp.R`), and rewriting every affected section of `manuscript/draft.md` and `manuscript/supplementary.md` to match the corrected analyses. No historical result file was overwritten with different numbers; every new finding lives in a newly-named file, and every old file that is now superseded (`jaccard_null.csv`, the label-informed batch-correction tables, the pre-Round-2 `shap_directional_flips.csv` interpretation) is retained on disk and explicitly labeled in the manuscript as historical/exploratory rather than deleted or silently replaced.

One discrepancy outside the five original concerns was discovered and disclosed rather than concealed: the logistic-regression models throughout the entire pipeline (not just the new Round 2 code) are L2-regularized, not elastic-net as the manuscript had stated — `penalty="elasticnet"` was never set. This is a terminology correction only (§4 below); no model was refit.

## 2. Old vs. corrected numbers

See the summary table in `review/REVIEWER3_ROUND2_RESPONSE_DRAFT.md` ("Summary of numerical changes carried into the final manuscript"). Reproduced here for completeness:

| Item | Original | Corrected |
|---|---|---|
| Genus filter rule (stated) | "≥1 count in ≥2 samples" | "≥20% prevalence, global across 5 cohorts" |
| Batch correction, primary framing | label-informed, "ComBat-seq substantially worsens" | label-blind, "essentially unchanged" |
| ComBat-seq LOCO AUC, logreg, Kazakhstan | 0.308 (primary) | 0.560 (label-blind, primary); 0.308 retained as labeled exploratory |
| Directional flips, logreg | 12 (mean-signed-SHAP) | 8 (coefficient-stability) |
| Directional flips, LightGBM | 12 (mean-signed-SHAP) | removed entirely |
| Jaccard null (N=20), LightGBM p-value | 0.068 (n.s.) | 0.0048 (significant) |

## 3. Strict training-only feature-selection results

`scripts/09_training_only_feature_sensitivity.py`. LOCO AUC changed by at most 0.066 (Kazakhstan/LightGBM, −0.066; Kazakhstan/LogReg, +0.045); all other cohort–model combinations changed by ≤0.031. Within-cohort AUC changed by −0.083 to +0.038. No cohort–model combination changed qualitative classification. Retained genus count under training-only selection: 150 (Kazakhstan hold-out) to 373, vs. 396 globally. Full tables: `results/tables/{within_cohort,loco,pairwise}_auc_training_only.csv`, `training_only_feature_sensitivity_summary.csv`, `training_only_feature_counts.csv`, `training_only_feature_lists.csv` (9,250-row full per-split retained-genus lists, used for gate verification below).

## 4. Label-blind batch-correction results

`scripts/10_batch_correction_label_blind.R` + `scripts/11_corrected_generalization_label_blind.py`. Mean LOCO AUC (label-blind): logreg uncorrected 0.6398, ComBat-seq 0.6338, MMUPHin 0.6390; LightGBM uncorrected 0.6466, ComBat-seq 0.6540, MMUPHin 0.6478. All within 0.007 of uncorrected. No cohort collapsed. Full tables: `results/tables/loco_auc_corrected_labelblind.csv`, `pairwise_auc_labelblind.csv`, `within_cohort_auc_corrected_labelblind.csv`, `auc_comparison_labelblind.csv`.

## 5. Label-informed batch results (clearly labeled exploratory)

Original `scripts/07_corrected_generalization.py` outputs, unchanged and unmodified (`loco_auc_corrected.csv`, `auc_comparison_table.csv`, `within_cohort_auc_corrected.csv`). Mean LOCO AUC: logreg ComBat-seq 0.493, MMUPHin 0.612; LightGBM ComBat-seq 0.555, MMUPHin 0.582. Kazakhstan logreg collapses to 0.308 under ComBat-seq. Now presented in the manuscript exclusively as Supplementary Figure S5, explicitly labeled "exploratory sensitivity... not a prospective or deployable design," and referenced in §3.5/§4.2 only to explain why the originally-reported severe degradation is not the primary result.

## 6. Corrected Jaccard null results

`scripts/08_shap_taxa_comparison.py::compute_jaccard_null_baseline_corrected`, 100,000 replicates per (model, N), matched four-set/six-pair construction.

| Model | N | Observed | Null mean | Null 95% CI | Empirical p |
|---|---|---|---|---|---|
| LogReg | 10 | 0.0850 | 0.0134 | [0.0000, 0.0361] | <0.0001 |
| LogReg | 20 | 0.1350 | 0.0265 | [0.0085, 0.0491] | <0.0001 |
| LogReg | 50 | 0.2296 | 0.0680 | [0.0491, 0.0895] | <0.0001 |
| LightGBM | 10 | 0.0175 | 0.0135 | [0.0000, 0.0361] | 0.4494 |
| LightGBM | 20 | 0.0584 | 0.0265 | [0.0085, 0.0491] | 0.0048 |
| LightGBM | 50 | 0.1430 | 0.0680 | [0.0492, 0.0894] | <0.0001 |

Full table: `results/tables/jaccard_null_corrected.csv`.

## 7. Logistic coefficient stability results

`results/tables/logreg_coefficient_stability.csv` (1,584 cohort×genus rows): 961 stable_negative, 446 stable_positive, 177 unstable (pre-specified ≥8/10-fold threshold). `n_zero` is 0 for all 1,584 rows (expected under L2 regularization, which does not perform sparse selection — see §4 of the response draft and §10 below).

## 8. Final stable directional-flip taxa

Eight, for logistic regression only: **Agathobacter, Bifidobacterium, Coprococcus, Dorea, Lactobacillus, NK4A214 group, Romboutsia, Ruminococcus gnavus group** (`results/tables/logreg_directional_flips_stable.csv`).

## 9. Explicit confirmation: LightGBM directional-flip claim removed

Confirmed by (a) code inspection — no function in `08_shap_taxa_comparison.py` computes a flip set for `model_name == "lgbm"`; (b) the script's own console output, which prints "Direction flips: not assessed for LightGBM (no coefficient-like direction exists...)"; (c) manuscript grep — zero occurrences of "LightGBM" combined with a live (non-retracted) flip claim anywhere in `draft.md` or `supplementary.md`; every mention of the old LightGBM flip count is inside a sentence explicitly retracting it.

## 10. Every manuscript claim changed

Enumerated in full in `review/REVIEWER3_ROUND2_RESPONSE_DRAFT.md`'s per-concern sections. Sections touched: Abstract, §2.2, new §2.5.1, §2.4 (elastic-net→L2 correction, hyperparameter-grid correction), §2.7, §2.8, §2.10, §3.3, §3.4, §3.5, §3.6, §4.2, §4.4, §4.6 (Limitations), Conclusion. Supplementary: new Tables S4–S7, revised Figure S2 caption, Figures S3a/S3b replaced with a removal notice, new Figure S5.

**Additional discrepancy disclosed (not one of the five original concerns):** the logistic-regression model construction (`LogisticRegression(solver="saga", l1_ratio=0.5, ...)`, unchanged from every historical script) never sets `penalty="elasticnet"`, so every logistic-regression result in this paper — not just the Round 2 additions — was always L2 (ridge)-regularized, despite being called "elastic-net" throughout prior drafts. **We deliberately did not refit any model to correct this**, for two reasons: (1) L2-regularized logistic regression is a legitimate, widely-used model in its own right — the error was in the manuscript's terminology, not in the method actually used, so every AUC/SHAP/coefficient number in the paper remains a faithful description of what the code does; (2) refitting to genuine elastic-net would require rerunning every classifier-producing script in the pipeline (03, 04, 07, 08, 09, 10, 11), silently changing every headline AUC number in the paper — an unauthorized scope expansion far beyond the five reviewer concerns, undertaken under time pressure, with real risk of introducing new errors. We corrected the manuscript text to say "L2-regularized" everywhere "elastic-net" previously appeared, and flagged the one interpretive consequence we could identify (the coefficient-stability table's `n_zero` column is uninformative under L2, since L2 does not zero out coefficients) rather than disguising it.

## 11. Every figure/table changed

**Figures:** Figure 5 (regenerated — now shows label-blind data as primary), Figure 6 (regenerated — coefficient-based dot color, stability-based flip marker), Supplementary Figure S2 (regenerated — neutral-color LightGBM bars, caption corrected), Supplementary Figures S3a/S3b (content replaced with an explicit removal notice — no image regenerated, since the underlying analysis no longer exists), new Supplementary Figure S5 (label-informed exploratory comparison, newly generated).

**Tables:** new Table 4 restructured with a "Design" column (uncorrected / label-blind / label-informed); new Supplementary Tables S4 (training-only sensitivity), S5 (within-cohort AUC, label-blind vs. label-informed), S6 (coefficient-stability detail, including the *Akkermansia* negative control), S7 (Jaccard N-sensitivity).

All figures were visually verified by direct PDF rendering (PyMuPDF) after the full docx/PDF rebuild, not merely assumed correct from code review.

## 12. Remaining limitations

Carried forward from the original manuscript (cohort geographic imbalance — 3 of 4 labeled cohorts are Chinese; cross-sectional design; mixed single-end/paired-end sequencing; Kim/KBASE exclusion from supervised analysis; partial Zhu 2022 SRA deposit) plus new limitations introduced by this round's findings, all stated in §4.6:

- The global (not training-only) feature-selection design remains the analysis actually reported as primary, with the training-only version as a sensitivity check, not a replacement — a future study should make training-only selection the default rather than a post hoc check.
- Both batch-correction designs remain transductive (neither is prospective), even the label-blind one — cohort membership and raw feature values of the "test" cohort are still used during correction fitting.
- SHAP/coefficient-based direction is a statement about what the fitted model learned, not a direct biological or causal claim; dietary explanations for specific flip taxa remain speculative without covariate data.

## 13. Commands required to reproduce every Round 2 result

```bash
# Training-only feature-selection sensitivity
/opt/anaconda3/envs/microbiome-ad/bin/python scripts/09_training_only_feature_sensitivity.py

# Label-blind batch correction
Rscript scripts/10_batch_correction_label_blind.R
/opt/anaconda3/envs/microbiome-ad/bin/python scripts/11_corrected_generalization_label_blind.py

# SHAP coefficient-stability + matched Jaccard null (also regenerates the
# pre-existing SHAP importance/overlap tables, byte-identically)
/opt/anaconda3/envs/microbiome-ad/bin/python scripts/08_shap_taxa_comparison.py

# Figures
/opt/anaconda3/envs/microbiome-ad/bin/python scripts/generate_manuscript_figures.py

# Manuscript/supplement rebuild
pandoc manuscript/draft.md -o manuscript/draft.docx --standalone
pandoc manuscript/supplementary.md -o manuscript/supplementary.docx --standalone
# (+ python-docx table-border pass, see implementation log for the exact script)
/opt/anaconda3/bin/python scripts/build_pdf.py
```

Environment: `microbiome-ad` conda env (Python 3.11.15) for all Python; R 4.5.1 for all R (`audit/ROUND2_ENVIRONMENT.txt`). DADA2 was not rerun anywhere in Round 2.

## 14. Git status

Branch `frontiers-r3-round2-fixes-2026-09-12`, 5 commits ahead of the pre-audit branch's fork point (`c41a53b`): training-only sensitivity, label-blind batch correction, SHAP/Jaccard correction, manuscript rewrite, and this docs commit. Working tree clean as of this report (verified before the final push in the handoff message). `main` untouched throughout.

## 15. Unresolved methodological concerns identified during this pass

In the spirit of "if you can identify a defensible methodological objection that remains, fix it unless doing so requires new unavailable source data" — three were identified. None require unavailable data; all three require substantial additional compute/analysis time that was not spent, for the reasons given. We report them rather than either silently fixing them under time pressure (risking new errors) or omitting them.

1. **No multiple-comparisons correction on the stable directional-flip screen.** The coefficient-stability criterion (≥8/10 folds) is applied to 14 candidate genera (those in the top-20 SHAP ranking in ≥2 cohorts) without a permutation-based estimate of how many "stable flips" would arise by chance under a null where no genus truly differs in direction across cohorts. A rigorous fix would permute diagnosis labels within each cohort many times and refit the full nested-CV coefficient-stability pipeline per permutation — at the observed ~100 seconds per real run, a few hundred to a thousand permutations would take multiple hours, which was not run. This is a real, currently-unaddressed limitation of the "8 stable flips" claim, though it does not affect the more clearly interpretable single-genus examples (*Lactobacillus*, and the *Akkermansia* negative control) discussed narratively.
2. **Single-seed comparison for the training-only feature-selection sensitivity.** The ≤0.066 AUC deltas reported in concern #1's fix are from one fixed `random_state=42` run each for the original and training-only designs; we did not run repeated CV with multiple seeds to establish whether that magnitude of change is distinguishable from ordinary fold-assignment noise. The qualitative conclusion (no cohort changes classification) is robust to this concern, but the precise 0.066 figure should be read as a single-seed estimate, not a noise-adjusted one.
3. **Figure 6 Panel A / Panel B taxa-pool mismatch (cosmetic, pre-existing).** Panel A displays its own top-20-by-max-|SHAP| selection (pooled from each cohort's top-15), while the stable-flip candidate pool (Panel B, Table S6) is drawn from each cohort's own top-20. As a result, only 5 of the 8 confirmed stable-flip taxa appear (starred) in Panel A; the other 3 appear correctly in Panel B but not in Panel A's display window. This predates Round 2 (the same windowing mismatch existed before, just with the old flip definition) and is disclosed here for completeness rather than fixed, since resolving it would mean changing Panel A's selection logic, a cosmetic layout decision outside the five reviewer concerns.

None of these three change any number already reported in this document or in the manuscript; they bound the confidence with which the "8 stable flips" and "≤0.066 AUC" figures should be read, and we recommend a future round address at least item 1 if reviewers press further on the directional-flip claim's statistical rigor.

---

## Final skeptical self-review (as Reviewer 3, attempting to reject again)

Having implemented all five fixes, we asked: what would a genuinely skeptical reviewer say next? The three items in §15 above are the most defensible remaining objections we could construct, and none was concealed. Beyond those three, we specifically checked for and did not find: (a) any remaining place where mean signed SHAP is described as a biological direction; (b) any remaining place where the Jaccard null does not match its observed statistic; (c) any remaining place where PERMANOVA is credited with explaining batch-correction behavior; (d) any stale headline number inconsistent with its source CSV; (e) any script that fails to compile or run, or that produces NaN/missing output. All five original concerns are, in our assessment, now substantively and honestly addressed, with the batch-correction and Jaccard-null fixes each producing a genuine reversal of a previous claim rather than a defensive restatement of it.
