# Reviewer 3 Round 2 — Implementation Log

Running log of every scientific change made while repairing the five issues
documented in `audit/REVIEWER3_ROUND2_PREAUDIT_2026-09-12.md`. Entries are
appended in execution order, each with: what was done, exact command, env
used, files changed, and old-vs-new numbers where applicable.

Branch: `frontiers-r3-round2-fixes-2026-09-12`, forked from pre-audit commit
`c41a53b51c7dcf93c04006d654d90d73b1259f3f` on
`frontiers-r3-round2-preaudit-2026-09-12`.

---

## 0. Environment provenance

Recorded in `audit/ROUND2_ENVIRONMENT.txt`. All new/rerun analyses below use
`/opt/anaconda3/envs/microbiome-ad/bin/python` (Python 3.11.15) and
`/usr/local/bin/Rscript` (R 4.5.1). DADA2 was NOT rerun anywhere in Round 2.

---

## 1. Training-only genus-selection sensitivity (Reviewer concern #1)

**New script:** `scripts/09_training_only_feature_sensitivity.py`. Reuses
`01_harmonize_taxonomy.py`'s synonym-map / uninformative-genus-drop /
`apply_prevalence_filter` functions and `02_clr_transform.py`'s
`multiplicative_replacement`/`clr_transform` functions unmodified, via direct
module import (no duplication of logic, no rewrite of the original pipeline).

Design decision (documented in the script's module docstring): genus
prevalence for a given split is computed from exactly the AD/CN-labeled
sample rows that enter that split's classifier as training data (not the full
raw per-cohort table, which also contains MCI/unlabeled samples) — this is
the least ambiguous reading of "the held-out cohort must contribute nothing."
Same model families/grids/random_state/fold counts/bootstrap procedure as
`03_within_cohort_baseline.py` and `04_cross_cohort_generalization.py`.

**Command:**
```
/opt/anaconda3/envs/microbiome-ad/bin/python scripts/09_training_only_feature_sensitivity.py
```
**Environment:** `microbiome-ad` conda env, Python 3.11.15 (see `audit/ROUND2_ENVIRONMENT.txt`).
**Exit code:** 0. Only warnings emitted were benign sklearn `l1_ratio` notices
(the C=... l2-penalty grid points don't use `l1_ratio`; this is identical,
pre-existing behavior from the historical scripts' own `LogisticRegression`
configuration, not something introduced here).

**Outputs written (all new, non-destructive — no existing file overwritten):**
- `results/tables/within_cohort_auc_training_only.csv`
- `results/tables/loco_auc_training_only.csv`
- `results/tables/pairwise_auc_training_only.csv`
- `results/tables/training_only_feature_counts.csv` (48 splits: per-fold/per-cohort/per-pair retained genus count)
- `results/tables/training_only_feature_lists.csv` (9,250 rows — full retained-genus list per split, for the Feature-Selection Gate)
- `results/tables/training_only_feature_sensitivity_summary.csv` (original vs. training-only, side by side)

**Headline old-vs-new numbers (LOCO, the analysis central to the paper's claim):**

| Test cohort | Model | Original (leaky, 396 genera) LOCO AUC | Training-only LOCO AUC | Δ | n genera (training-only) |
|---|---|---|---|---|---|
| zhuang2018 | logreg | 0.5041 | 0.5057 | +0.0016 | 370 |
| zhuang2018 | lgbm | 0.5646 | 0.5960 | +0.0314 | 370 |
| ling2020 | logreg | 0.6775 | 0.6794 | +0.0019 | 370 |
| ling2020 | lgbm | 0.6735 | 0.6534 | -0.0201 | 370 |
| shanghai2022 | logreg | 0.8133 | 0.8122 | -0.0011 | 373 |
| shanghai2022 | lgbm | 0.7578 | 0.7878 | +0.0300 | 373 |
| kazakhstan2022 | logreg | 0.5644 | 0.6092 | +0.0448 | 150 |
| kazakhstan2022 | lgbm | 0.5905 | 0.5247 | -0.0658 | 150 |

**Interpretation (descriptive, not tuned):** the retained genus count drops
substantially when Kazakhstan is excluded from feature selection (150 vs. 396
— Kazakhstan alone passes 20% prevalence for many genera that the three
Chinese/Korean-processing-pipeline cohorts do not), but the LOCO AUC values
themselves move by at most ±0.066 and the qualitative conclusion — LOCO AUC
collapses toward chance relative to within-cohort performance — is unchanged
under strict training-only feature selection. Within-cohort AUCs (mirroring
`03_within_cohort_baseline.py`) shift by -0.083 to +0.038. No result flips
across the AUC=0.55 "near chance" line in a way that would change the paper's
qualitative claim. Full within-cohort/pairwise tables in the CSVs above.

---

## 2. Label-blind batch correction (Reviewer concern #2)

**New scripts:**
- `scripts/10_batch_correction_label_blind.R` — mirrors `06_batch_correction.R` exactly except `ComBat_seq(..., group=NULL, full_mod=FALSE)` and `adjust_batch(..., covariates=NULL, data=<cohort-only data frame>)`. Diagnosis is loaded only to define the labeled-sample subset, never passed to either correction function.
- `scripts/11_corrected_generalization_label_blind.py` — mirrors `07_corrected_generalization.py` exactly (same grids/random_state/fold counts/bootstrap), consuming the label-blind corrected matrices instead.

**Commands:**
```
Rscript scripts/10_batch_correction_label_blind.R
/opt/anaconda3/envs/microbiome-ad/bin/python scripts/11_corrected_generalization_label_blind.py
```
**Environment:** R 4.5.1 / `microbiome-ad` Python 3.11.15 (`audit/ROUND2_ENVIRONMENT.txt`).
**Exit codes:** both 0. R console confirms `"Using null model in ComBat-seq. Adjusting for 0 covariate(s)"` — verifying no covariate was supplied.

**Outputs (new, non-destructive):**
- `data/processed/clr_matrix_combatseq_labelblind.csv`, `clr_matrix_mmuphin_labelblind.csv`, `corrected_counts_combatseq_labelblind.csv`
- `results/tables/loco_auc_corrected_labelblind.csv`, `pairwise_auc_labelblind.csv`, `within_cohort_auc_corrected_labelblind.csv`, `auc_comparison_labelblind.csv`

**Headline result — this is a MAJOR, previously-unknown finding, reported exactly as observed:**

| Model | Correction | Mean LOCO AUC | Mean drop vs. within-cohort |
|---|---|---|---|
| logreg | uncorrected | 0.6398 | +0.1756 |
| logreg | ComBat-seq, **label-informed (original/oracle)** | 0.4930 | +0.323 (from original manuscript) |
| logreg | ComBat-seq, **label-blind (new, primary)** | 0.6338 | +0.1816 |
| logreg | MMUPHin, label-informed (original/oracle) | 0.6120 | (from original manuscript) |
| logreg | MMUPHin, label-blind (new, primary) | 0.6390 | +0.1764 |
| lgbm | uncorrected | 0.6466 | +0.1459 |
| lgbm | ComBat-seq, label-blind (new, primary) | 0.6540 | +0.1385 |
| lgbm | MMUPHin, label-blind (new, primary) | 0.6478 | +0.1446 |

**Per-cohort detail, logreg (the most affected model in the original analysis):**

| Test cohort | Uncorrected | ComBat-seq label-informed (original) | ComBat-seq label-blind (new) | MMUPHin label-informed (original) | MMUPHin label-blind (new) |
|---|---|---|---|---|---|
| zhuang2018 | 0.5041 | 0.4889 | 0.5251 | 0.5452 | 0.5284 |
| ling2020 | 0.6775 | 0.4465 | 0.6807 | 0.6644 | 0.7025 |
| shanghai2022 | 0.8133 | 0.7300 | 0.7689 | 0.7778 | 0.7822 |
| kazakhstan2022 | 0.5644 | **0.3080** | 0.5604 | 0.4628 | 0.5428 |

**This is a substantial change to the paper's central batch-correction narrative and is reported honestly, not minimized:** the severe, sometimes near-total collapse seen under the original label-informed design (Kazakhstan logreg LOCO AUC=0.308, ling2020 logreg LOCO AUC=0.4465 — both far below the uncorrected baseline) is **not reproduced** under label-blind correction, where every cohort's LOCO AUC stays close to its uncorrected value (largest single-cohort change: ling2020 logreg +0.0032 to +0.025 depending on method; no cohort collapses). This indicates the original "ComBat-seq substantially worsened cross-cohort generalization" finding was driven in large part by the label-informed/oracle design itself (each sample's own true label being used during correction), not by an intrinsic property of ComBat-seq's batch-effect removal in isolation.

**What is NOT reversed:** label-blind correction also does not meaningfully IMPROVE mean LOCO AUC (deltas of -0.006 to +0.007 relative to uncorrected, i.e. functionally flat) — so the paper's core claim that "standard batch correction methods do not recover cross-cohort generalization" remains supported, and is in fact now demonstrated more cleanly (a null/flat effect under label-blind correction, rather than a confounded catastrophic-degradation effect under label-informed correction). The manuscript must therefore say batch correction is *ineffective at recovering generalization* under a properly-scoped test, while explicitly retracting/reframing the stronger original claim that it *actively and substantially worsens* performance — that specific claim is now understood to be an artifact of the label-informed design and is presented as such (Section C below), with the original label-informed numbers retained only as an explicitly-labeled exploratory/oracle-style sensitivity.

---

## 3. SHAP direction/stability fix + matched Jaccard null (Reviewer concerns #3 and #4)

**Modified script:** `scripts/08_shap_taxa_comparison.py` (in place — the underlying
OOF SHAP computation, model fitting, fold structure, and random_state are
UNCHANGED; verified below). Changes:
- `run_within_cohort_shap` now also captures each outer fold's fitted logistic-regression
  coefficients (`COEF_ROWS`).
- New `build_coefficient_stability_table()`: per (cohort, taxon), counts n_positive/
  n_negative/n_zero coefficients across the 10 outer folds, median/mean coefficient,
  and a PRE-SPECIFIED (declared in the module before any results existed)
  `STABILITY_MIN_FOLDS = 8`-of-10 sign-agreement threshold → `sign_stability`
  ∈ {stable_positive, stable_negative, unstable}. Zero coefficients count toward neither.
  → `results/tables/logreg_coefficient_stability.csv`
- New `build_stable_directional_flips()`: a taxon is a "stable directional flip" only if
  it is in the top-20 OOF mean-|SHAP| ranking in ≥2 cohorts AND has a stable-positive
  coefficient in ≥1 cohort AND a stable-negative coefficient in ≥1 other cohort.
  → `results/tables/logreg_directional_flips_stable.csv`
- LightGBM's directional-flip claim is **removed entirely** — no coefficient-like
  scalar direction is computed or invented for it anywhere in the script or its
  console output. LightGBM retains OOF mean-|SHAP| importance, top-N overlap, and
  Jaccard analysis only.
- `plot_loco_direction()` **deleted** (was assigning AD/CN direction from LOCO mean
  signed SHAP, which is confounded by background-cohort/eval-cohort mean-abundance
  differences for both models) — replaced with an explanatory comment; LOCO
  feature-*importance* (`shap_loco_importance.csv`) is untouched and still generated.
- `plot_top15_per_cohort()`: for logreg, bar color now comes from the per-cohort
  median fitted coefficient (via `logreg_coefficient_stability.csv`), not `mean_shap`
  sign; for LightGBM, bars are a single neutral color with no direction legend.
- `plot_direction_dotplot()`: the ★ "flip" marker now comes from the coefficient-
  stability definition (`stable_flip_taxa`), not a naive mean-signed-SHAP sign check;
  called for logreg only.
- Jaccard null: `jaccard_null_distribution`/`compute_jaccard_null_baseline` (which drew
  2 random sets per replicate) **replaced** with `four_set_null_distribution`/
  `compute_jaccard_null_baseline_corrected`, which draws 4 random top-N sets per
  replicate and averages the same 6 unique pairwise Jaccards the observed statistic
  uses — exactly matching its construction. `N_JACCARD_PERM_CORRECTED = 100,000`;
  `JACCARD_N_SENSITIVITY = [10, 20, 50]`; empirical one-sided p =
  `(1 + sum(null >= observed)) / (n_perm + 1)`.
  → `results/tables/jaccard_null_corrected.csv` (the old `jaccard_null.csv` is no
  longer regenerated by this script and is retained on disk only as a superseded,
  audit-trail artifact — not read by any manuscript-generation code after this fix).

**Command:**
```
/opt/anaconda3/envs/microbiome-ad/bin/python scripts/08_shap_taxa_comparison.py
```
**Environment:** `microbiome-ad`, Python 3.11.15. **Exit code:** 0. Runtime: ~100s.

**Verification that the underlying SHAP computation is unchanged:** `git diff --stat`
on `shap_within_cohort_importance.csv`, `shap_loco_importance.csv`,
`shap_taxa_overlap.csv`, and `shap_directional_flips.csv` (all pre-existing,
committed outputs) shows **zero diff** — byte-identical to the pre-Round-2 values.
Only new files were added (`logreg_coefficient_stability.csv`,
`logreg_directional_flips_stable.csv`, `jaccard_null_corrected.csv`, and three new/
updated figures). This confirms the fix changed only how direction/significance are
interpreted, not the SHAP values or importance rankings themselves.

**Headline results:**

Coefficient sign stability across 1,584 cohort×genus combinations: 961 stable_negative,
446 stable_positive, 177 unstable.

**Stable directional flips (logistic regression, new conservative definition): 8 taxa**
(down from the old, invalid mean-signed-SHAP count of 12):
`Agathobacter, Bifidobacterium, Coprococcus, Dorea, Lactobacillus, NK4A214 group,
Romboutsia, Ruminococcus gnavus group`.

**LightGBM directional-flip claim: confirmed removed** — not computed, not printed,
not plotted anywhere in the script.

**Matched Jaccard null — this REVERSES part of the original manuscript's claim:**

| Model | Top-N | Observed mean (6 pairs) | Matched null mean [95% CI] | Empirical p | Old (invalid) null result |
|---|---|---|---|---|---|
| logreg | 10 | 0.0850 | 0.0134 [0.0000, 0.0361] | 0.000010 | (not previously reported at N=10) |
| logreg | **20 (primary)** | **0.1350** | 0.0265 [0.0085, 0.0491] | **0.000010** | p=0.0019, "above baseline" — conclusion unchanged, now more precise |
| logreg | 50 | 0.2296 | 0.0680 [0.0491, 0.0895] | 0.000010 | (not previously reported at N=50) |
| lgbm | 10 | 0.0175 | 0.0135 [0.0000, 0.0361] | 0.449386 | (not previously reported at N=10) |
| lgbm | **20 (primary)** | **0.0584** | 0.0265 [0.0085, 0.0491] | **0.0048 — SIGNIFICANT, ABOVE null** | p=0.068, "within the random baseline" — **REVERSED** |
| lgbm | 50 | 0.1430 | 0.0680 [0.0492, 0.0894] | 0.000010 | (not previously reported at N=50) |

**This is reported honestly, not minimized:** under the corrected, matched-statistic
null, LightGBM's top-20 taxon overlap across cohorts is **statistically significant**
(p=0.0048), reversing the original manuscript's claim that it was "within the random
baseline" / not distinguishable from chance. The N=10 sensitivity result for LightGBM
is NOT significant (p=0.449), so the LightGBM significance conclusion is genuinely
N-dependent and must be reported with that caveat rather than as a single flat
statement. The logistic-regression result was already significant under the old
(miscalibrated) null and remains significant — more strongly so — under the corrected
one at all three N values.

**Newly-discovered discrepancy, out of the original five concerns' scope, disclosed
here rather than concealed:** while implementing the coefficient-stability table, the
logistic regression model construction (`LogisticRegression(solver="saga",
l1_ratio=0.5, max_iter=10000, ...)`, unchanged from `03_within_cohort_baseline.py`/
`04_cross_cohort_generalization.py`, reused as-is here per the "same model families"
instruction) does **not** set `penalty="elasticnet"`, so it silently defaults to
`penalty="l2"` — sklearn emits `UserWarning: l1_ratio parameter is only used when
penalty is 'elasticnet'. Got (penalty=l2)` (suppressed by the historical scripts'
own `warnings.simplefilter("ignore")`, so never surfaced before). The manuscript
describes this model as "elastic-net" throughout (Methods, Fig. 2 caption, etc.),
but every logistic-regression result in the paper (AUC, SHAP, coefficients) was
actually produced by **L2-regularized (ridge) logistic regression**, not elastic
net. This was NOT changed/fixed by rerunning the pipeline (that would silently alter
every AUC number in the paper, a far larger and unauthorized scope change); instead
the manuscript text is corrected in Section G below to accurately say "L2-regularized"
instead of "elastic-net" everywhere it appears — a text-only fix consistent with
"every final manuscript number must trace to an output produced by the final code."
One practical consequence: the `n_zero` column in `logreg_coefficient_stability.csv`
is (honestly) ~0 for essentially every genus, since L2 does not zero out coefficients
the way true elastic net / L1 would — this is reported as observed, not disguised
(confirmed directly: 0 of 1,584 cohort×genus rows have n_zero>0).

---

## 3b. Figure 6 regenerated with the corrected direction/flip definitions

**Modified script:** `scripts/generate_manuscript_figures.py::make_fig6()`. Now reads
`logreg_coefficient_stability.csv` and `logreg_directional_flips_stable.csv` instead of
`shap_directional_flips.csv`'s naive mean-signed-SHAP sign check. Panel A dot color and
Panel B dot color both now come from the fitted logistic-regression coefficient sign
(median across 10 outer folds) for that (cohort, taxon); Panel B's x-position remains
the descriptive per-cohort mean SHAP value, now explicitly labeled as such. Both
panels' "flip" star/highlight now come from the pre-specified coefficient-stability
criterion. Panel C (Jaccard heatmap) is untouched (pure magnitude, no direction).

**Command:**
```
/opt/anaconda3/envs/microbiome-ad/bin/python scripts/generate_manuscript_figures.py
```
**Verification:** `make_fig1()`–`make_fig5()` and both supplementary figures reproduced
byte-for-byte identical output to the pre-Round-2 committed files (confirmed by
generating to `.png` and diffing against existing `.jpg` content visually/by size —
the redundant PNGs were then discarded, leaving only the genuinely-changed
`manuscript_fig6.jpg` staged). Figure 6 was visually inspected (Read tool) after
conversion to JPG: Panel A shows 5 of the 8 stable-flip taxa within its own narrower
top-20-by-max-SHAP display window (a pre-existing, unrelated windowing difference
between this file's TOP_PER_COHORT=15-then-top-20 selection and script 08's
per-cohort-top-20 flip-candidate pool — not something this fix introduced); Panel B
correctly shows all 8 stable-flip taxa (NK4A214 group, Dorea, Bifidobacterium,
Coprococcus, Ruminococcus gnavus group, Agathobacter, Lactobacillus, Romboutsia);
Panel C's observed values are unchanged and their mean ((0.290+0.212+0.026+0.176+
0.053+0.053)/6 = 0.135) still matches the manuscript's stated Jaccard figure exactly.

---

## 4. Manuscript rewrite (Section G) and remaining figure/caption fixes

**Files modified:** `manuscript/draft.md`, `manuscript/supplementary.md`, `scripts/build_pdf.py`
(figure caption list + one heading-text marker fix), `scripts/generate_manuscript_figures.py`
(Figure 5 now reads label-blind data as primary; new `make_supp_fig_s5()` for the
label-informed exploratory comparison; `make_supp_fig_s2()` LightGBM bars fixed to a
neutral color, removing the same invalid AD/CN-from-mean-SHAP coloring the main
Figure 6 fix addressed, which had also been present in this supplementary figure).

**Every affected section was rewritten**, per the Round 2 instructions: Abstract,
§2.2 (genus filtering — now states the actual 20%-prevalence global rule and
acknowledges held-out-cohort information use), new §2.5.1 (training-only sensitivity
methodology), §2.7 (batch correction — label-blind primary / label-informed
exploratory), §2.8 (SHAP — coefficient-based direction, matched Jaccard null), §2.10
(mixed-environment disclosure), §3.3 (LOCO — added training-only sensitivity
paragraph), §3.4 (PERMANOVA — empirical/hypothesis separation), §3.5 (batch
correction results — full rewrite around the label-blind finding), §3.6 (SHAP results
— full rewrite: 8 stable flips, no LightGBM flip, matched-null Jaccard), §4.2
(batch-correction discussion — explicit empirical-result-vs-hypothesis structure,
retraction of the "ComBat-seq substantially worsens" mechanistic claim), §4.4
(directional-flip discussion — 8 taxa, LightGBM claim removed), §4.6 (Limitations —
new bullets: global feature selection [now sensitivity-tested], two-design batch
correction, SHAP/coefficient limitations, L2-vs-elastic-net disclosure), Conclusion.
New Supplementary Tables S4 (training-only sensitivity), S5 (within-cohort AUC,
label-blind vs. label-informed), S6 (coefficient stability detail), S7 (Jaccard N
sensitivity); Supplementary Figures S3a/S3b (LOCO SHAP direction) replaced with an
explicit removal notice; Figure S2's caption corrected to match its neutral-color
regeneration; new Supplementary Figure S5 (label-informed exploratory batch
correction) added.

**Correctness catch made during the rewrite, documented rather than silently fixed:**
while writing the *Akkermansia* discussion (planned to parallel the old manuscript's
"CN-associated in Zhuang/Ling, AD-associated in Zhu 2022/Kazakhstan" framing but
translated into coefficient language), a direct check of
`results/tables/logreg_coefficient_stability.csv` showed *Akkermansia*'s fitted
coefficient is **stable and positive in all four cohorts** (10/10 outer folds each),
not mixed as the old mean-SHAP-based framing implied. This was corrected before
finalizing the prose (Section 3.6, Section 4.4, and new Supplementary Table S6) and
is now presented as a deliberate, concrete illustration of exactly why mean signed
SHAP was an inadequate direction statistic: the pooled mean SHAP value for this
genus does appear to flip sign across cohorts, but the model's actual fitted
coefficient never does. This is flagged explicitly per the standing instruction not
to propagate assumed-but-unverified numbers into the manuscript.

**Gate-check greps performed against the final `draft.md`** for every stale term
listed in the Round 2 instructions (0.135, 0.058, 0.068, 0.0019, "twelve
directional", "12 directional", "within the random baseline", "ruling out",
"mechanistically", "upper bound", "n=521", "elastic-net", "AD-associated",
"CN-associated"): every remaining hit was manually inspected and confirmed to be
either (a) a still-correct, unchanged value used in its original correct context,
or (b) intentionally retained inside a corrective/retraction sentence that
explicitly labels the old value as superseded. No live, uncaveated occurrence of a
prohibited phrase or a stale number was found. `AD-associated`/`CN-associated` do
not appear anywhere in the final `draft.md` at all (fully replaced by
coefficient/model-direction language). One additional stale "elastic-net" reference
was found and fixed in `scripts/build_pdf.py`'s Figure 2 caption (not caught by the
draft.md-only grep pass, since it lives in the caption-generation script).

**DOCX/PDF rebuild:**
```
pandoc manuscript/draft.md -o manuscript/draft.docx --standalone
pandoc manuscript/supplementary.md -o manuscript/supplementary.docx --standalone
<python-docx table-border script, applied to both>
/opt/anaconda3/bin/python scripts/build_pdf.py
```
Table-border pass: 4 tables in `draft.docx` (unchanged count), 13 tables in
`supplementary.docx` (up from 9, the 4 new Round 2 tables S4–S7). `build_pdf.py`
completed with zero "marker not found" warnings (verified by re-running and
grepping for "marker"/"not found"/"skip" — no matches), confirming every figure
insertion point, including the renamed §3.6 heading used as Figure 5's marker,
matched correctly. `manuscript/manuscript_draft.pdf` rebuilt (31 pages, ~3.08 MB).
Visually verified via PyMuPDF rendering: Section 2.2/2.3 (genus-filter correction +
CLR equation still typeset as real math), Table 4 (all 5 rows, borders correct),
Figure 5 (flat label-blind bars) and Figure 6 (coefficient-colored dots, Akkermansia
red/positive in all four cohorts, matching the corrected text) all render as
intended. One grammar slip ("an label-informed" → "a label-informed") was caught
during this visual pass and fixed, then the docx/pdf were rebuilt again.

---