# Reviewer 4 Revision — Final Validation Report

Branch `frontiers-r4-revision-2026-09-14`, forked from `main` at
`9111938e12281344c2d3d41e5c3ebf301cec2c5a` (tag `frontiers-round2-ready-2026-09-12`).

---

## M.1 Code validation

**Python:** `python -m py_compile scripts/*.py` — all Python scripts compile
cleanly, including the four new scripts (`12_zhu_robustness.py`,
`13_auc_difference_bootstrap.py`, `15_jaccard_prevalence_null.py`, plus all
pre-existing scripts).

**R:** `Rscript -e "parse('<file>.R')"` run against every `.R` script,
including the new `14_permanova_dispersion_sensitivity.R`. All parse without
syntax errors. (A naive grep for the string "error" against the parse output
of `06_batch_correction.R` and `10_batch_correction_label_blind.R` matches
their own legitimate `tryCatch(..., error = function(e) ...)` code, not a
parse failure; manually confirmed no `Error:` diagnostic was emitted by
`Rscript` for any file.)

**Execution:** all four new scripts were run to completion (not just
syntax-checked): `12_zhu_robustness.py` (background, ~24 min total: 61.4s
repeated-CV + 1348.0s permutation test), `13_auc_difference_bootstrap.py`,
`14_permanova_dispersion_sensitivity.R`, `15_jaccard_prevalence_null.py`, all
exiting 0 with no unhandled exceptions. Full commands, timings, and console
output are logged in `audit/REVIEWER4_IMPLEMENTATION_LOG_2026-09-14.md`.

## M.2 Output file existence

All 13 expected new/updated output files confirmed present on disk:
`results/tables/zhu_repeated_nested_cv.csv`,
`zhu_repeated_nested_cv_summary.csv`, `zhu_label_permutation.csv`,
`zhu_label_permutation_summary.csv`,
`within_vs_loco_auc_difference_bootstrap.csv`,
`within_cohort_auc_ci_stratified_bootstrap.csv`,
`loco_auc_ci_stratified_bootstrap.csv`, `betadisper_pairwise.csv`,
`permanova_chinese_miseq_sensitivity.csv`,
`betadisper_chinese_miseq_sensitivity.csv`,
`betadisper_chinese_miseq_pairwise.csv`,
`jaccard_prevalence_restricted_null.csv`, and
`results/figures/supp_zhu_robustness.png`.

## M.3 NaN / headline-result check

Checked every new summary CSV for unexpected NaNs. All are 0 except
`permanova_chinese_miseq_sensitivity.csv`, whose 2 NaNs are the `F_stat`/
`p_value` cells of the `Residual` row — expected and correct, since a
PERMANOVA residual term has no F-test of its own. No headline numeric result
(AUC, R², F, p) used anywhere in the manuscript is NaN, null, or a fallback
placeholder.

## M.4 Numeric cross-check: manuscript text vs. source CSVs

Every new number quoted in the manuscript and supplement was checked against
its source CSV. Spot-checks recorded here (full set checked, not just these):

- Zhu repeated-CV: logreg mean 0.9955/SD 0.0027, lgbm mean 0.9760/SD 0.0103 —
  matches `zhu_repeated_nested_cv_summary.csv` exactly.
- Zhu permutation: p=0.000999 for both models, null mean 0.4731 (logreg) /
  0.4719 (lgbm) — matches `zhu_label_permutation_summary.csv` exactly.
- Paired AUC-difference bootstrap: all 8 rows (observed delta, bootstrap mean,
  95% CI, proportion≤0) transcribed directly from
  `within_vs_loco_auc_difference_bootstrap.csv` with rounding only.
- Chinese-MiSeq PERMANOVA: cohort R²=0.058796910933659206→0.0588,
  diagnosis R²=0.022912163151868066→0.0229, F=10.049840991825823→10.05,
  F=7.832506599363772→7.83 — matches `permanova_chinese_miseq_sensitivity.csv`.
- Betadisper pairwise: all 10 rows (F, p, Holm, BH) transcribed directly from
  `betadisper_pairwise.csv`.
- Jaccard prevalence-restricted null: all 6 rows (observed, null mean, CI, p)
  transcribed directly from `jaccard_prevalence_restricted_null.csv`.
- Every within-cohort/LOCO AUC 95% CI in Table 2 and Table 3 was regenerated
  from `within_cohort_auc_ci_stratified_bootstrap.csv` /
  `loco_auc_ci_stratified_bootstrap.csv` and updated in the manuscript text
  and tables (previously reported CIs differed by ≤0.01 at 2-decimal
  rounding from the newly regenerated ones; both prose and Table 2/Table 3
  were updated for consistency, including the Section 3.2 prose sentence that
  had been missed in an initial editing pass and was caught and fixed during
  this validation).

No discrepancy between a manuscript-quoted number and its source CSV was
found during this check.

## M.5 Cross-reference integrity

Verified every `Table S#`/`Figure S#` citation in `draft.md` resolves to an
existing header in `supplementary.md`, and that no supplementary table/figure
number collides with a pre-existing one: new tables S8–S13 and new Figure S6
extend (not overwrite) the prior S1–S7/S1–S5 set. Verified new Methods
subsection numbers (2.4.1, 2.5.2) do not collide with the pre-existing 2.5.1,
and that all forward/backward section cross-references (e.g., "Section
2.5.2", "Section 3.4") point to sections that exist under those numbers.

## M.6 Prohibited-phrase sweep

Grepped `manuscript/draft.md` and `manuscript/supplementary.md` for the
following phrases, per the specified checklist, after the full edit pass:
"statistically reliable drop", "structurally much larger", "technical batch
effects versus biological heterogeneity", "AD-associated", "CN-associated",
"fundamental generalization problem", "cannot recover", "earlier version",
"we retract", "minimal assistance" — **zero matches for all ten**. Extended
the sweep (per Sections G/H/I of the task) to "readers of an earlier draft",
"mislabeled in earlier", "an earlier draft", "directionally consistent",
"association with AD status", "universal population effect", "cross-population
failure", "Round 2" — **zero matches for all**, except one benign remaining
use of "biological association" (Section 4.6, Limitations: "...is not a
direct estimate of biological association...") — this is a negation
disclaiming the association, the correct scientific framing, not a live
claim, and was intentionally left as-is.

## M.7 Figure regeneration

Checked whether any main-text figure embeds a numeric value that changed as
part of this revision. `scripts/generate_manuscript_figures.py` hardcodes
`r2_cohort = 0.1722` for Figure 4 — this is the full four-cohort model's
cohort R², which is unchanged by this revision (only a new, separate
Chinese-MiSeq-subset sensitivity value was added, reported in text/supplement,
not overwriting the main Figure 4 value). No other main-text figure encodes a
number that changed. **No main-text figure required regeneration.** The one
new figure, `results/figures/supp_zhu_robustness.png` (Supplementary Figure
S6), was generated fresh by `12_zhu_robustness.py` and confirmed present
(102 KB PNG).

## M.8 Document rebuild

- `manuscript/manuscript_draft.pdf` rebuilt via `scripts/build_pdf.py`
  (pandoc → HTML → Chrome headless print), 3,160 KB, built successfully with
  no script errors.
- `manuscript/draft.docx` and `manuscript/supplementary.docx` rebuilt via
  `pandoc --from markdown+raw_html --to docx --standalone --mathml`,
  56,619 bytes and 28,311 bytes respectively, no pandoc errors or warnings.
- **Limitation, disclosed honestly:** this environment has no PDF-rendering
  tool available (`pdftoppm`/poppler not installed, and no Python PDF library
  installed) to produce page-by-page raster images for visual inspection, and
  installing new system/Python packages was not attempted in order to avoid
  unrequested environment changes. Visual inspection was therefore performed
  at the markdown-source level (the actual content fed into the build
  pipeline, read in full via the file-reading tool across this entire
  revision) plus build-log verification (no pandoc/Chrome errors, output file
  sizes consistent with prior successful builds) rather than rendered-page
  inspection. If the author has poppler or a PDF viewer available locally, a
  final visual pass over `manuscript_draft.pdf` before submission is
  recommended as an added check beyond what this automated pass could do.

## M.9 Reviewer-3 points re-checked against the more rigorous standard

Per the task instruction to re-check every Reviewer 3 point, not just add
new Reviewer 4 content: the Reviewer 3 fixes (training-only feature
selection, label-blind batch correction, coefficient-based SHAP direction,
matched-statistic Jaccard null) were re-examined against Reviewer 4's more
rigorous statistical standard and found to still hold as methodologically
sound; the newly added Reviewer 4 analyses in this pass (formal paired AUC
bootstrap, prevalence-restricted Jaccard null, PERMANOVA/dispersion
sensitivity, Zhu repeated-CV/permutation) are additive extensions of, not
corrections to, the Reviewer 3 fixes, with one exception: the Reviewer
3-era "non-overlap of CIs indicates a statistically reliable drop" sentence
(Section 2.5) was itself a residual statistical-rigor gap that Reviewer 4's
point 2 (Section 3 of the response draft) directly identified and this pass
corrected.

## M.10 Skeptical self-audit: could Reviewer 4 reasonably repeat any original objection?

Going point-by-point through the task's original A–O specification and the
response draft:

- **Zhu robustness (B):** addressed with real repeated-CV and permutation
  results, both run under the leakage-safe training-only pipeline. A
  skeptical reviewer could ask why repeated-CV alone (without permutation)
  wouldn't have sufficed — but we anticipated this and ran both, with the
  manuscript explicit that repeated-CV is a stability check and permutation
  is the formal significance test. **Not repeatable.**
- **Formal AUC-difference test (C):** addressed with a real paired bootstrap,
  and — critically — the honest 3/8 non-significant result is reported
  rather than hidden. A skeptical reviewer could still ask why n=10,000
  bootstrap replicates and not a larger number, or why a percentile CI rather
  than BCa; these are defensible standard choices but not literally the only
  possible choices. **Largely addressed; a reviewer could ask for BCa CIs as
  a further refinement, which would be a reasonable minor follow-up, not a
  repeat of the original objection** (which was specifically about the
  CI-overlap heuristic, now removed).
- **PERMANOVA/dispersion (D):** addressed with pairwise post-hoc and the
  Chinese-MiSeq sensitivity subset, and the ~12x language is now qualified
  everywhere it appears. **Not repeatable** as originally framed; a
  remaining, inherent limitation (not something more code could fix) is that
  no sensitivity subset can fully isolate "technical" from "biological"
  without covariate-level metadata (diet, exact protocol details) that these
  public cohorts do not provide — this is now stated as a limitation, not
  hidden.
- **Jaccard prevalence sensitivity (E):** addressed, including the honest
  LightGBM reversal. **Not repeatable.**
- **Preprocessing leakage clarification (F):** addressed with an explicit
  per-step classification in Methods. **Not repeatable.**
- **SHAP terminology (G):** addressed via full-text sweep (Section M.6, zero
  matches on all prohibited terms) and shortened dietary speculation.
  **Not repeatable.**
- **Scope/geography (H):** addressed; Abstract/Introduction/Conclusion now
  use cohort/dataset-scoped language. **Not repeatable.**
- **Revision-changelog language (I):** removed from the manuscript body;
  confirmed via grep sweep. **Not repeatable.**
- **AI disclosure (J):** rewritten to name the actual tool/version/model and
  describe its actual role, per Frontiers policy. **Not repeatable** as
  originally framed ("minimal assistance" no longer appears); a reviewer
  could still note the disclosure cannot describe AI use outside the
  author's own visible record (e.g., undocumented ChatGPT use) — this is
  disclosed as an explicit limitation of the disclosure itself, which is the
  honest position available.
- **Manuscript/supplement completeness (K):** Abstract, Introduction,
  Methods 2.2–2.10, Results 3.2–3.6, Discussion 4.1–4.6, Limitations,
  Conclusion, and supplementary.md were all updated; new supplementary
  tables S8–S13 and Figure S6 added. **Not repeatable.**

**Remaining issue a reviewer could reasonably still raise, not fully
resolvable within this pass:** the underlying dataset remains five public
cohorts with no shared covariate metadata (diet, medication, detailed
protocol parameters) — every sensitivity analysis in this revision quantifies
*how much* a result depends on cohort composition or null-model choice, but
none of them can add the missing covariate data itself. This is stated
directly in the Limitations section and is not something further
re-analysis of the existing public data can fix; it would require new,
prospectively collected, covariate-rich multi-cohort data, which is exactly
what Section 4.5's clinical-translation recommendations call for.

---

**Overall assessment:** every A–O task item has been implemented with real
executed analyses (no fabricated or estimated numbers), the manuscript and
supplement have been rewritten to incorporate all results honestly —
including the two findings that reverse or soften a previously reported
result (LightGBM Jaccard significance, and the ~12x PERMANOVA disparity) —
and the prohibited-phrase and cross-reference sweeps are clean. The one
disclosed gap is the inability to visually render PDF pages in this
environment (M.8); everything else in Section M has been completed and
verified against source data.
