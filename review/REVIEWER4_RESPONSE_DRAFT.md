# Response to Reviewer 4

Manuscript 1918251, *Frontiers in Aging Neuroscience*. We thank Reviewer 4 for a
careful, favorable, and statistically rigorous independent review. Reviewer 4's
central assessment — that the main empirical finding of reduced cross-cohort
performance is supported and robust — is one we share, and every point raised
below has been addressed with either a new quantitative robustness analysis or
a manuscript wording correction. Where a point overlaps with changes already
made in the unpublished, interim revision prepared in response to Reviewer 3,
we say so explicitly below, together with what has been added or tightened
specifically in response to Reviewer 4. Reviewer 4 evaluated the version of
the manuscript available at the time of their review; nothing below implies
they missed anything in that version.

---

## 1. Zhu 2022: is the near-perfect within-cohort AUC (n=60) robust?

**Reviewer request (paraphrased):** With only 60 participants, a single
nested-CV partition's AUC of 0.998/0.979 could be an artifact of a favorable
fold assignment or of leakage; the manuscript should assess sensitivity to
fold assignment and formally test the result against a permutation null.

**Response.** Agreed, and this was not addressed in the interim Reviewer-3
revision. We added two new analyses (`scripts/12_zhu_robustness.py`,
Section 2.4.1), both run under the strict training-only feature-selection
pipeline (not the global feature set) so neither could leak held-out-fold
information:

- **Repeated nested cross-validation** (50 repetitions, different
  pre-specified seed each time): mean pooled-OOF AUC 0.9955 (SD 0.0027) for
  logistic regression and 0.9760 (SD 0.0103) for LightGBM — closely matching
  the single-partition value and showing the result is not an artifact of one
  fold assignment. We are explicit in the manuscript that this is a
  fold-assignment stability check on the same 60 participants, not
  independent replication.
- **Label permutation test** (1,000 permutations, same strict pipeline,
  fixed partition): null AUC centered near chance for both models (mean
  ≈0.47, 95% range ≈[0.27–0.69]); no permutation among the 1,000 reached the
  observed AUC for either model, giving p at the resolution floor for this
  permutation count (p ≈ 0.000999 for both models). We report this p-value at
  its actual resolution and do not claim greater precision than that.

**Manuscript location:** Section 2.4.1 (Methods), Section 3.2 (Results),
Abstract, Supplementary Table S8/S9, Supplementary Figure S6.

**Did the scientific conclusion change?** No — the near-perfect Zhu 2022
within-cohort AUC is now shown to be robust to fold assignment and clearly
distinguishable from chance, strengthening rather than altering the reported
result. The manuscript's cautionary framing that this figure is based on a
small sample (Section 4.3) is retained regardless.

---

## 2. Formal statistical test for the within-cohort vs. LOCO AUC difference

**Reviewer request (paraphrased):** Non-overlap of two independently computed
bootstrap confidence intervals is not a formal statistical test of a nonzero
difference; the manuscript should not describe it as indicating a
"statistically reliable drop," and should instead report a formal, paired
comparison.

**Response.** Agreed. This was not addressed in the interim Reviewer-3
revision, which used the CI-overlap heuristic (Section 2.5's original text:
"the non-overlap of within-cohort and LOCO CIs serves as an indicator of a
statistically reliable drop"). We implemented a paired, diagnosis-stratified,
participant-level bootstrap (`scripts/13_auc_difference_bootstrap.py`,
10,000 replicates, Section 2.5.2): for each cohort–model pair, the same
resampled participant indices are used to compute both the within-cohort and
LOCO AUC in each replicate, so the resulting delta distribution is a genuine
paired comparison. We report the observed delta, the bootstrap mean, the 95%
CI on the delta, and the proportion of replicates with delta≤0 — the last
explicitly described as an empirical bootstrap quantity, not an exact
parametric p-value.

**Result, reported honestly rather than selectively:** five of eight
cohort–model pairs (Ling 2020 both models, Zhu 2022 both models, Kazakhstan
logistic regression) have a 95% CI on the delta that excludes zero. The
remaining three (Zhuang 2018 both models, Kazakhstan LightGBM) show the same
directional pattern (within-cohort higher than LOCO) but a CI that includes
zero, and are not formally significant at this sample size. We revised the
manuscript to state this distinction explicitly (Section 3.3, Table 3,
Section 4.1, Conclusion) rather than treating the mean LOCO drop as equally
well-supported for every individual pair. We also removed the "non-overlap
... indicates a statistically reliable drop" wording from Methods (Section
2.5) and regenerated every AUC 95% CI in the paper (Table 2, Table 3, Table
S1) using this same diagnosis-stratified participant bootstrap, so every CI
in the manuscript now derives from one consistent, explicitly specified
procedure.

**Manuscript location:** Section 2.5.2 (Methods), Section 3.3 (Results),
Table 3, Abstract, Conclusion, Limitations, Supplementary Table S10.

**Did the scientific conclusion change?** Partially and honestly: the
aggregate finding (classifiers trained in one cohort under-perform on
held-out cohorts) remains well-supported, but we no longer imply that all
eight individual cohort–model comparisons are equally statistically robust —
three are not, at this sample size, and the manuscript now says so.

---

## 3. PERMANOVA/dispersion robustness and the "~12-fold" cohort-vs-diagnosis figure

**Reviewer request (paraphrased):** "Cohort" bundles technical and biological
differences that PERMANOVA cannot separate; the manuscript's Introduction
should not describe PERMANOVA as distinguishing "technical batch effects
versus biological heterogeneity," and the ~12-fold cohort-vs-diagnosis
variance disparity should be checked for sensitivity to which cohorts are
compared and which cohort pairs drive the significant dispersion difference.

**Response.** Agreed, and this was not addressed in the interim Reviewer-3
revision, which retained the "technical batch effects versus biological
heterogeneity" framing in the Introduction and presented the ~12-fold figure
without a cohort-composition sensitivity check. We made two additions
(`scripts/14_permanova_dispersion_sensitivity.R`, Section 2.6):

- **Pairwise post-hoc beta-dispersion comparisons** (all 10 cohort pairs,
  Holm/BH-adjusted): the overall significant dispersion difference is driven
  predominantly by Ling 2020, which differs significantly from every other
  cohort (all Holm-adjusted p≤0.0063); two pairs (Zhu 2022 vs. Zhuang 2018;
  Kazakhstan vs. Kim/KBASE) are not significantly different from one another.
- **A "more technically comparable" sensitivity subset** (Zhuang 2018 + Ling
  2020 + Zhu 2022 — all Chinese, paired-end MiSeq V3–V4 studies, still not
  technically identical): marginal PERMANOVA in this subset gives cohort
  R²=0.0588 vs. diagnosis R²=0.0229 — **an approximately 2.6-fold disparity,
  not 12-fold.** We report this directly as evidence that the ~12-fold figure
  is not a fixed, generalizable property of "cohort effects" in this kind of
  data; it depends substantially on which cohorts are compared, and in
  particular on the inclusion of Kazakhstan (a different sequencing platform)
  and Kim/KBASE (a different country).

We rewrote the Introduction's PERMANOVA framing (removing the "technical
batch effects versus biological heterogeneity" language), added an explicit
statement in Methods (Section 2.6) that cohort is a composite study-of-origin
variable PERMANOVA cannot decompose, and added this caveat to every place the
~12-fold figure appears (Abstract, Section 3.4, Section 4.2, Conclusion,
Limitations).

**Manuscript location:** Section 2.6 (Methods), Section 3.4 (Results),
Section 1 (Introduction), Limitations, Abstract, Conclusion, Supplementary
Table S11/S12.

**Did the scientific conclusion change?** Yes, honestly: we no longer present
the ~12-fold cohort-vs-diagnosis disparity as a stable, generalizable figure.
The qualitative finding that cohort explains substantially more variance than
diagnosis in this dataset is unchanged, but its magnitude is now explicitly
shown to be composition-dependent, and PERMANOVA is no longer described as
separating technical from biological sources anywhere in the manuscript.

---

## 4. Jaccard overlap null model: is the uniform-pool null realistic?

**Reviewer request (paraphrased):** The Jaccard-overlap null draws candidate
genera as if every genus were equally available in every cohort, but per-
cohort genus prevalence varies substantially; a more realistic,
prevalence-restricted null should be checked, and the manuscript should not
assume the previously reported significance survives.

**Response.** Agreed. A different null-model correction was already made in
the interim Reviewer-3 revision — the statistic-matching fix ensuring the
null draws four random top-N sets and averages the same six pairwise Jaccards
as the observed statistic, rather than a mismatched single-pair null. That
fix is distinct from, and does not address, the prevalence-availability issue
Reviewer 4 raises. We added a second, prevalence-restricted null
(`scripts/15_jaccard_prevalence_null.py`, Section 2.8): each cohort's top-N
ranking and random draws are now restricted to genera reaching ≥20%
prevalence within that cohort's own AD/CN samples (eligible universe sizes:
115–131 for the three MiSeq cohorts, 356 for Kazakhstan), rather than the
shared 396-genus pool.

**Result, reported as a genuine reversal:** logistic regression remains
significant at all three tested N (10, 20, 50) under this more realistic
null. **LightGBM's overlap, which was significant at N=20/50 under the
shared-pool null, is not significant at any of N=10/20/50 under the
prevalence-restricted null** (p=0.80, 0.24, 0.061). We do not preserve the
earlier significant-LightGBM finding; the manuscript now reports both null
models side by side and states plainly that LightGBM's apparent cross-cohort
top-predictor overlap does not survive this more rigorous check.

**Manuscript location:** Section 2.8 (Methods), Section 3.6 (Results),
Abstract, Limitations, Supplementary Table S13 (new), Supplementary Table S7
(updated to cross-reference S13).

**Did the scientific conclusion change?** Yes, for LightGBM specifically: we
no longer claim LightGBM shows a statistically detectable cross-cohort
top-predictor overlap. Logistic regression's overlap finding is unchanged and
robust to both null designs.

---

## 5. Preprocessing leakage: is every step's leakage status stated explicitly?

**Reviewer request (paraphrased):** The Methods should explicitly classify
every preprocessing operation's leakage status (which steps are label-free
and cohort-independent, which use information from the eventual test
cohort, and which have been checked with a training-only sensitivity
analysis) rather than leaving this to be inferred.

**Response.** Agreed. We added an explicit classification of every
preprocessing step to Section 2.2: DADA2 (cohort-specific, label-free, no
leakage pathway), global genus selection (uses all cohorts, transparently
identified as not leakage-free, bounded by the training-only sensitivity
analysis), zero-replacement/CLR (deterministic, sample-local, no leakage
pathway), and model tuning/fitting (training-side samples only). No new
preprocessing was run for this point; the existing training-only sensitivity
analysis (Section 2.5.1, unchanged from the interim revision) already
quantifies the impact of the one step that is not leakage-free.

**Manuscript location:** Section 2.2 (Methods).

**Did the scientific conclusion change?** No — this is a clarity/completeness
correction to already-accurate Methods text.

---

## 6. SHAP terminology: "AD-associated," "biological association," "directionally consistent"

**Reviewer request (paraphrased):** Language implying a validated biological
association with AD status (e.g., "AD-associated," "CN-associated,"
"directionally consistent... association with AD status") should not be used
for a fitted-model coefficient pattern that has not been independently
validated or covariate-adjusted.

**Response.** Agreed. The move away from mean-signed-SHAP-based direction
language to fitted-coefficient-based language was made in the interim
Reviewer-3 revision (Section 2.8's coefficient-stability definition). In
response to Reviewer 4, we did a full sweep of the manuscript and
supplementary text for "AD-associated," "CN-associated," "association with
AD status," "biological association" (as a live claim), and "directionally
consistent," and rewrote every remaining instance to use only "positive/
negative fitted coefficient" and "higher/lower predicted AD log-odds
conditional on the fitted model." We also substantially shortened the
Akkermansia/Lactobacillus dietary-speculation passages (Sections 3.6, 4.4),
retaining only the fitted-coefficient pattern, a brief statement that diet is
a plausible unmeasured confounder, and an explicit statement that no dietary
covariate data exist and no mechanistic attribution is made.

**Manuscript location:** Section 3.6, Section 4.4 (rewritten and shortened).

**Did the scientific conclusion change?** No — this is a terminology and
scope-of-claim correction; the underlying coefficient values are unchanged.

---

## 7. Scope and geography: overclaiming beyond the cohorts studied

**Reviewer request (paraphrased):** Language such as "cohort effect is
structurally much larger than the disease effect" and framing the finding as
a "universal population effect" or "cross-population failure" overstates what
five cohorts from three countries can establish; the conclusion should be
scoped to the cohorts and datasets actually evaluated.

**Response.** Agreed. We removed "structurally much larger" from the Abstract
and replaced it with wording that separates the empirical variance-magnitude
finding from the PERMANOVA composite-variable and dispersion-heterogeneity
caveats (see point 3 above). The Abstract, Introduction, and Conclusion now
state the finding as "substantial cross-cohort generalization limitations
among the evaluated publicly available 16S datasets" rather than as a
universal or population-level claim.

**Manuscript location:** Abstract, Introduction, Conclusion.

**Did the scientific conclusion change?** No — this is a scope-of-claim
correction consistent with the cohorts actually studied.

---

## 8. Revision-changelog language in the manuscript body

**Reviewer request (paraphrased; also a general manuscript-hygiene point):**
Phrases like "an earlier version of this manuscript," "we retract," "readers
of an earlier draft should disregard," and "Round 2 revision" belong in a
reviewer-response letter, not in the manuscript itself, which should simply
state the current, correct methodology.

**Response.** Agreed. We removed this language throughout the manuscript and
supplementary materials (Sections 2.4, 2.9, 2.10, 3.6, 4.2, 4.4, Limitations,
and multiple supplementary table/figure headers), restating each point as a
direct description of the current method or finding. The full history of
what changed and why remains fully transparent in this response document and
in `review/REVIEWER3_ROUND2_RESPONSE_DRAFT.md` and the audit logs, which are
not part of the manuscript itself.

**Manuscript location:** throughout; see `audit/REVIEWER4_IMPLEMENTATION_LOG_2026-09-14.md`
for the specific line-level changes.

**Did the scientific conclusion change?** No — wording only.

---

## 9. AI tool disclosure accuracy

**Reviewer request (paraphrased; also a Frontiers policy requirement):** The
disclosure should accurately name and describe the generative-AI tools used,
including model/version, rather than describing the assistance as "minimal."

**Response.** Agreed. The disclosure previously stated "minimal assistance
from Claude (Anthropic)," which understated the actual use of AI tools in
drafting code, debugging, drafting text, and preparing reviewer responses. We
rewrote Section 2.9 to name Claude Code (Anthropic; CLI version 2.1.104,
Claude Sonnet 5 model) and describe its actual role, state that we have no
evidence of ChatGPT/OpenAI tool use and are not fabricating such a
disclosure, and state that all code executed in the author's own environment,
all results were checked against source outputs by the author, all
methodological decisions were made and are owned by the author, and no AI
tool is credited as an author.

**Manuscript location:** Section 2.9 (Methods), Acknowledgements.

**Did the scientific conclusion change?** No.

---

## Summary of net effect on conclusions

Two findings changed as a direct, honest result of these more rigorous
checks, and we want these flagged explicitly rather than left for a reader to
notice on their own:

1. **LightGBM's cross-cohort top-predictor Jaccard overlap is no longer
   reported as statistically significant** under a more realistic,
   prevalence-restricted null model (point 4). Logistic regression's overlap
   finding is unchanged.
2. **The ~12-fold cohort-vs-diagnosis PERMANOVA variance disparity is now
   explicitly shown to be sensitive to cohort composition** (falling to
   ~2.6-fold in a more technically comparable three-cohort subset) and is no
   longer presented as a fixed, generalizable figure (point 3).

Additionally, **three of eight within-cohort-vs-LOCO AUC comparisons do not
reach formal statistical significance** under the new paired bootstrap test,
though all eight show the same directional pattern (point 2). None of these
three changes alters the paper's central, Reviewer-4-endorsed conclusion that
cross-cohort generalization is substantially and measurably limited in this
dataset; if anything, the paper is now more conservative and more precisely
scoped than the version Reviewer 4 reviewed.
