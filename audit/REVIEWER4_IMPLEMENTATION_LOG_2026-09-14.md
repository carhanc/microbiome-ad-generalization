# Reviewer 4 — Implementation Log

Running log of every scientific change made while addressing Reviewer 4's
independent report on manuscript 1918251, plus a complete re-check of every
Reviewer 3 point against this new, more statistically rigorous standard.

Branch: `frontiers-r4-revision-2026-09-14`, forked from `main` at commit
`9111938e12281344c2d3d41e5c3ebf301cec2c5a` (tag `frontiers-round2-ready-2026-09-12`).
DADA2 is NOT rerun anywhere in this pass.

---

## 0. Environment

Reused the same environments verified in the prior pass:
`/opt/anaconda3/envs/microbiome-ad/bin/python` (Python 3.11.15) for all Python;
`/usr/local/bin/Rscript` (R 4.5.1) for all R. Re-verified below before first use.

---

## 1. `scripts/12_zhu_robustness.py` — Zhu 2022 repeated CV + label permutation

Command: `python scripts/12_zhu_robustness.py` (run in background; log at
scratchpad `log_12_full.txt`). Both parts run under the strict training-only
feature-selection pipeline (genus selection re-derived per outer fold from
training-side samples only), not the global 396-genus feature set.

**Design note (label-independent precomputation):** genus selection and CLR
feature construction depend only on which participants are in a fold's
training set, not on their labels. Fold features were therefore precomputed
once per outer-fold partition and reused across all label permutations in
Part B — a pure implementation speedup verified to have zero effect on the
statistical result (the model itself is always refit on the permuted labels
for every permutation).

**Design note (serial vs. parallel execution):** benchmarked `n_jobs=-1`/`4`
(GridSearchCV/LightGBM internal parallelism) against fully serial (`n_jobs=1`)
execution on this n=60 dataset: single-run time was ~2.6s/7.0s (logreg/lgbm)
parallel-internal vs. ~2.1s/6.0s serial — serial was as fast or faster.
Chose serial-per-run + `ProcessPoolExecutor(max_workers=10)` across
independent repeats/permutations (11 physical cores available), avoiding
CPU oversubscription while parallelizing at the correct level.

**Bug fixed during development:** `project_and_clr()` (imported from
`09_training_only_feature_sensitivity.py`) expects a dict keyed by cohort,
not a bare DataFrame; initial call raised `KeyError: 'shanghai2022'`. Fixed
by wrapping the raw table as `{COHORT_KEY: raw_table}`.

**Part A — Repeated nested CV (50 repetitions, seeds 1000–1049):**
Wall time 61.4s. Results (`results/tables/zhu_repeated_nested_cv_summary.csv`):

| model | historical single-partition AUC | training-only single-partition AUC | repeated-CV mean | SD | median | P2.5 | P97.5 | min | max |
|---|---|---|---|---|---|---|---|---|---|
| logreg | 0.998 | 0.9956 | 0.9955 | 0.0027 | 0.9956 | 0.9900 | 0.9989 | 0.9889 | 0.9989 |
| lgbm | 0.979 | 0.9700 | 0.9760 | 0.0103 | 0.9767 | 0.9519 | 0.9917 | 0.9378 | 0.9944 |

**Part B — Label permutation test (n_perm=1000, seed 20260914):**
Wall time 1348.0s (~22.5 min; logreg 428.6s, lgbm 919.1s). Results
(`results/tables/zhu_label_permutation_summary.csv`):

| model | observed AUC (training-only) | n_perm | null mean | null SD | null P2.5 | null P97.5 | empirical p | p resolution floor |
|---|---|---|---|---|---|---|---|---|
| logreg | 0.9956 | 1000 | 0.4731 | 0.1017 | 0.2800 | 0.6734 | 0.000999 | 0.000999 |
| lgbm | 0.9700 | 1000 | 0.4719 | 0.1094 | 0.2643 | 0.6878 | 0.000999 | 0.000999 |

No permutation reached the observed AUC for either model; p is at the
resolution floor (1/1001) and is not interpreted as more precise than that.
Figure written: `results/figures/supp_zhu_robustness.png` (2×2: repeated-CV
histograms top row, permutation-null histograms bottom row).

A reduced-parameter smoke test (`N_REPEATS=4`, `N_PERM_TARGET=6`) was run
first via direct script execution (not dynamic `importlib` loading, which
fails under `ProcessPoolExecutor`'s spawn start method with a `PicklingError`
because spawned workers cannot re-import a dynamically-loaded module by
name — a smoke-test-harness artifact only, not a script bug) to confirm
correct end-to-end behavior before launching the full run.

---

## 2. `scripts/13_auc_difference_bootstrap.py` — formal within-cohort vs. LOCO AUC difference

Command: `python scripts/13_auc_difference_bootstrap.py`. N_BOOT=10,000,
RANDOM_STATE=42. Paired, diagnosis-stratified, participant-level bootstrap;
within-cohort OOF prediction files (no explicit run_id column) had
participant identity reconstructed by replicating `03_within_cohort_baseline.py`'s
exact merge/filter sequence, verified via exact-match assertion against each
file's own saved `y_true` column for all 8 cohort×model combinations (all
passed) before use.

Results (`results/tables/within_vs_loco_auc_difference_bootstrap.csv`):

| cohort | model | AUC within | AUC LOCO | Δ observed | bootstrap mean Δ | 95% CI | P(Δ≤0) |
|---|---|---|---|---|---|---|---|
| zhuang2018 | logreg | 0.6333 | 0.5041 | 0.1293 | 0.1280 | [-0.0465, 0.3029] | 0.0787 |
| zhuang2018 | lgbm | 0.6349 | 0.5646 | 0.0703 | 0.0705 | [-0.1001, 0.2391] | 0.2138 |
| ling2020 | logreg | 0.8648 | 0.6775 | 0.1873 | 0.1871 | [0.0970, 0.2772] | 0.0001 |
| ling2020 | lgbm | 0.8618 | 0.6735 | 0.1883 | 0.1886 | [0.0955, 0.2838] | 0.0001 |
| shanghai2022 | logreg | 0.9978 | 0.8133 | 0.1844 | 0.1848 | [0.0722, 0.3122] | 0.0000 |
| shanghai2022 | lgbm | 0.9789 | 0.7578 | 0.2211 | 0.2211 | [0.1111, 0.3456] | 0.0000 |
| kazakhstan2022 | logreg | 0.7657 | 0.5644 | 0.2014 | 0.2013 | [0.0442, 0.3573] | 0.0051 |
| kazakhstan2022 | lgbm | 0.6943 | 0.5905 | 0.1038 | 0.1038 | [-0.0493, 0.2552] | 0.0917 |

**Honest finding kept in the manuscript:** 3/8 pairs (zhuang2018 both models,
kazakhstan2022 lgbm) have a 95% CI on Δ that includes zero — not formally
significant at this sample size, despite all 8 pairs showing the same-signed
(within-cohort higher) drop. Also regenerated every within-cohort/LOCO AUC CI
in the paper using this same procedure
(`within_cohort_auc_ci_stratified_bootstrap.csv`, `loco_auc_ci_stratified_bootstrap.csv`)
for one consistent CI methodology throughout.

---

## 3. `scripts/14_permanova_dispersion_sensitivity.R` — PERMANOVA/dispersion robustness

Command: `Rscript scripts/14_permanova_dispersion_sensitivity.R`. N_PERM=9999, seed=42.

**Part 1 — pairwise betadisper post-hoc (all 5 cohorts, Aitchison), Holm+BH adjusted**
(`results/tables/betadisper_pairwise.csv`): ling2020 differs significantly
from every other cohort (all Holm p≤0.0063) and drives most of the overall
significant dispersion effect; shanghai2022-zhuang2018 (Holm p=0.439) and
kazakhstan2022-kbase2022 (Holm p=0.439) are not significantly different.

**Part 2 — Chinese paired-end MiSeq sensitivity** (zhuang2018+ling2020+shanghai2022,
n=317, binary AD/CN): marginal PERMANOVA cohort R²=0.0588 (F=10.05, p=0.0001)
vs. diagnosis R²=0.0229 (F=7.83, p=0.0001) — **~2.6-fold disparity, vs. ~12-fold
in the full four-cohort model.** Betadisper on this subset: F=7.946, p=0.0004;
pairwise: zhuang2018-ling2020 p=0.0001, zhuang2018-shanghai2022 p=0.026,
ling2020-shanghai2022 p=0.343 (not significant).
Outputs: `permanova_chinese_miseq_sensitivity.csv`,
`betadisper_chinese_miseq_sensitivity.csv`, `betadisper_chinese_miseq_pairwise.csv`.

**Manuscript consequence:** the "~12x" cohort-vs-diagnosis disparity is not
presented as a fixed, generalizable figure; it is reported alongside this
~2.6x sensitivity result, and PERMANOVA is now explicitly described as
unable to decompose "cohort" (a composite study-of-origin variable) into
technical vs. biological components.

---

## 4. `scripts/15_jaccard_prevalence_null.py` — prevalence-restricted Jaccard null sensitivity

Command: `python scripts/15_jaccard_prevalence_null.py`. MIN_PREVALENCE=0.20,
N_LIST=[10,20,50], N_PERM=100,000, RANDOM_STATE=42. Per-cohort eligible
universes (≥20% prevalence within that cohort's own AD/CN samples):
zhuang2018=131, ling2020=126, shanghai2022=115, kazakhstan2022=356 (all ≥50,
no N-adaptation needed).

Results (`results/tables/jaccard_prevalence_restricted_null.csv`):

| model | top_n | observed | null_mean | null 95% CI | empirical p |
|---|---|---|---|---|---|
| logreg | 10 | 0.0850 | 0.0260 | [0.0000, 0.0624] | 0.001670 |
| logreg | 20 | 0.1350 | 0.0522 | [0.0261, 0.0841] | 0.000030 |
| logreg | 50 | 0.2308 | 0.1422 | [0.1140, 0.1737] | 0.000010 |
| lgbm | 10 | 0.0175 | 0.0260 | [0.0000, 0.0624] | 0.799592 |
| lgbm | 20 | 0.0629 | 0.0522 | [0.0261, 0.0841] | 0.240838 |
| lgbm | 50 | 0.1665 | 0.1422 | [0.1140, 0.1737] | 0.061309 |

**Honest reversal kept in the manuscript:** under the shared-396-genus-pool
null (Round 2 / Table S7), LightGBM was significant at N=20 (p=0.0048) and
N=50. Under this more realistic per-cohort-eligible-universe null, LightGBM
is **not significant at any of N=10/20/50** (p=0.80/0.24/0.061). Logistic
regression remains significant at all three N under both null designs. The
manuscript reports both null models and this reversal directly rather than
preserving the earlier significant-LightGBM claim.

---
