# kbase2022: Per-Sample Diagnosis Labels — IRB-Restricted, Not Publicly Available

**Date documented:** 2026-06-21  
**Cohort:** Kim/KBASE 2022 (PRJEB50447, ENA)  
**Impact on analysis:** Excluded from supervised classifier analyses (Phases 2–3); included in unsupervised analyses (Phase 4+).

---

## What was searched

All publicly accessible metadata sources were exhaustively checked:

1. **ENA Portal filereport** (`https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB50447&...`)  
   Fields available: `run_accession`, `sample_accession`, `sample_alias`, `sample_title`, `instrument_model`, `library_layout`, `library_strategy`  
   Result: No diagnosis or amyloid-status field. `sample_alias` values follow a sequential SP-number pattern (e.g., SP5100_1_R1) with no encoded label.

2. **ENA BioSample XML** (checked for individual sample accessions)  
   Result: No `disease`, `diagnosis`, `condition`, or `amyloid` attributes present.

3. **NCBI SRA runinfo** (alternative source for the same accessions)  
   Result: Same metadata, no clinical labels.

4. **Sample title strings** in ENA filereport  
   Values are uninformative instrument/library descriptions, not subject identifiers.

5. **SP-number pattern analysis** (sample_alias like `SP5100_1_R1`)  
   The numeric range spans SP5100–SP5178 (roughly 78 samples). No pattern in these numbers maps to a binary AD/CN split — the numbers are sequential and do not encode diagnosis.

---

## What the paper reports

Kim et al. 2022 (PRJEB50447) reports **18 amyloid-PET positive preclinical AD** and **60 amyloid-PET negative cognitively normal controls** at the group level. Individual sample-level amyloid-PET labels, MMSE scores, and other clinical attributes are not deposited in any public archive. This is consistent with IRB protocols that permit group-level reporting in publications but restrict individual-level clinical data from public databases.

---

## Consequence for this project

| Analysis type | kbase2022 included? | Reason |
|---|---|---|
| Phase 2: Within-cohort baseline classifiers | **NO** | Requires per-sample diagnosis labels for training/testing |
| Phase 3: Cross-cohort generalization (LOCO) | **NO** | Cannot be a training or test set without labels |
| Phase 4: PERMANOVA variance decomposition | **YES** | Only requires cohort identity, not per-sample labels |
| Phase 5: Batch correction + re-test | **Partial** — correction applied, generalization test excluded | Can include in batch correction to improve the correction model; excluded from classifier re-evaluation |
| Phase 6: Manuscript | **YES** — documented as limitation | Report group-level statistics (n=18 preclinAD, n=60 CN) with explicit note that per-sample labels were not available for classifier training |

---

## Manuscript language (draft, for Methods/Limitations)

> The kbase2022 cohort (Kim et al. 2022; PRJEB50447; n=78) was included in compositional analyses but excluded from supervised classifier training and evaluation. Per-sample amyloid-PET staging labels are not deposited in any public archive (ENA, NCBI SRA, or linked BioSample records), consistent with IRB restrictions on individual-level clinical data. Group-level counts (18 amyloid-positive preclinical AD, 60 amyloid-negative cognitively normal controls) are reported in the original publication. This cohort is included in PERMANOVA-based variance decomposition (Phase 4) where cohort identity alone is sufficient.

---

## Do not attempt to impute labels

Do not attempt to infer per-sample diagnosis from sequencing metadata, read counts, or microbiome composition. Any such inference would constitute label leakage and would invalidate the cross-cohort generalization analysis. The 18/60 split is known at the group level; individual assignments are not.
