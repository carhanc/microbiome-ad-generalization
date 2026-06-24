# Dataset Registry — Gut Microbiome / Alzheimer's Disease / MCI Cross-Cohort Project

**Phase 0 output. Written before any modeling code. Date: 2026-06-20.**
**Phase 1 accession verification amendment: 2026-06-20. See amendment notes inline.**

This file documents every cohort investigated during data reconnaissance and Phase 1 verification. It is the authoritative record of what was found, what accession numbers exist, and what the data-access situation is for each study.

---

## Final Active Cohort Set (locked after Phase 1 verification)

**5 cohorts | 599 samples | 3 geographic regions: China (×3), South Korea, Kazakhstan**

| # | Cohort | Country | N (public SRA) | Accession | Verification |
|---|--------|---------|----------------|-----------|--------------|
| 1 | Zhuang 2018 | China | 86 (43 AD / 43 CN) | PRJNA554111 | ✓ 86 SRA runs confirmed |
| 2 | Ling 2020 | China | 171 (100 AD / 71 CN) | PRJNA633959 | ✓ 171 SRA runs confirmed |
| 3 | Shanghai 2022 | China | 180 of 302 participants | PRJNA489760 | ✓ 180 SRA runs confirmed (partial — see notes) |
| 4 | Kim/KBASE 2022 | South Korea | 78 (18 preclinAD / 60 CN) | PRJEB50447 | ✓ 78 ENA samples, V3-V4 confirmed |
| 5 | Kazakhstan 2022 | Kazakhstan | 84 (41 AD / 43 CN) | PRJNA811324 | ✓ 84 SRA runs confirmed |

**Scope:** 16S rRNA amplicon only. Shotgun arm (AlzBiom/Germany, Ferreiro/USA, Taiwan 2025) is out of scope for v1 — see CLAUDE.md Section 0 for rationale.

**Two originally-planned cohorts were excluded after Phase 1 accession verification:**

- **Liu 2019 (EXCLUDED — no public accession):** PRJNA496408 was listed in Phase 0 reconnaissance as the Liu 2019 deposit, but live verification confirmed it is a Zhejiang University *Mus musculus* mouse study, completely unrelated to Liu 2019 (PMID 31063846). An NCBI elink query from PMID 31063846 returns zero SRA records. Raw reads for Liu 2019 are not publicly deposited. Do not use PRJNA496408 for any purpose in this project.

- **ALBION 2025 / Greece (EXCLUDED — data not yet public; possible future add):** Accession PRJNA1297934 is confirmed correct — the paper at PMC12473069 explicitly cites it, and the study enrolled 99 samples (50 MCI / 49 CN, 16S V3-V4). However, NCBI BioProject reports PRJNA1297934 as "not public" as of 2026-06-20 verification. Data has been deposited but not released. Excluded for now. **If the data is released before manuscript submission, ALBION should be revisited as a 6th cohort** — it is the only European 16S cohort identified and would materially improve geographic diversity. Check PRJNA1297934 periodically.

**curatedMetagenomicData (Bioconductor):** Contains no Alzheimer's/MCI cohorts. Not a viable data source for this project.

---

## Active Cohort Entries

---

### Cohort 1 — Zhuang 2018 (China — Chongqing)

| Field | Value |
|---|---|
| **First author / year** | Zhuang et al. 2018 |
| **Full citation** | Zhuang Z-Q et al. "Gut Microbiota is Altered in Patients with Alzheimer's Disease." *Journal of Alzheimer's Disease* 63(4):1337–1346, 2018. |
| **PMID** | 29758946 |
| **DOI** | 10.3233/JAD-180176 |
| **Accession** | **PRJNA554111** (NCBI SRA) |
| **Verification** | ✓ 86 SRA runs confirmed via NCBI E-utilities on 2026-06-20 |
| **Country** | China (Chongqing; Third Military Medical University) |
| **N cases** | 43 AD |
| **N controls** | 43 cognitively normal (age- and gender-matched) |
| **N MCI** | 0 (AD vs. CN only) |
| **Sequencing type** | 16S rRNA amplicon |
| **16S region** | Not explicitly stated in abstract; verify in SRA record Methods — likely V3-V4 |
| **Platform** | Illumina (confirm from SRA record on download) |
| **Covariates available** | Age, sex (matched); clinical AD diagnosis confirmed; APOE4 not confirmed |
| **Processed tables in supplement** | No |
| **Data access** | Public; `fasterq-dump PRJNA554111` |
| **Notes** | This is one of the two Chinese cohorts at the center of the Bacteroides directional-flip finding in the 2023 meta-analysis. A 2025 reanalysis (PMC12439993) reprocessed these exact reads with DADA2 + SILVA 138.2 — use as a benchmark when verifying our own DADA2 output. |

---

### Cohort 2 — Ling 2020 (China — Lishui, Zhejiang)

| Field | Value |
|---|---|
| **First author / year** | Ling et al. 2020 |
| **Full citation** | Ling Z et al. "Structural and Functional Dysbiosis of Fecal Microbiota in Chinese Patients with Alzheimer's Disease." *Frontiers in Cell and Developmental Biology* 8:634069, 2021. |
| **PMID** | 33604329 |
| **PMC** | PMC7889981 |
| **DOI** | 10.3389/fcell.2020.634069 |
| **Accession** | **PRJNA633959** (NCBI SRA; SRP identifier: SRP262626) |
| **Verification** | ✓ 171 SRA runs confirmed via NCBI E-utilities on 2026-06-20 |
| **Country** | China (Lishui, Zhejiang province; Zhejiang University affiliated) |
| **N cases** | 100 AD |
| **N controls** | 71 cognitively normal (age- and gender-matched) |
| **N MCI** | 0 |
| **Sequencing type** | 16S rRNA amplicon |
| **16S region** | V3–V4 |
| **Platform** | Illumina MiSeq (300×2 bp paired-end) |
| **Covariates available** | Age, sex (matched), BMI, smoking, alcohol, MMSE, WAIS, Barthel index; comorbidities (hypertension, diabetes, hypercholesterolemia). APOE4 not confirmed. |
| **Processed tables in supplement** | No downloadable OTU/abundance table; supplementary figures only |
| **Data access** | Public; `fasterq-dump PRJNA633959` |
| **Notes** | Largest Chinese clinical AD cohort (N=171) among the confirmed public 16S datasets. Good covariate coverage including BMI and comorbidities — useful for the PERMANOVA variance decomposition. |

---

### Cohort 3 — Shanghai Aging / Memory Study 2022 (China — Shanghai)

| Field | Value |
|---|---|
| **First author / year** | Gao et al. 2022 |
| **Full citation** | Gao Y et al. "Altered Gut Microbiota and Its Clinical Relevance in Mild Cognitive Impairment and Alzheimer's Disease: Shanghai Aging Study and Shanghai Memory Study." *Nutrients* 14(19):3959, 2022. |
| **PMID** | 36235612 |
| **PMC** | PMC9570603 |
| **DOI** | 10.3390/nu14193959 |
| **Accession** | **PRJNA489760** (NCBI SRA) |
| **Verification** | ✓ 180 SRA runs confirmed on 2026-06-20; BioProject title matches ("Fecal and blood microbiota of AD, MCI and normal controls") |
| **Country** | China (Shanghai; Ruijin Hospital, Shanghai Jiao Tong University) |
| **N AD (paper)** | 83 |
| **N MCI (paper)** | 125 |
| **N CN (paper)** | 94 |
| **N total (paper)** | 302 |
| **N in public SRA** | **180** (122 participants' data not deposited or in restricted access) |
| **Sequencing type** | 16S rRNA amplicon |
| **16S region** | V3–V4 |
| **Platform** | Illumina MiSeq PE300 |
| **Covariates available** | Age, sex, education, APOE genotype, hypertension, diabetes, stroke history |
| **Processed tables in supplement** | "Data availability: Not applicable" — raw reads in SRA only |
| **Data access** | Public; `fasterq-dump PRJNA489760` |
| **Notes** | **IMPORTANT: Only 180 of 302 participants are in the public SRA deposit.** The paper-reported sample sizes (83 AD / 125 MCI / 94 CN) are larger than what is accessible. The breakdown of the 180 public samples by diagnosis group is not specified in the SRA metadata — must be determined from sample-level metadata when downloading. Treat n=180 (unknown AD/MCI/CN split) as the working sample size; report the 302 figure from the paper in the manuscript with a note on restricted-access samples. This is the only three-group dataset (CN/MCI/AD) in the active cohort set. |

---

### Cohort 4 — Kim et al. 2022 (South Korea — Seoul, KBASE study)

| Field | Value |
|---|---|
| **First author / year** | Kim et al. 2022 (KBASE — Korean Brain Aging Study for Early Diagnosis and Prediction of AD) |
| **Full citation** | Kim M et al. "Gut microbiome alterations in preclinical Alzheimer's disease." *PLOS ONE* 17(12):e0278276, 2022. |
| **PMID** | 36525397 |
| **PMC** | PMC9707757 |
| **DOI** | 10.1371/journal.pone.0278276 |
| **Accession** | **PRJEB50447** (EBI-ENA) |
| **Verification** | ✓ 78 ENA samples confirmed; V3-V4; study title "Gut microbiome alterations in preclinical Alzheimer's disease" confirmed on 2026-06-20 |
| **Country** | South Korea (Seoul National University Hospital; KBASE cohort) |
| **N preclinical AD** | 18 (cognitively normal, Aβ+ by amyloid PET) |
| **N controls** | 60 (cognitively normal, Aβ−) |
| **Sequencing type** | 16S rRNA amplicon |
| **16S region** | V3–V4 (confirmed in ENA metadata) |
| **Platform** | Illumina MiSeq (2×300 bp paired-end) |
| **Covariates available** | Age, sex, BMI, APOE4 status, diabetes, hypertension |
| **Processed tables in supplement** | Clinical metadata restricted (IRB); raw reads public |
| **Data access** | Public reads via EBI-ENA PRJEB50447: `enaGroupGet -g read -f fastq PRJEB50447` |
| **Notes** | **LABEL CAVEAT:** This cohort uses biomarker-confirmed preclinical AD (amyloid PET-positive but cognitively normal), not a clinical AD or MCI diagnosis. This is a fundamentally different phenotype from the other four cohorts (which use clinical diagnosis). Must be flagged explicitly anywhere results from this cohort are pooled with or compared to clinically-diagnosed cohorts. Consider analyzing as a distinct group or clearly labeling as "preclinical AD" throughout. The small case N (n=18) limits this cohort's utility for training; it is most valuable as a test set. |

---

### Cohort 5 — Kazakhstan AD Study 2022

| Field | Value |
|---|---|
| **First author / year** | Auyezbayeva et al. 2022 |
| **Full citation** | Auyezbayeva A et al. "Study of gut microbiota alterations in Alzheimer's dementia patients from Kazakhstan." *Frontiers in Aging Neuroscience* 14:938672, 2022. |
| **PMC** | PMC9448737 |
| **DOI** | 10.3389/fnagi.2022.938672 |
| **Accession** | **PRJNA811324** (NCBI SRA) |
| **Verification** | ✓ 84 SRA runs confirmed via NCBI E-utilities on 2026-06-20 |
| **Country** | Kazakhstan (Nur-Sultan/Astana; Center for Life Sciences, National Laboratory Astana) |
| **N cases** | 41 AD |
| **N controls** | 43 cognitively healthy |
| **Sequencing type** | 16S rRNA amplicon |
| **16S region** | Not confirmed in abstract — verify from SRA record on download |
| **Platform** | Illumina (sequenced at Novogene, Beijing) |
| **Covariates available** | Age, sex, ethnicity, MMSE, CDT, APOE4, diabetes, hypertension, biochemical labs (adiponectin, glucose, lipids, CRP) |
| **Processed tables** | Supplementary PowerPoint with tables (journal page) |
| **Data access** | Public; `fasterq-dump PRJNA811324` |
| **Notes** | Geographically unique — the only Central Asian cohort in the active set. Provides non-East-Asian, non-Western comparison. Rich clinical covariates including APOE4 and biochemical markers are valuable for the PERMANOVA covariate partitioning. AD vs. CN only (no MCI arm). **SEQUENCING: Single-end reads confirmed** — fasterq-dump produced only `_1.fastq` files (no `_2.fastq`), 84 files at 8.0 GB total. DADA2 must use the single-end workflow for this cohort; all other 4 cohorts are paired-end. This must be handled explicitly in `01_harmonize_taxonomy.py`. |

---

## Excluded Cohorts

### Excluded: Liu 2019 (China — Zhejiang)

| Field | Value |
|---|---|
| **Citation** | Liu P et al. "Altered microbiomes distinguish Alzheimer's disease from amnestic mild cognitive impairment and health in a Chinese cohort." *Brain, Behavior, and Immunity* 80:633–643, 2019. PMID 31063846. |
| **Reason for exclusion** | **No public raw sequencing data.** The accession PRJNA496408 identified during Phase 0 reconnaissance was verified on 2026-06-20 to be an unrelated Zhejiang University *Mus musculus* mouse 16S study (117 runs, mouse organism confirmed). An NCBI elink query from PMID 31063846 returns zero SRA links. Raw reads for Liu 2019 are not deposited in NCBI SRA. Author-contact to obtain data was not pursued (out of scope for v1 timeline). |
| **Scientific value if data ever becomes available** | Three-group design (33 AD / 32 aMCI / 32 CN) with AUC 0.940 (AD vs. CN) and 0.925 (AD vs. aMCI) reported in-cohort — a compelling benchmark. From First Affiliated Hospital, Zhejiang University. Worth revisiting for a future multi-cohort extension. |

---

### Excluded: ALBION 2025 (Greece — Athens) — Possible Future Add

| Field | Value |
|---|---|
| **Citation** | Maraki M et al. "Gut Microbiome Alterations in Mild Cognitive Impairment: Findings from the ALBION Greek Cohort." *Journal of Alzheimer's Disease*, 2025. PMC12473069. |
| **Accession** | PRJNA1297934 (NCBI SRA) — correct accession confirmed from paper |
| **Reason for exclusion** | **Data deposited but not yet publicly released.** NCBI BioProject PRJNA1297934 is confirmed as the deposit site in the paper, but as of 2026-06-20 NCBI reports it as "not public." Author-contact for early release was not pursued. |
| **Action if data releases** | Add as Cohort 6. It is the only European 16S cohort identified, V3-V4, 99 samples (50 MCI / 49 CN). Note: within-cohort platform difference (MiSeq for 50 samples, NextSeq2000 for 49) must be handled as a covariate if included. Check PRJNA1297934 periodically before submission. |

---

## Out-of-Scope Cohorts (Shotgun — v1 exclusion)

The following shotgun cohorts were identified during Phase 0 but are excluded from v1 per the scope decision in CLAUDE.md Section 0. They are documented here for future reference only.

| Cohort | Country | N | Accession | Notes |
|--------|---------|---|-----------|-------|
| AlzBiom 2022 | Germany | 175 (75 AD / 100 CN) | PRJEB47976 (ENA) | Shotgun; GitHub: UliSchopp/AlzBiom |
| Ferreiro/Dantas 2023 | USA | 164 (49 preclinAD / 115 CN) | PRJNA798058 | Shotgun; WashU Knight ADRC |
| Taiwan 2025 | Taiwan | 439 (119 MCI / 320 CN) | PRJNA1258384 | Shotgun; processed tables on request |

---

## Studies Without Public Data (Reference Only)

| Study | Country | N | Reason Not Usable |
|---|---|---|---|
| Vogt 2017 | USA | 50 (25 AD / 25 CN) | Raw reads not in SRA; supplementary differential-abundance tables only (Sci Rep S2/S3) |
| Haran 2019 | USA | 108 | Data on request only (John.Haran@umassmemorial.org) |
| Nagpal 2019 | USA | 17 | Tiny N; no confirmed public accession |
| Guo 2021 | China | 56 | No confirmed SRA accession |
| Zhou 2021 | China | 92 | No confirmed SRA accession |
| Liu 2021 | China | 42 | No confirmed SRA accession |
| Sheng 2021 | China | 105 | No confirmed SRA accession |
| Zhang 2021 | China | 127 | No confirmed SRA accession |
| SILCODE (Beijing) | China | ~88 | On-request only |
| Li 2019 | China | 90 | No confirmed SRA accession |

---

## Data Download Instructions (Phase 1)

### NCBI SRA cohorts (Zhuang 2018, Ling 2020, Shanghai 2022, Kazakhstan 2022)

```bash
# Install SRA toolkit
conda install -c bioconda sra-tools

# Download each cohort — split-files gives R1/R2 for paired-end reads
fasterq-dump --split-files --outdir data/raw/zhuang2018/ PRJNA554111
fasterq-dump --split-files --outdir data/raw/ling2020/   PRJNA633959
fasterq-dump --split-files --outdir data/raw/shanghai2022/ PRJNA489760
fasterq-dump --split-files --outdir data/raw/kazakhstan2022/ PRJNA811324

# Parallel download with prefetch first (more reliable for large projects):
prefetch --option-file <(echo -e "PRJNA554111\nPRJNA633959\nPRJNA489760\nPRJNA811324")
```

### EBI-ENA cohort (Kim/KBASE 2022)

```bash
pip install enaBrowserTools
enaGroupGet -g read -f fastq -d data/raw/kbase2022/ PRJEB50447
```

---

## Key Technical Notes for Phase 1 Harmonization

1. **16S region:** Most confirmed cohorts used V3-V4 (Ling 2020, Shanghai 2022, Kim/KBASE confirmed). Zhuang 2018 and Kazakhstan 2022 regions must be confirmed from SRA sample records on download. Log the region as a per-cohort metadata field — it is a covariate in the PERMANOVA model (Phase 4).

2. **Processing pipeline:** DADA2 (v1.28+) with SILVA 138.2 taxonomy for all cohorts. The 2025 reanalysis (PMC12439993) used exactly this combination on Zhuang 2018 reads and is a direct benchmark for our DADA2 output on that cohort.

3. **Genus-level harmonization:** After DADA2, collapse ASVs to genus level. Retain genera present in at least 2 cohorts and in at least 20% of samples in at least one cohort (standard prevalence filter to remove rare sporadic taxa).

4. **Zero-inflation / CLR:** Use multiplicative replacement (`zCompositions::multRepl()`) for zero handling before CLR transformation. Do not use naive +1 pseudocount.

5. **Shanghai 2022 partial deposit:** The 180 public SRA samples are a subset of 302 enrolled participants. Determine the AD/MCI/CN breakdown of the 180 from SRA sample metadata on download. Report in manuscript methods that the publicly deposited subset was used.

6. **Kim/KBASE phenotype label:** The 18 "cases" are amyloid-PET positive but cognitively normal — not clinical AD/MCI. Use label "preclinAD" rather than "AD" throughout the pipeline for this cohort. Flag explicitly in any analysis that pools or compares them with clinically-diagnosed cohorts.

---

*Registry compiled 2026-06-20 (Phase 0). Amended 2026-06-20 (Phase 1 verification: Liu 2019 dropped, ALBION excluded pending data release, Shanghai 2022 partial deposit confirmed, final active set locked at 5 cohorts / 599 samples).*
