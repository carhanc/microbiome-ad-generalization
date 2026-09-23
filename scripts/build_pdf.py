#!/usr/bin/env python3

import subprocess, sys, os, shutil, textwrap

BASE   = "/Users/arhan/Desktop/microbiome-ad-generalization"
DRAFT  = f"{BASE}/manuscript/draft.md"
FIGS   = f"{BASE}/results/figures"
OUT_MD = f"{BASE}/manuscript/_draft_with_figs.md"
OUT_HTML = f"{BASE}/manuscript/_draft_with_figs.html"
OUT_PDF  = f"{BASE}/manuscript/manuscript_draft.pdf"

FIGURES = [
    (
        "**Table 1.** Cohort characteristics.",
        1,
        "**Figure 1.** Study design. "
        "(**A**) Overview of the five 16S rRNA amplicon sequencing cohorts assembled for this study: "
        "four with per-sample diagnosis labels supporting supervised classification (Phases 2–3, "
        "5–6) and a fifth (Kim/KBASE 2022) contributing only to unsupervised compositional analysis "
        "(Phase 4). Columns give country of origin, sample size, sequencing mode, 16S variable "
        "region, diagnosis label availability, and analytic phases in which each cohort participates. "
        "(**B**) Summary of the six-phase analysis pipeline, from DADA2 processing through "
        "SHAP feature-importance and logistic-regression coefficient-stability analyses. Phases 2–3 "
        "(within-cohort baseline and cross-cohort generalization) use strict training-only genus "
        "selection, re-derived independently within every split, as the primary classifier-"
        "performance pipeline (Section 2.5.1); Phases 4–6 use the fixed common 396-genus universe "
        "(Section 2.2), which those specific analyses require as a single, shared feature space "
        "across cohorts."
    ),
    (
        "**Table 2.** Within-cohort nested cross-validation performance",
        2,
        "**Figure 2.** Within-cohort classification performance under nested cross-validation "
        "(10-fold outer / 5-fold inner CV), using strict training-only genus selection as the "
        "primary classifier-performance pipeline (Section 2.5.1). Bars show out-of-fold AUC-ROC; "
        "error bars represent 95% CIs from the diagnosis-stratified, participant-level bootstrap "
        "(10,000 replicates; Section 2.5.2), matching Table 2. Dashed line = chance (AUC = 0.50). "
        "LogReg = L2-regularized logistic regression; LGBM = LightGBM. Matched values under the "
        "fixed common-universe (396-genus) schema are given in Supplementary Table S4."
    ),
    (
        "### 3.4 PERMANOVA Variance Decomposition",
        3,
        "**Figure 3.** Cross-cohort generalization limitations, using strict training-only genus "
        "selection as the primary classifier-performance pipeline (Section 2.5.1). "
        "(**A**) Within-cohort (solid bars) versus LOCO cross-cohort (hatched bars) AUC-ROC for "
        "logistic regression (blue) and LightGBM (pink). Error bars on the LOCO bars = 95% CI from "
        "the diagnosis-stratified, participant-level bootstrap (10,000 replicates; Section 2.5.2), "
        "matching Table 3; within-cohort bars are shown without error bars in this panel. "
        "(**B**, **C**) Pairwise single-cohort-train / single-cohort-test AUC heatmaps for "
        "logistic regression and LightGBM respectively. Diagonal (within-cohort) values are "
        "omitted; colour scale: red = lower AUC, green = higher AUC. Matched values under the fixed "
        "common-universe (396-genus) schema are given in Supplementary Table S1; one cell "
        "(Kazakhstan→Ling 2020, LightGBM) differs materially between the two pipelines (0.58 vs. "
        "0.44) and is discussed in Section 3.3."
    ),
    (
        "### 3.5 Label-Blind Transductive Batch Adjustment Does Not Improve Cross-Cohort Transferability",
        4,
        "**Figure 4.** PERMANOVA variance decomposition and beta-diversity structure. "
        "(**A**) Marginal R² values from a two-variable PERMANOVA model on the four labeled "
        "cohorts (n=401): cohort R²=0.172 and diagnosis R²=0.014, shown as independent bars. "
        "These Type III marginal R² values are estimated independently and do not sum to 1. "
        "(**B**) PCoA of Aitchison distance matrix (all five cohorts; n=509). Cohort R²=0.193 "
        "(F=30.13, p=0.0001, 9,999 permutations). Ellipses = 95% normal confidence ellipses. "
        "(**C**) Beta-dispersion: distance from each sample to its cohort centroid "
        "(F=13.93, p=0.0001, 9,999 permutations)."
    ),
    (
        "### 3.6 SHAP Analysis Reveals Cohort-Specific Taxonomic Signatures and a Descriptive Coefficient-Stability Screen",
        5,
        "**Figure 5.** Effect of label-blind batch correction on LOCO cross-cohort AUC-ROC "
        "(primary analysis; no diagnosis information used during correction). "
        "(**A**) Logistic regression and (**B**) LightGBM LOCO AUC for each held-out cohort "
        "under three conditions: uncorrected (green), ComBat-seq label-blind corrected (magenta), and "
        "MMUPHin label-blind corrected (salmon). Error bars = 95% bootstrap CI (1,000 non-stratified "
        "resamples on held-out test predictions; Section 2.5). Dashed line = chance (0.50). "
        "Under this label-blind design, both methods left mean LOCO AUC essentially unchanged "
        "relative to uncorrected (within 0.007 AUC for both methods and models), neither "
        "recovering nor substantially worsening cross-cohort generalization. An exploratory "
        "label-informed transductive sensitivity comparison, in which correction used each sample's own "
        "true diagnosis label and which produced substantially different, more severe "
        "degradation for some cohorts, is reported separately in Supplementary Figure S4 "
        "and is not treated as evidence about prospective performance (Section 2.7, Section 3.5)."
    ),
    (
        "---\n\n## 4. Discussion",
        6,
        "**Figure 6.** SHAP feature importance and coefficient-based direction analysis, logistic regression. "
        "(**A**) Top-20 SHAP genera by maximum |SHAP| across cohorts. Dot size ∝ mean |SHAP|; "
        "feature importance only — no direction is encoded in this panel. "
        "(**B**) Fitted-coefficient direction by cohort for the eight genera meeting a pre-specified "
        "descriptive coefficient-sign-stability screen (top-20 in ≥2 cohorts; a stable positive "
        "coefficient sign in ≥8/10 outer folds in ≥1 cohort and a stable negative sign in ≥8/10 folds "
        "in ≥1 other). Cell color = median fitted logistic-regression coefficient across the 10 outer "
        "folds; cell text = fold-count and direction. Bordered cells meet the ≥8/10 stability "
        "criterion. This is a descriptive screen, not a multiplicity-corrected hypothesis test. "
        "(**C**) Pairwise Jaccard similarity of top-20 SHAP taxa between cohort pairs; "
        "mean of the six unique pairwise values = 0.135, significantly above a null matched to this "
        "statistic's construction (100,000 replicates of four random top-20 sets; p<0.0001). "
        "A prevalence-restricted null sensitivity analysis (Supplementary Table S13) shows logistic "
        "regression remains significant while LightGBM does not."
    ),
]

CSS = textwrap.dedent("""\
    body {
        font-family: "Times New Roman", Times, serif;
        font-size: 11pt;
        line-height: 1.5;
        max-width: 170mm;
        margin: 20mm auto;
        color: #111;
    }
    h1 { font-size: 14pt; margin-top: 1.5em; }
    h2 { font-size: 12pt; margin-top: 1.4em; border-bottom: 1px solid #ccc; padding-bottom: 2px; }
    h3 { font-size: 11pt; margin-top: 1.2em; }
    p  { text-align: justify; margin: 0.5em 0; }
    table {
        border-collapse: collapse;
        width: 100%;
        font-size: 9pt;
        margin: 1em 0;
    }
    th, td {
        border: 1px solid #aaa;
        padding: 3px 6px;
        text-align: left;
    }
    th { background: #e8e8e8; font-weight: bold; }
    tr:nth-child(even) { background: #f7f7f7; }
    .figure-block {
        text-align: center;
        margin: 1.5em 0;
        page-break-inside: avoid;
    }
    .figure-block img {
        max-width: 100%;
        height: auto;
        display: block;
        margin: 0 auto 0.5em auto;
    }
    .figure-caption {
        font-size: 9pt;
        text-align: left;
        color: #333;
        margin-top: 0.3em;
    }
    hr { border: none; border-top: 1px solid #ccc; margin: 1.5em 0; }
    /* Pandoc's standalone template renders --metadata title=... as a visible
       H1 title block above the actual manuscript content. The manuscript
       must start directly with its own title/authors/Abstract, so this
       build artifact is hidden here rather than removing the metadata
       title itself (which still sets the invisible PDF/HTML document title). */
    #title-block-header { display: none; }
    blockquote { margin-left: 2em; color: #555; }
    @media print {
        body { margin: 15mm; }
        .figure-block { page-break-inside: avoid; }
    }
""")

def build_figure_block(fig_num: int, caption: str) -> str:
    img_path = f"{FIGS}/manuscript_fig{fig_num}.jpg"
    if not os.path.exists(img_path):
        print(f"  figure {fig_num} not found at {img_path}, skipping")
        return ""
    return (
        f'\n\n<div class="figure-block">\n'
        f'![Figure {fig_num}]({img_path}){{width=100%}}\n'
        f'<p class="figure-caption">{caption}</p>\n'
        f'</div>\n\n'
    )

def insert_figures(text: str) -> str:
    # markers are literal strings from draft.md — fragile but works
    for marker, fig_num, caption in FIGURES:
        block = build_figure_block(fig_num, caption)
        if not block:
            continue
        if marker not in text:
            print(f"  marker for figure {fig_num} not found in draft")
            continue
        text = text.replace(marker, block + marker, 1)
    return text

def run(cmd, desc=""):
    print(f"  {desc or ' '.join(cmd[:3])}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  STDERR: {result.stderr[:800]}")
        sys.exit(1)
    return result

def main():
    print("\nBuilding manuscript PDF\n")

    with open(DRAFT, encoding="utf-8") as f:
        text = f.read()

    print("Inserting figures…")
    text_with_figs = insert_figures(text)
    with open(OUT_MD, "w", encoding="utf-8") as f:
        f.write(text_with_figs)
    print(f"  written {OUT_MD}")

    css_path = f"{BASE}/manuscript/_style.css"
    with open(css_path, "w") as f:
        f.write(CSS)

    print("\nConverting to HTML (pandoc)…")
    pandoc_cmd = [
        "pandoc",
        OUT_MD,
        "--from", "markdown+raw_html+implicit_figures",
        "--to", "html5",
        "--standalone",
        "--embed-resources",
        "--mathjax",
        "--css", css_path,
        "--metadata", "title=Microbiome AD Generalization",
        "-o", OUT_HTML,
    ]
    run(pandoc_cmd, "pandoc markdown → HTML")
    print(f"  written {OUT_HTML}")

    print("\nPrinting PDF (Chrome headless)…")
    chrome = "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
    chrome_cmd = [
        chrome,
        "--headless=new",
        "--disable-gpu",
        "--no-sandbox",
        "--run-all-compositor-stages-before-draw",
        # --mathjax renders equations via async JS after page load; without
        # giving Chrome time to execute and finish typesetting before the
        # print snapshot, the PDF can capture raw/unrendered LaTeX.
        "--virtual-time-budget=8000",
        "--print-to-pdf-no-header",
        f"--print-to-pdf={OUT_PDF}",
        f"--no-pdf-header-footer",
        OUT_HTML,
    ]
    run(chrome_cmd, "Chrome headless PDF print")
    print(f"  written {OUT_PDF}")

    print("\nVerifying output…")
    if os.path.exists(OUT_PDF):
        size_kb = os.path.getsize(OUT_PDF) // 1024
        print(f"  {OUT_PDF} ({size_kb} KB)")
    else:
        print("  PDF not found")
        sys.exit(1)

    for tmp in [OUT_MD, OUT_HTML, css_path]:
        try:
            os.remove(tmp)
        except FileNotFoundError:
            pass

    print("\nDone.")

if __name__ == "__main__":
    main()
