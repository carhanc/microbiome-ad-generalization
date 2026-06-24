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
        "(**A**) Overview of the five 16S rRNA amplicon sequencing cohorts assembled for this study, "
        "including country of origin, sample size, sequencing mode, 16S variable region, diagnosis "
        "label availability, and analytic phases in which each cohort participates. "
        "(**B**) Summary of the six-phase analysis pipeline, from DADA2 processing through "
        "SHAP feature importance analysis."
    ),
    (
        "**Table 2.** Within-cohort nested cross-validation performance",
        2,
        "**Figure 2.** Within-cohort classification performance under nested cross-validation "
        "(10-fold outer / 5-fold inner CV). Bars show out-of-fold AUC-ROC; error bars represent "
        "95% bootstrap confidence intervals (1,000 resamples). Dashed line = chance (AUC = 0.50). "
        "LogReg = logistic regression with elastic-net regularization; LGBM = LightGBM."
    ),
    (
        "### 3.4 PERMANOVA Variance Decomposition",
        3,
        "**Figure 3.** Cross-cohort generalization failure. "
        "(**A**) Within-cohort (solid bars) versus LOCO cross-cohort (hatched bars) AUC-ROC for "
        "logistic regression (blue) and LightGBM (pink). Error bars = 95% bootstrap CI. "
        "(**B**, **C**) Pairwise single-cohort-train / single-cohort-test AUC heatmaps for "
        "logistic regression and LightGBM respectively. Diagonal (within-cohort) values are "
        "omitted; colour scale: red = near chance, green = high AUC."
    ),
    (
        "### 3.5 Batch Correction Does Not Recover Cross-Cohort Generalization",
        4,
        "**Figure 4.** PERMANOVA variance decomposition and beta-diversity structure. "
        "(**A**) Stacked bar showing the proportion of total Aitchison variance explained by "
        "cohort (R²=0.172, marginal) versus residual in the four labeled cohorts. "
        "(**B**) PCoA of Aitchison distance matrix (all five cohorts; n=509). Cohort R²=0.193 "
        "(F=30.13, p<0.0001, 9,999 permutations). Ellipses = 95% normal confidence ellipses. "
        "(**C**) Beta-dispersion: distance from each sample to its cohort centroid "
        "(F=13.93, p<0.0001, 9,999 permutations)."
    ),
    (
        "### 3.6 SHAP Analysis Reveals Cohort-Specific and Contradictory Taxonomic Signatures",
        5,
        "**Figure 5.** Effect of batch correction on LOCO cross-cohort AUC-ROC. "
        "(**A**) Logistic regression and (**B**) LightGBM LOCO AUC for each held-out cohort "
        "under three conditions: uncorrected (green), ComBat-seq corrected (magenta), and "
        "MMUPHin corrected (salmon). Error bars = 95% bootstrap CI. Dashed line = chance (0.50). "
        "Both methods failed to improve cross-cohort generalization; ComBat-seq substantially worsened it."
    ),
    (
        "---\n\n## 4. Discussion",
        6,
        "**Figure 6.** SHAP feature importance analysis. "
        "(**A**) Top-20 SHAP genera by maximum |SHAP| across cohorts for logistic regression. "
        "Dot size ∝ mean |SHAP|; red = AD-associated; blue = CN-associated; "
        "★ = directional flip taxon; yellow rows highlight flip taxa. "
        "(**B**) Directional flip taxa (logistic regression): genera appearing in the top-20 for "
        "≥2 cohorts with opposite AD/CN associations. Letters indicate cohort initials. "
        "(**C**) Pairwise Jaccard similarity of top-20 SHAP taxa between cohort pairs "
        "(logistic regression); mean off-diagonal Jaccard = 0.101."
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
    blockquote { margin-left: 2em; color: #555; }
    @media print {
        body { margin: 15mm; }
        .figure-block { page-break-inside: avoid; }
    }
""")

def build_figure_block(fig_num: int, caption: str) -> str:
    img_path = f"{FIGS}/manuscript_fig{fig_num}.png"
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
        "--css", css_path,
        "--metadata", "title=Microbiome AD Generalization",
        "--toc",
        "--toc-depth=2",
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
