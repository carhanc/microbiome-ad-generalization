#!/usr/bin/env python3
"""Builds manuscript/supplementary.pdf from supplementary.md.

Unlike build_pdf.py, no figure-marker insertion is needed here: the actual
Figure S1-S5 images are embedded directly in supplementary.md via ordinary
markdown image syntax, immediately above each figure's "Source figure" line.
"""

import subprocess, sys, os, textwrap

BASE = "/Users/arhan/Desktop/microbiome-ad-generalization"
SUPP = f"{BASE}/manuscript/supplementary.md"
OUT_HTML = f"{BASE}/manuscript/_supplementary.html"
OUT_PDF = f"{BASE}/manuscript/supplementary.pdf"

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
        font-size: 8.5pt;
        margin: 1em 0;
    }
    th, td {
        border: 1px solid #aaa;
        padding: 3px 6px;
        text-align: left;
    }
    th { background: #e8e8e8; font-weight: bold; }
    tr:nth-child(even) { background: #f7f7f7; }
    img { max-width: 100%; height: auto; display: block; margin: 0.8em auto;
          page-break-inside: avoid; break-inside: avoid; }
    hr { border: none; border-top: 1px solid #ccc; margin: 1.5em 0; }
    #title-block-header { display: none; }
    blockquote { margin-left: 2em; color: #555; }
    @media print {
        body { margin: 15mm; }
    }
""")


def run(cmd, desc=""):
    print(f"  {desc or ' '.join(cmd[:3])}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  STDERR: {result.stderr[:800]}")
        sys.exit(1)
    return result


def main():
    print("\nBuilding supplementary PDF\n")

    css_path = f"{BASE}/manuscript/_supp_style.css"
    with open(css_path, "w") as f:
        f.write(CSS)

    print("Converting to HTML (pandoc)...")
    pandoc_cmd = [
        "pandoc",
        SUPP,
        "--from", "markdown+raw_html+implicit_figures",
        "--to", "html5",
        "--standalone",
        "--embed-resources",
        "--mathjax",
        "--css", css_path,
        "--metadata", "title=Supplementary Materials",
        "-o", OUT_HTML,
    ]
    run(pandoc_cmd, "pandoc markdown -> HTML")
    print(f"  written {OUT_HTML}")

    print("\nPrinting PDF (Chrome headless)...")
    chrome = "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
    chrome_cmd = [
        chrome,
        "--headless=new",
        "--disable-gpu",
        "--no-sandbox",
        "--run-all-compositor-stages-before-draw",
        "--virtual-time-budget=8000",
        "--print-to-pdf-no-header",
        f"--print-to-pdf={OUT_PDF}",
        "--no-pdf-header-footer",
        OUT_HTML,
    ]
    run(chrome_cmd, "Chrome headless PDF print")
    print(f"  written {OUT_PDF}")

    print("\nVerifying output...")
    if os.path.exists(OUT_PDF):
        size_kb = os.path.getsize(OUT_PDF) // 1024
        print(f"  {OUT_PDF} ({size_kb} KB)")
    else:
        print("  PDF not found")
        sys.exit(1)

    for tmp in [OUT_HTML, css_path]:
        try:
            os.remove(tmp)
        except FileNotFoundError:
            pass

    print("\nDone.")


if __name__ == "__main__":
    main()
