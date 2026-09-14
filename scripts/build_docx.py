#!/usr/bin/env python3
"""Builds the Frontiers-submission manuscript DOCX from draft.md.

Frontiers requires the Word source file to carry Tables 1-4 as editable
Word tables and Figure 1-6 captions, consolidated at the end of the
document, with main figures NOT embedded (submitted as separate image
files). draft.md itself keeps Tables 1-4 interleaved in Results (needed
for the PDF build, manuscript_draft.pdf, via build_pdf.py, which is left
untouched by this script) and does not carry a figure-caption list of its
own (build_pdf.py's FIGURES list is the single source of truth for exact
caption text, reused here unchanged).

This script does not alter draft.md. It produces a DOCX-only derived
markdown: extracts the four Table N blocks (caption + table + any
immediately-following footnote paragraph) out of their in-text Results
position, and appends a "Tables" section (the four extracted blocks) and
a "Figure Captions" section (Figure 1-6 captions, from build_pdf.py's
FIGURES list) at the end of the document, after References. No prose,
number, citation, or table content is changed -- only position.
"""

import importlib.util
import re
import subprocess
import sys

BASE = "/Users/arhan/Desktop/microbiome-ad-generalization"
DRAFT = f"{BASE}/manuscript/draft.md"
OUT_MD = f"{BASE}/manuscript/_draft_docx_source.md"
OUT_DOCX = f"{BASE}/manuscript/draft.docx"

TABLE_CAPTIONS = [
    "**Table 1.**",
    "**Table 2.**",
    "**Table 3.**",
    "**Table 4.**",
]


def load_build_pdf_module():
    spec = importlib.util.spec_from_file_location("build_pdf", f"{BASE}/scripts/build_pdf.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def extract_table_blocks(text: str):
    """Extract each **Table N.** caption + markdown table + trailing
    footnote paragraph (if present, i.e. non-blank text before the next
    blank-line-delimited block boundary), in document order. Returns
    (text_with_tables_removed, [block_text, ...])."""
    lines = text.split("\n")
    blocks = []
    out_lines = []
    i = 0
    n = len(lines)
    while i < n:
        line = lines[i]
        if any(line.startswith(c) for c in TABLE_CAPTIONS):
            start = i
            j = i + 1
            # blank line after caption
            while j < n and lines[j].strip() == "":
                j += 1
            # markdown table rows
            while j < n and lines[j].strip().startswith("|"):
                j += 1
            # optional single blank line + one trailing footnote paragraph
            # (ends at the next blank line, i.e. the paragraph boundary).
            # A genuine table footnote is plain or italic-led text (e.g.
            # "*Amyloid-PET...", "N=60 for..."); a paragraph starting with
            # "**" is a bold subsection lead-in (this paper's convention for
            # substantive Results prose, e.g. "**Sensitivity analyses**...")
            # and must NOT be swallowed as a footnote (Table 4 has no
            # footnote at all -- the next paragraph is real Results text).
            if j < n and lines[j].strip() == "":
                k = j + 1
                if (
                    k < n
                    and lines[k].strip() != ""
                    and not lines[k].startswith("#")
                    and not lines[k].startswith("**")
                ):
                    while k < n and lines[k].strip() != "":
                        k += 1
                    j = k
            block = "\n".join(lines[start:j]).strip("\n")
            blocks.append(block)
            i = j
            # skip the blank separator line(s) left behind
            while i < n and lines[i].strip() == "":
                i += 1
            continue
        out_lines.append(line)
        i += 1
    return "\n".join(out_lines), blocks


def build_figure_captions_section(figures) -> str:
    parts = ["## Figure Captions", ""]
    for _marker, _fig_num, caption in figures:
        parts.append(caption)
        parts.append("")
    return "\n".join(parts).rstrip() + "\n"


def build_tables_section(blocks) -> str:
    parts = ["## Tables", ""]
    for block in blocks:
        parts.append(block)
        parts.append("")
    return "\n".join(parts).rstrip() + "\n"


def run(cmd, desc=""):
    print(f"  {desc or ' '.join(cmd[:3])}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  STDERR: {result.stderr[:800]}")
        sys.exit(1)
    return result


def main():
    print("\nBuilding manuscript DOCX (Frontiers source-file layout)\n")

    bp = load_build_pdf_module()

    with open(DRAFT, encoding="utf-8") as f:
        text = f.read()

    print("Extracting Table 1-4 blocks out of Results...")
    text_no_tables, table_blocks = extract_table_blocks(text)
    if len(table_blocks) != 4:
        print(f"  expected 4 table blocks, found {len(table_blocks)}")
        sys.exit(1)
    for blk in table_blocks:
        cap = blk.split("\n", 1)[0]
        print(f"  extracted: {cap}")

    # Insert Tables + Figure Captions sections after References, at the
    # very end of the document (before the trailing "---" separator, if
    # any survives at EOF).
    tables_section = build_tables_section(table_blocks)
    figcap_section = build_figure_captions_section(bp.FIGURES)

    final_text = text_no_tables.rstrip("\n") + "\n\n---\n\n" + tables_section + "\n---\n\n" + figcap_section + "\n"

    with open(OUT_MD, "w", encoding="utf-8") as f:
        f.write(final_text)
    print(f"  written {OUT_MD}")

    print("\nConverting to DOCX (pandoc, no figure embedding, matching Frontiers separate-figure convention)...")
    pandoc_cmd = [
        "pandoc",
        OUT_MD,
        "--from", "markdown+raw_html",
        "--to", "docx",
        "--standalone",
        "--mathml",
        "-o", OUT_DOCX,
    ]
    run(pandoc_cmd, "pandoc markdown -> docx")
    print(f"  written {OUT_DOCX}")

    print("\nApplying continuous line numbers + single PAGE-field footer...")
    run([sys.executable, f"{BASE}/scripts/add_docx_line_and_page_numbers.py"], "line/page numbering")

    import os
    try:
        os.remove(OUT_MD)
    except FileNotFoundError:
        pass

    print("\nDone.")


if __name__ == "__main__":
    main()
