#!/usr/bin/env python3
"""Post-processing pass: adds continuous line numbers (left margin) and page
numbers (bottom-center) to an already-built manuscript PDF.

This does not touch the HTML/CSS/JS build pipeline in build_pdf.py at all --
it operates purely on the final rendered PDF via PyMuPDF, reading existing
text-line geometry and stamping small annotations into the page margins.
No scientific content, figure, or table value is read, computed, or changed.

Body-text lines are identified by a left-edge heuristic: regular paragraph
and heading text in this document consistently starts at the same x0
(~70pt, the CSS body margin), while table cell contents and multi-column
figure-internal text start at other x-offsets. Restricting numbering to
lines at the standard margin gives every prose paragraph, heading, and
caption a number while naturally leaving table cell contents unnumbered.
"""

import os
import sys
import pymupdf as fitz

IN_PDF = "/Users/arhan/Desktop/microbiome-ad-generalization/manuscript/manuscript_draft.pdf"
TMP_PDF = IN_PDF + ".tmp"
OUT_PDF = IN_PDF  # final destination; written via a temp file, then renamed into place

MARGIN_X0_TARGET = 70.0
MARGIN_X0_TOL = 1.5
LINE_NUM_X = 40.0
LINE_NUM_FONTSIZE = 6.5
LINE_NUM_COLOR = (0.55, 0.55, 0.55)

PAGE_NUM_FONTSIZE = 8
PAGE_NUM_COLOR = (0.2, 0.2, 0.2)


def main():
    doc = fitz.open(IN_PDF)
    running_line = 0

    for pno in range(doc.page_count):
        page = doc[pno]
        blocks = page.get_text("dict")["blocks"]

        # Collect qualifying lines with their vertical position, then sort
        # top-to-bottom to preserve reading order.
        candidate_lines = []
        for b in blocks:
            if "lines" not in b:
                continue
            for l in b["lines"]:
                x0, y0, x1, y1 = l["bbox"]
                if abs(x0 - MARGIN_X0_TARGET) <= MARGIN_X0_TOL:
                    candidate_lines.append((y0, y1))
        candidate_lines.sort(key=lambda t: t[0])

        for (y0, y1) in candidate_lines:
            running_line += 1
            baseline_y = y1 - 1.5  # align near the visual text baseline
            page.insert_text(
                (LINE_NUM_X, baseline_y),
                str(running_line),
                fontsize=LINE_NUM_FONTSIZE,
                color=LINE_NUM_COLOR,
                fontname="helv",
            )

        # Page number, bottom-center, well within the existing 15mm/~42.5pt
        # bottom margin (page content on this document ends well above that).
        page_rect = page.rect
        label = str(pno + 1)
        text_width = fitz.get_text_length(label, fontname="helv", fontsize=PAGE_NUM_FONTSIZE)
        x = (page_rect.width - text_width) / 2
        y = page_rect.height - 22
        page.insert_text(
            (x, y),
            label,
            fontsize=PAGE_NUM_FONTSIZE,
            color=PAGE_NUM_COLOR,
            fontname="helv",
        )

    page_count = doc.page_count
    doc.save(TMP_PDF)
    doc.close()
    os.replace(TMP_PDF, OUT_PDF)
    print(f"Added {running_line} continuous line numbers across {page_count} pages.")
    print(f"Saved: {OUT_PDF}")


if __name__ == "__main__":
    main()
