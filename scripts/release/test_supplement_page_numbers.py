#!/usr/bin/env python3
"""Small regression fixtures for the assembled supplement footer gate."""

import fitz

from supplement_page_numbers import check_numbers, number_pages


def must_fail(pdf, label):
    try:
        check_numbers(pdf)
    except ValueError:
        return
    raise AssertionError(f"{label} passed the gate")


def numbered_page(extra=None):
    """One correctly numbered portrait page, optionally with extra footer text."""
    pdf = fitz.open()
    pdf.new_page()
    number_pages(pdf)
    check_numbers(pdf)
    if extra:
        pdf[0].insert_text((50, 826), extra, fontsize=11)
    return pdf


# A per-table footer must fail; numbering the page then passes (portrait).
pdf = fitz.open()
page = pdf.new_page()
page.insert_text((220, 826), "Table S2 --- Page 1", fontsize=11)
must_fail(pdf, "per-table footer")
number_pages(pdf)
check_numbers(pdf)

# Letter-suffixed split-table footers (Table S2a/S2b) must fail, then number cleanly.
for label in ("Table S2a — Page 1", "Table S2b — Page 2", "Table S2a --- Page 3"):
    split = fitz.open()
    split.new_page().insert_htmlbox(fitz.Rect(200, 815, 500, 840), label, css="* {font-size: 11px}")
    must_fail(split, label)
    number_pages(split)
    check_numbers(split)

# A correct centred 1 plus any extra number or Page marker in the footer fails.
for extra in ("99", "99 100", "Page 99", "--- Page 99", "— Page 99"):
    must_fail(numbered_page(extra), f"centred 1 plus {extra!r} at x=50")

# Landscape: the running number is upright and passes.
landscape = pdf.new_page(width=842, height=595)
number_pages(pdf)
check_numbers(pdf)
footer = [line for block in landscape.get_text("dict")["blocks"]
          for line in block.get("lines", []) if line["spans"][0]["text"] == "2"]
assert len(footer) == 1 and footer[0]["dir"] == (1.0, 0.0)

# Scientific integers in the body above the footer region are tolerated.
body = fitz.open()
page = body.new_page()
page.insert_text((72, 700), "GH5_2 n = 99, 100 sites; 37 genomes", fontsize=11)
page.insert_text((72, 780), "12 34", fontsize=11)
number_pages(body)
check_numbers(body)
print("PAGE_NUMBER_FIXTURE_PASS")
