"""Apply and validate one running footer on the assembled supplement."""

import re

import fitz


TABLE_FOOTER = re.compile(r"Table\s+S\d+[a-z]?\s*(?:---|—)\s*Page\s+\d+")
PAGE_MARKER = re.compile(r"(?:—|-{2,3})\s*Page")
FOOTER_PAGE = re.compile(r"\bPage\s+\d+")
NUMBER = re.compile(r"\d+")


def footer_lines(page):
    """Return text lines in the bottom 45 bp, in displayed coordinates."""
    for block in page.get_text("dict")["blocks"]:
        for line in block.get("lines", []):
            value = "".join(span["text"] for span in line["spans"]).strip()
            rect = fitz.Rect(line["bbox"])
            shown = rect * page.rotation_matrix
            if value and shown.y0 >= page.rect.height - 45:
                yield value, rect, shown


def number_pages(pdf):
    for index, page in enumerate(pdf, 1):
        for block in page.get_text("dict")["blocks"]:
            for line in block.get("lines", []):
                value = "".join(span["text"] for span in line["spans"]).strip()
                rect = fitz.Rect(line["bbox"])
                shown = rect * page.rotation_matrix
                old_number = value.isdecimal() and shown.y0 >= page.rect.height - 45 and abs(shown.x0 + shown.x1 - page.rect.width) < 50
                if TABLE_FOOTER.fullmatch(value) or old_number:
                    page.add_redact_annot(rect + (-1, -1, 1, 1), fill=(1, 1, 1))
        page.apply_redactions(images=0, graphics=0)
        shown = fitz.Point(page.rect.width / 2, page.rect.height - 18)
        point = shown * page.derotation_matrix
        page.insert_text(point, str(index), fontsize=11, fontname="tiro", rotate=(-page.rotation) % 360)


def check_numbers(pdf):
    for index, page in enumerate(pdf, 1):
        text = page.get_text()
        if TABLE_FOOTER.search(text):
            raise ValueError(f"page {index}: per-table footer remains")
        if PAGE_MARKER.search(text):
            raise ValueError(f"page {index}: per-page marker remains")
        footer = " ".join(value for value, _, _ in footer_lines(page))
        if FOOTER_PAGE.search(footer):
            raise ValueError(f"page {index}: per-page marker in footer")
        # Every numeric token anywhere in the footer region counts, so a
        # correct centred number cannot hide an extra one elsewhere.
        numbers = NUMBER.findall(footer)
        if numbers != [str(index)]:
            raise ValueError(f"page {index}: expected one running footer {index}, found {numbers}")
