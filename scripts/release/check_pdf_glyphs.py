"""Reject missing font glyphs and replacement characters in rendered prose PDFs."""

from pathlib import Path

import fitz


def check_pdf_glyphs(log: Path, pdf: Path) -> None:
    if "Missing character" in log.read_text(encoding="utf-8", errors="replace"):
        raise RuntimeError(f"missing glyph reported in {log}")
    with fitz.open(pdf) as document:
        extracted = "".join(page.get_text() for page in document)
    if "\ufffd" in extracted:
        raise RuntimeError(f"replacement character in {pdf}")
