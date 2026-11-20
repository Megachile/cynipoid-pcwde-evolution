"""Focused regression for superscript rendering and missing-glyph rejection."""

import tempfile
from pathlib import Path

import fitz

from build_manuscript import tex_escape
from check_pdf_glyphs import check_pdf_glyphs


def main() -> None:
    escaped = tex_escape("P < 10⁻¹⁰; 3 × 10⁻⁵")
    assert r"10\textsuperscript{\ensuremath{-10}}" in escaped, escaped
    assert r"10\textsuperscript{\ensuremath{-5}}" in escaped, escaped
    with tempfile.TemporaryDirectory() as root:
        log = Path(root) / "xelatex.log"
        log.write_text("Missing character: There is no ⁻ in font", encoding="utf-8")
        try:
            check_pdf_glyphs(log, Path(root) / "unused.pdf")
        except RuntimeError as exc:
            assert "missing glyph" in str(exc)
        else:
            raise AssertionError("missing-glyph log was accepted")
        log.write_text("XeLaTeX completed without missing glyphs\n", encoding="utf-8")
        pdf = Path(root) / "positive.pdf"
        document = fitz.open()
        page = document.new_page()
        page.insert_text((72, 72), "P < 10")
        document.save(pdf)
        document.close()
        check_pdf_glyphs(log, pdf)
    print("PDF_GLYPH_TEST_PASS")


if __name__ == "__main__":
    main()
