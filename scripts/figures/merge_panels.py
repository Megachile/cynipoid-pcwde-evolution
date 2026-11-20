#!/usr/bin/env python3
"""Package validated single-page panels into a public multipage PDF."""
import argparse
from pathlib import Path

from pypdf import PdfReader, PdfWriter


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--panel", type=Path, action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--expected-pages", type=int, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise SystemExit(f"refusing to overwrite: {args.output}")
    writer = PdfWriter()
    if len(args.panel) != args.expected_pages:
        raise SystemExit("panel count does not match --expected-pages")
    for panel in args.panel:
        source = PdfReader(panel)
        if len(source.pages) != 1:
            raise SystemExit(f"{panel} is not exactly one page")
        writer.add_page(source.pages[0])
    with args.output.open("wb") as handle:
        writer.write(handle)
    if len(PdfReader(args.output).pages) != args.expected_pages:
        raise SystemExit("merged public figure has the wrong page count")


if __name__ == "__main__":
    main()
