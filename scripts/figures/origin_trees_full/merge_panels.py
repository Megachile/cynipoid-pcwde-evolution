#!/usr/bin/env python3
"""Permission-safe public wrapper for the established Figure S1 pdfpages composer."""

from __future__ import annotations

import argparse
import hashlib
import os
import shutil
import subprocess
from pathlib import Path

from pypdf import PdfReader, PdfWriter


ORDER = ("origin_trees_full_GH5_2_ggtree.pdf", "origin_trees_full_PL1_8_ggtree.pdf", "origin_trees_full_PL4_4_ggtree.pdf")
TEX = r"""\documentclass{article}\usepackage{pdfpages}\begin{document}
\includepdf[pages=-,fitpaper=true]{origin_trees_full_GH5_2_ggtree.pdf}
\includepdf[pages=-,fitpaper=true]{origin_trees_full_PL1_8_ggtree.pdf}
\includepdf[pages=-,fitpaper=true]{origin_trees_full_PL4_4_ggtree.pdf}
\end{document}
"""


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--panel", action="append", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--xelatex", default="xelatex")
    parser.add_argument("--panel-scale", type=float, default=0.75,
                        help="Uniform native-panel scale; complete S1 panels retain large readable pages")
    args = parser.parse_args()
    if not 0 < args.panel_scale <= 1:
        raise SystemExit("panel scale must be in (0,1]")
    names = tuple(path.name for path in args.panel)
    if names != ORDER:
        raise SystemExit(f"Figure S1 panel order drift: {names} != {ORDER}")
    if args.output.exists() or args.work_dir.exists():
        raise SystemExit("refusing pre-existing Figure S1 composer output/work directory")
    for panel in args.panel:
        if not panel.is_file() or len(PdfReader(panel).pages) != 1:
            raise SystemExit(f"missing or non-single-page Figure S1 panel: {panel}")
    args.work_dir.mkdir(parents=True)
    for panel in args.panel:
        page = PdfReader(panel).pages[0]
        page.scale_by(args.panel_scale)
        writer = PdfWriter(); writer.add_page(page)
        with (args.work_dir / panel.name).open("wb") as handle:
            writer.write(handle)
    (args.work_dir / "assemble.tex").write_text(TEX)
    version = subprocess.run(
        [args.xelatex, "--version"], text=True, stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT, check=True,
    ).stdout.splitlines()[0]
    subprocess.run(
        [args.xelatex, "-interaction=nonstopmode", "-halt-on-error", "assemble.tex"],
        cwd=args.work_dir, check=True,
    )
    built = args.work_dir / "assemble.pdf"
    if len(PdfReader(built).pages) != 3:
        raise SystemExit("established Figure S1 composer did not produce three pages")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    os.replace(built, args.output)
    print(f"S1_COMPOSER_PASS xelatex={version!r} sha256={digest(args.output)}")


if __name__ == "__main__":
    main()
