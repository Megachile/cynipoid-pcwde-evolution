#!/usr/bin/env python3
"""Render every supplement page to PNG with Ghostscript and assemble a contact sheet.

Inputs: the built supplement PDF. Outputs (in --out): page_NNN.png at --dpi,
contact_sheet.png, PAGES.tsv (page, file, sha256, pixel size) and RENDER.json
(Ghostscript version, resolution, source PDF sha256, page count). Fails unless
one image exists per PDF page.

Environment: project Python (Pillow, PyMuPDF from environment/python-lock.txt)
and a Ghostscript executable given by --gs (make variable GS).
"""
import argparse
import csv
import hashlib
import json
import math
import subprocess
from pathlib import Path

import fitz
from PIL import Image, ImageDraw

parser = argparse.ArgumentParser()
parser.add_argument("--gs", required=True)
parser.add_argument("--pdf", type=Path, required=True)
parser.add_argument("--out", type=Path, required=True)
parser.add_argument("--dpi", type=int, default=150)
args = parser.parse_args()
if args.out.exists():
    raise SystemExit(f"refusing existing output: {args.out}")
args.out.mkdir(parents=True)


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


version = subprocess.run([args.gs, "--version"], check=True, capture_output=True, text=True).stdout.strip()
subprocess.run([args.gs, "-q", "-dNOPAUSE", "-dBATCH", "-dSAFER", "-sDEVICE=png16m", f"-r{args.dpi}",
                "-dTextAlphaBits=4", "-dGraphicsAlphaBits=4", f"-sOutputFile={args.out / 'page_%03d.png'}",
                str(args.pdf)], check=True)
pages = sorted(args.out.glob("page_*.png"))
expected = len(fitz.open(args.pdf))
if len(pages) != expected or [p.name for p in pages] != [f"page_{i:03d}.png" for i in range(1, expected + 1)]:
    raise SystemExit(f"rendered {len(pages)} page images for a {expected}-page PDF")

rows = []
cols, thumb = 8, 220
thumbs = []
for number, path in enumerate(pages, 1):
    image = Image.open(path).convert("RGB")
    rows.append(dict(page=number, file=path.name, sha256=sha256(path), width_px=image.width, height_px=image.height))
    image.thumbnail((thumb, thumb * 2))
    thumbs.append(image)
cell = max(im.height for im in thumbs) + 18
sheet = Image.new("RGB", (cols * (thumb + 10) + 10, math.ceil(len(thumbs) / cols) * (cell + 10) + 10), "white")
draw = ImageDraw.Draw(sheet)
for i, image in enumerate(thumbs):
    x, y = 10 + (i % cols) * (thumb + 10), 10 + (i // cols) * (cell + 10)
    sheet.paste(image, (x, y))
    draw.rectangle([x, y, x + image.width - 1, y + image.height - 1], outline="grey")
    draw.text((x, y + image.height + 3), f"p{i + 1}", fill="black")
sheet.save(args.out / "contact_sheet.png")

with (args.out / "PAGES.tsv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, list(rows[0]), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
(args.out / "RENDER.json").write_text(json.dumps({
    "renderer": "ghostscript", "ghostscript_version": version, "device": "png16m", "dpi": args.dpi,
    "source_pdf": args.pdf.name, "source_pdf_sha256": sha256(args.pdf), "pages": expected,
    "contact_sheet_sha256": sha256(args.out / "contact_sheet.png"),
}, indent=2) + "\n")
print(f"SUPPLEMENT_PAGES_PASS pages={expected} gs={version} dpi={args.dpi}")
