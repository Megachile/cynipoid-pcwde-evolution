#!/usr/bin/env python3
"""Rebuild the accepted public Figure S4 B-lectin display.

The admitted historical producer called this object Figure S3. This portable
interface uses the final public number throughout while retaining the accepted
68-domain membership (the accepted 69 domains less one CB-derived Bassettia
pallida transcript domain, pruned after inference), rooted tree, sequences,
alignment, complete architectures, and display-only taxonomy substitutions.
The drawing and its QC are regenerated from the tree and architecture tables.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import html
import re
import shutil
import os
import subprocess
from pathlib import Path

import cairosvg
from Bio import Phylo, SeqIO
from PIL import Image


EXPECTED_RECORDS = 68
TAXONOMY_REPLACEMENTS = {
    "Eschatocerus acacia": "Eschatocerus acaciae",
    "Isocolus centaureae": "Eubothrus centaureae",
    "Ceroptres masudai": "Ceroptres kovalevi",
}
FORBIDDEN_BYTES = (
    b"/" + b"carc" + b"/scratch",
    b"/" + b"easley" + b"/scratch",
    b"dryocosmus_" + b"cannoni",
    b"annotation:" + b"12",
    b"source-" + b"12",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fasta_ids(path: Path) -> list[str]:
    values = [record.id for record in SeqIO.parse(path, "fasta")]
    if not values or len(values) != len(set(values)):
        raise AssertionError(f"empty or duplicate FASTA identifiers: {path}")
    return values


def tree_tips(path: Path) -> list[str]:
    values = [tip.name for tip in Phylo.read(path, "newick").get_terminals()]
    if not values or None in values or len(values) != len(set(values)):
        raise AssertionError(f"empty or duplicate tree tips: {path}")
    return values


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def normalize_text(value: str) -> str:
    for old, new in TAXONOMY_REPLACEMENTS.items():
        value = value.replace(old, new)
    return value


def copy_science(input_dir: Path, out: Path) -> None:
    tree_source = input_dir / "B_lectin_domains.rooted.treefile"
    sequence_source = input_dir / "B_lectin_domains.faa"
    alignment_source = input_dir / "B_lectin_domains.aln.faa"
    tips = set(tree_tips(tree_source))
    sequences = set(fasta_ids(sequence_source))
    alignment = set(fasta_ids(alignment_source))
    if len(tips) != EXPECTED_RECORDS or tips != sequences or tips != alignment:
        raise AssertionError("Figure S4 tree/sequence/alignment membership mismatch")
    destinations = (
        (tree_source, out / "trees/Figure_S4_B_lectin_domains.rooted.treefile"),
        (sequence_source, out / "sequences/Figure_S4_B_lectin_domains.faa"),
        (alignment_source, out / "alignments/Figure_S4_B_lectin_domains.aln.faa"),
    )
    for source, target in destinations:
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        if sha256(source) != sha256(target):
            raise AssertionError(f"science copy checksum mismatch: {target}")


def build_sidecars(input_dir: Path, out: Path) -> None:
    for source_name, target_name in (
        ("domain_manifest.tsv", "Figure_S4_domain_manifest.tsv"),
        ("architecture.tsv", "Figure_S4_architecture.tsv"),
    ):
        rows = read_tsv(input_dir / source_name)
        if len(rows) != EXPECTED_RECORDS:
            raise AssertionError(f"expected {EXPECTED_RECORDS} rows: {source_name}")
        for row in rows:
            for field in ("species", "display_label"):
                if field in row:
                    row[field] = normalize_text(row[field])
        write_tsv(out / "sidecars" / target_name, rows, list(rows[0]))


def build_visuals(input_dir: Path, out: Path) -> None:
    target = out / "figures/Figure_S4.svg"
    target.parent.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, S3_INTERFACE=str(out / "sidecars/Figure_S4_architecture.tsv"),
               S3_TREE_FILE=str(input_dir / "B_lectin_domains.rooted.treefile"),
               S3_OUTPUT_DIR=str(target.parent), S3_EXPECTED_TIPS=str(EXPECTED_RECORDS))
    subprocess.run([os.environ.get("RENDER_RSCRIPT", "Rscript"),
                    str(Path(__file__).with_name("render_domains.R"))], env=env, check=True)
    (out / "qc").mkdir(exist_ok=True)
    for name in ("render_qc", "terminal_order", "rendered_domains", "plot_data"):
        shutil.move(str(target.parent / f"Figure_S4_{name}.tsv"),
                    str(out / "qc" / f"Figure_S4_{name}.tsv"))
    svg = target.read_text()
    tip_count = 0

    def italicize_tip(match: re.Match[str]) -> str:
        nonlocal tip_count
        attrs, escaped = match.groups()
        if "font-size: 8.96px" not in attrs:
            return match.group(0)
        tip_count += 1
        label = html.unescape(escaped)
        ex = re.fullmatch(r"([A-Z][A-Za-z-]+) sp\. ex ([A-Z][A-Za-z-]+ [a-z][a-z-]+)", label)
        if ex:
            pieces = [(ex.group(1), "italic"), (" sp. ex ", "normal"),
                      (ex.group(2), "italic")]
        else:
            binomial = re.fullmatch(r"([A-Z][A-Za-z-]+ [a-z][a-z-]+)(.*)", label)
            if not binomial:
                raise AssertionError(f"unhandled S6 tip label: {label}")
            pieces = [(binomial.group(1), "italic")]
            if binomial.group(2):
                pieces.append((binomial.group(2), "normal"))
        body = "".join(f'<tspan style="font-style: {style};">{html.escape(text)}</tspan>'
                       for text, style in pieces)
        return f"<text{attrs}>{body}</text>"

    svg = re.sub(r"<text\b([^>]*)>([^<]*)</text>", italicize_tip, svg)
    if tip_count != EXPECTED_RECORDS:
        raise AssertionError(f"expected {EXPECTED_RECORDS} S6 labels, got {tip_count}")
    target.write_text(svg)
    cairosvg.svg2pdf(url=str(target), write_to=str(target.with_suffix(".pdf")))
    cairosvg.svg2png(
        url=str(target), write_to=str(target.with_suffix(".png")), output_width=4200
    )
    with Image.open(target.with_suffix(".png")) as image:
        if image.width < 1000 or image.height < 1000 or image.convert("RGB").getbbox() is None:
            raise AssertionError(f"implausible Figure S4 raster geometry: {image.size}")


def validate_sanitization(out: Path) -> None:
    for path in sorted(item for item in out.rglob("*") if item.is_file()):
        data = path.read_bytes().lower()
        for forbidden in FORBIDDEN_BYTES:
            if forbidden in data:
                raise AssertionError(f"forbidden public token in {path}")


def build_payload(input_dir: Path, out: Path) -> None:
    if out.exists():
        raise AssertionError(f"refusing to overwrite output: {out}")
    out.mkdir(parents=True)
    copy_science(input_dir, out)
    build_sidecars(input_dir, out)
    build_visuals(input_dir, out)
    write_tsv(
        out / "qc/Figure_S4_membership_qc.tsv",
        [{
            "tree_tips": EXPECTED_RECORDS,
            "sequence_records": EXPECTED_RECORDS,
            "alignment_records": EXPECTED_RECORDS,
            "rooting": "accepted two-domain bacterial MRCA rooting retained",
            "status": "PASS",
        }],
        ["tree_tips", "sequence_records", "alignment_records", "rooting", "status"],
    )
    validate_sanitization(out)


def compare_rebuild(primary: Path, secondary: Path) -> None:
    first = {str(path.relative_to(primary)): path for path in primary.rglob("*") if path.is_file()}
    second = {str(path.relative_to(secondary)): path for path in secondary.rglob("*") if path.is_file()}
    if set(first) != set(second):
        raise AssertionError("Figure S4 deterministic rebuild file-set mismatch")
    rows = []
    for relative_path in sorted(first):
        first_hash = sha256(first[relative_path])
        second_hash = sha256(second[relative_path])
        identical = first_hash == second_hash
        if not identical and Path(relative_path).suffix.lower() != ".pdf":
            raise AssertionError(f"nondeterministic non-PDF output: {relative_path}")
        rows.append({
            "relative_path": relative_path,
            "primary_sha256": first_hash,
            "rebuild_sha256": second_hash,
            "byte_identical": "yes" if identical else "no",
            "status": "PASS",
        })
    write_tsv(primary / "qc/deterministic_rebuild_qc.tsv", rows, list(rows[0]))


def finalize_manifest(out: Path) -> None:
    files = sorted(
        path for path in out.rglob("*")
        if path.is_file() and path.name != "PUBLIC_MANIFEST.tsv"
    )
    rows = [{
        "relative_path": str(path.relative_to(out)),
        "size_bytes": path.stat().st_size,
        "sha256": sha256(path),
    } for path in files]
    write_tsv(out / "PUBLIC_MANIFEST.tsv", rows, list(rows[0]))
    validate_sanitization(out)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--rebuild", type=Path, required=True)
    args = parser.parse_args()
    input_dir = args.input_dir.resolve()
    if not input_dir.is_dir():
        raise SystemExit(f"input directory does not exist: {input_dir}")
    build_payload(input_dir, args.out)
    build_payload(input_dir, args.rebuild)
    compare_rebuild(args.out, args.rebuild)
    finalize_manifest(args.out)
    print(json.dumps({"status": "PASS", "figure": "S4", "records": EXPECTED_RECORDS}, sort_keys=True))


if __name__ == "__main__":
    main()
