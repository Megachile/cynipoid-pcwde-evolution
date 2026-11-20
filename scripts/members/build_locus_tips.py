#!/usr/bin/env python3
"""Build the locus-tree tip table from the frozen tip crosswalk and the built copy table.

Inputs: repository_inputs/members/locus_tip_crosswalk.tsv (frozen with the accepted
2026-10-23 locus trees: each display tip, its raw tree tip, source kind and source
protein; species only for transcript tips) and the copy table built by
build_copy_table.py. Genome-copy tips take species and BK chromosome from the copy
table through their source protein; transcript tips keep their recorded species and
have no chromosome. Output columns and row order match the table the locus-tree
renderers and the count check read.

Environment: project Python 3 (standard library only).
"""
import argparse
import csv
from pathlib import Path

COLUMNS = ["family", "display_tree_tip", "raw_tip", "species", "source_kind", "source_protein_id", "corrected_chromosome"]

parser = argparse.ArgumentParser()
parser.add_argument("--data-root", type=Path, required=True)
parser.add_argument("--copy-table", type=Path, required=True)
parser.add_argument("--out", type=Path, required=True)
args = parser.parse_args()
with args.copy_table.open(newline="") as handle:
    copies = {}
    for row in csv.DictReader(handle, delimiter="\t"):
        for protein in row["protein_ids"].split(","):
            copies[protein.strip()] = row
with (args.data_root / "repository_inputs/members/locus_tip_crosswalk.tsv").open(newline="") as handle:
    tips = list(csv.DictReader(handle, delimiter="\t"))
rows, kinds = [], {"genome_copy": 0, "transcript_assembly": 0}
for tip in tips:
    kind = tip["source_kind"]
    if kind == "genome_copy":
        copy = copies.get(tip["source_protein_id"])
        if copy is None:
            raise SystemExit(f"locus tip {tip['display_tree_tip']}: source protein is not a copy-table copy")
        if copy["family"] != tip["family"]:
            raise SystemExit(f"locus tip {tip['display_tree_tip']}: family {tip['family']} != copy table {copy['family']}")
        species, chromosome = copy["species"], copy["bk_chromosome"]
    elif kind == "transcript_assembly":
        species, chromosome = tip["transcript_species"], ""
    else:
        raise SystemExit(f"locus tip {tip['display_tree_tip']}: unknown source kind {kind}")
    kinds[kind] += 1
    rows.append({"family": tip["family"], "display_tree_tip": tip["display_tree_tip"], "raw_tip": tip["raw_tip"],
                 "species": species, "source_kind": kind, "source_protein_id": tip["source_protein_id"],
                 "corrected_chromosome": chromosome})
args.out.parent.mkdir(parents=True, exist_ok=True)
with args.out.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
print(f"LOCUS_TIPS_PASS tips={len(rows)} genome_copy={kinds['genome_copy']} transcript={kinds['transcript_assembly']} out={args.out}")
