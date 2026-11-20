#!/usr/bin/env python3
"""Build the gene-copy table from the frozen member layer.

Inputs (deposit, repository_inputs/members/): member_table.tsv (one row per gene
copy) and member_order.tsv (each three-family copy's positions in the ordered
gene list). Three-family copies are grouped into loci by scripts/members/locus_rule.py;
curated B-lectin units are one-copy rows as frozen. Placement fields come from the
member rows; every copy of a locus must agree. Output columns and row order match
the copy table the figures and tables read (loci sorted by locus id, then curated units).

Environment: project Python 3 (standard library only).
"""
import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from locus_rule import group_loci  # noqa: E402

COLUMNS = ["species", "genome", "family", "locus_protein_id", "protein_ids", "annotation_source", "scaffold",
           "placement_route", "bk_chromosome", "reason", "provenance"]
LOCUS_FIELDS = ("species", "annotation_source", "placement_route", "bk_chromosome", "placement_reason")

parser = argparse.ArgumentParser()
parser.add_argument("--data-root", type=Path, required=True)
parser.add_argument("--out", type=Path, required=True)
args = parser.parse_args()
members_dir = args.data_root / "repository_inputs/members"
with (members_dir / "member_table.tsv").open(newline="") as handle:
    members = list(csv.DictReader(handle, delimiter="\t"))
with (members_dir / "member_order.tsv").open(newline="") as handle:
    order = list(csv.DictReader(handle, delimiter="\t"))

rule = {(m["genome"], m["protein_id"]): m for m in members if m["unit_kind"] == "rule"}
curated = [m for m in members if m["unit_kind"] == "curated_unit"]
if len(rule) + len(curated) != len(members):
    raise SystemExit("member_table.tsv has rows that are neither rule copies nor curated units, or duplicates")
loci = group_loci([(o["genome"], o["scaffold"], int(o["rank"]), int(o["position"]), o["protein_id"]) for o in order],
                  {key: m["family"] for key, m in rule.items()})
rows = []
for locus in loci:
    copies = [rule[(locus["genome"], p)] for p in locus["protein_ids"]]
    values = {f: {c[f] for c in copies} for f in LOCUS_FIELDS}
    split = [f for f, v in values.items() if len(v) != 1]
    if split:
        raise SystemExit(f"copies of locus {locus['locus_id']} disagree on {split}")
    one = {f: v.pop() for f, v in values.items()}
    rows.append({"species": one["species"], "genome": locus["genome"], "family": locus["family"],
                 "locus_protein_id": locus["locus_id"], "protein_ids": ",".join(locus["protein_ids"]),
                 "annotation_source": one["annotation_source"], "scaffold": locus["scaffold"],
                 "placement_route": one["placement_route"], "bk_chromosome": one["bk_chromosome"],
                 "reason": one["placement_reason"], "provenance": "member_table.tsv"})
for m in curated:
    rows.append({"species": m["species"], "genome": m["genome"], "family": m["family"], "locus_protein_id": m["unit_id"],
                 "protein_ids": m["protein_id"], "annotation_source": m["annotation_source"], "scaffold": m["curated_scaffold"],
                 "placement_route": m["placement_route"], "bk_chromosome": m["bk_chromosome"], "reason": m["placement_reason"], "provenance": "member_table.tsv"})
# The validated raw PL1_8 locus has no Helixer protein model or ordered-gene
# position. Keep it as one copy unit, outside the annotation-derived locus rule.
with (args.data_root / "repository_inputs/shared/tblastn_presence.tsv").open(newline="") as handle:
    rescues = [r for r in csv.DictReader(handle, delimiter="\t")
               if r["species"] == "Eubothrus centaureae" and r["family"] == "PL1_8"
               and r["final_state"] == "admitted_validated_rescue"]
if len(rescues) != 1 or rescues[0]["reason"] != "fig5_validated_raw_locus_rescue.tsv":
    raise SystemExit("expected one pinned Eubothrus PL1_8 rescue")
with (args.data_root / "repository_inputs/shared/fig5_validated_raw_locus_rescue.tsv").open(newline="") as handle:
    raw = [r for r in csv.DictReader(handle, delimiter="\t")
           if r["species"] == "Isocolus centaureae" and r["pcwde_family"] == "PL1_8"
           and r["figure5_action"] == "add_validated_raw_locus_rescue"]
if len(raw) != 1 or raw[0]["bk_equivalent_chr"] != "unplaced" or raw[0]["best_contig"] != "NODE_38372_length_2218_cov_4.820177":
    raise SystemExit("validated rescue identity or placement changed")
rescue_id = f'{raw[0]["best_contig"]}:{raw[0]["sstart"]}-{raw[0]["send"]}'
if rescue_id not in rescues[0]["best_passing_hsp"]:
    raise SystemExit("rescue source and presence evidence disagree")
rows.append({"species": "Isocolus centaureae", "genome": "Isocolus_centaureae", "family": "PL1_8",
             "locus_protein_id": rescue_id, "protein_ids": rescue_id,
             "annotation_source": "validated raw tBLASTn locus", "scaffold": raw[0]["best_contig"],
             "placement_route": "RagTag", "bk_chromosome": "unplaced",
             "reason": "unassigned_ragtag_scaffold",
             "provenance": "fig5_validated_raw_locus_rescue.tsv;tblastn_presence.tsv"})
args.out.parent.mkdir(parents=True, exist_ok=True)
with args.out.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
copies = sum(len(r["protein_ids"].split(",")) for r in rows)
print(f"COPY_TABLE_PASS loci={len(loci)} curated_units={len(curated)} rows={len(rows)} copies={copies} out={args.out}")
