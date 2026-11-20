#!/usr/bin/env python3
"""Fail the build if Figure 1, Tables S2/S3 or the locus trees disagree with the shared copy table.

Environment: project Python 3 (standard library only).
"""
import argparse
import csv
import json
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from copy_table import FAMILIES, PAPER_NAMES, copy_cells, read_copies, read_tblastn_presence  # noqa: E402

parser = argparse.ArgumentParser()
parser.add_argument("--data-root", type=Path, required=True)
parser.add_argument("--build", type=Path, required=True)
parser.add_argument("--copy-table", type=Path, required=True)
parser.add_argument("--locus-tips", type=Path, required=True)
args = parser.parse_args()
cells = copy_cells(read_copies(args.copy_table))
failures = []

# Full copy inventory retained in the deposit, outside the reader table.
s2 = args.build / "table_s2/candidate/supplementary_tables/Table_S2.tsv"
columns = {"GH5_2": "gh5_2_copies", "PL1_8": "pl1_8_copies", "PL4_4": "pl4_4_copies", "B-lectin": "b_lectin_copies"}
with (args.build / "table_s2/candidate/data/accepted_copy_counts.tsv").open(newline="") as handle:
    part_b = list(csv.DictReader(handle, delimiter="\t"))
s2_species = {r["species_display"] for r in part_b}
for species, _ in cells:
    if species not in s2_species:
        failures.append(f"Deposit copy inventory lacks copy-table species {species}")
for row in part_b:
    for family, column in columns.items():
        expected = cells.get((row["species_display"], family), Counter())["total"]
        if int(row[column]) != expected:
            failures.append(f"Deposit copy inventory {row['species_display']} {family}: {row[column]} != {expected}")

# Table S3: total, named-chromosome and unplaced copies per genome and family.
with (args.build / "tables_s3_s4/candidate/supplementary_tables/Table_S3.tsv").open(newline="") as handle:
    s3 = list(csv.DictReader(handle, delimiter="\t"))
for row in s3:
    counts = cells.get((row["species_display"], row["family_id"]), Counter())
    expected = (counts["total"], counts["total"] - counts["unplaced"], counts["unplaced"])
    found = (int(row["accepted_target_count"]), int(row["named_chromosome_count"]), int(row["unplaced_count"]))
    if found != expected:
        failures.append(f"Table S3 {row['species_display']} {row['family_id']}: {found} != {expected}")
if sum(int(r["accepted_target_count"]) for r in s3) != sum(c["total"] for c in cells.values()):
    failures.append("Table S3 total differs from the copy table")

# Locus trees (Figure 4 / S3): every genome tip is a copy-table copy on the same chromosome.
copy_chromosome = {c["protein_id"]: c["chromosome"] for c in read_copies(args.copy_table)}
with args.locus_tips.open(newline="") as handle:
    genome_tips = [r for r in csv.DictReader(handle, delimiter="\t") if r["source_kind"] == "genome_copy"]
for tip in genome_tips:
    chromosome = copy_chromosome.get(tip["source_protein_id"])
    if chromosome is None:
        failures.append(f"Locus tree tip {tip['display_tree_tip']} is not a copy-table copy")
    elif chromosome != tip["corrected_chromosome"]:
        failures.append(f"Locus tree tip {tip['display_tree_tip']}: {tip['corrected_chromosome']} != copy table {chromosome}")

# Locus-tree figures: the chromosome text drawn for every tip equals the built tip table
# (genome copy: chromosome number, or U if unplaced; transcript: T).
def drawn_chromosome(tip):
    if tip["source_kind"] == "transcript_assembly":
        return "T"
    if tip["corrected_chromosome"] == "unplaced":
        return "U"
    number = tip["corrected_chromosome"].removeprefix("chr")
    return number if number.isdigit() and 1 <= int(number) <= 10 else ""


with args.locus_tips.open(newline="") as handle:
    tip_rows = {(r["family"], r["display_tree_tip"]): r for r in csv.DictReader(handle, delimiter="\t")}
drawn_files = sorted((args.build / "figures/locus_trees/evidence").glob("chromosome_labels_*.tsv")) + \
    sorted((args.build / "figures/complete_locus_trees/evidence").glob("chromosome_labels_*.tsv"))
drawn_views = set()
drawn_labels = 0
for path in drawn_files:
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            drawn_views.add((path.parent.parent.name, row["family"], row["view"]))
            tip = tip_rows.get((row["family"], row["display_tree_tip"]))
            if tip is None:
                failures.append(f"{path.name}: drawn tip {row['display_tree_tip']} is not in the built tip table")
            elif row["label"] != drawn_chromosome(tip):
                failures.append(f"{path.name}: {row['display_tree_tip']} drawn {row['label']!r}, tip table gives {drawn_chromosome(tip)!r}")
            drawn_labels += 1
for figure in ("locus_trees", "complete_locus_trees"):
    for family in ("GH5_2", "PL1_8", "PL4_4"):
        if not any(v[0] == figure and v[1] == family for v in drawn_views):
            failures.append(f"{figure}: no drawn chromosome labels recorded for {family}")

# Figure 5 (B-lectin tree): every drawn tip colour equals the copy table
# (genome copy: its BK chromosome or "unplaced"; transcript tip: "transcript"; bacterium: "outgroup").
b_lectin_chromosome = {c["protein_id"]: c["chromosome"] for c in read_copies(args.copy_table) if c["family"] == "B-lectin"}
with args.copy_table.open(newline="") as handle:
    for locus in csv.DictReader(handle, delimiter="\t"):
        if locus["family"] == "B-lectin":
            b_lectin_chromosome.setdefault(locus["locus_protein_id"], locus["bk_chromosome"])
drawn_b_lectin = args.build / "figures/b_lectin_tree/evidence/drawn_placements.tsv"
b_lectin_tips = 0
if not drawn_b_lectin.exists():
    failures.append("Figure 5: no drawn placements recorded")
else:
    with drawn_b_lectin.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            b_lectin_tips += 1
            if row["tip_role"] == "outgroup":
                expected = "outgroup"
            elif row["biological_unit_id"].startswith("transcript|"):
                expected = "transcript"
            else:
                expected = b_lectin_chromosome.get(row["one_best_tip_id"]) or b_lectin_chromosome.get(row["biological_unit_id"])
            if row["drawn_placement"] != expected:
                failures.append(f"Figure 5: {row['biological_unit_id']} drawn {row['drawn_placement']!r}, copy table gives {expected!r}")

# Figure 1 (full roster): placed copies per chromosome for every cynipoid row.
tip_to_paper = {k.replace(" ", "_"): v for k, v in PAPER_NAMES.items()}
figure = Counter()
shown = set()
with (args.build / "figures/tree_counts/full_combined_cells.tsv").open(newline="") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        if row["position"] == "genome":
            continue
        species = tip_to_paper.get(row["species"], row["species"]).replace("_", " ")
        family = "B-lectin" if row["family"] == "CBM13" else row["family"]
        figure[species, family, row["position"]] += int(row["value"])
        shown.add((species, family))
for (species, family), counts in cells.items():
    if (species, family) not in shown:
        failures.append(f"Figure 1 lacks {species} {family}")
        continue
    placed = counts["total"] - counts["unplaced"]
    drawn = sum(v for (sp, fam, _), v in figure.items() if (sp, fam) == (species, family))
    if drawn != placed:
        failures.append(f"Figure 1 {species} {family}: {drawn} placed copies drawn != {placed}")
    for chromosome in (f"chr{i}" for i in range(1, 11)):
        if figure[species, family, chromosome] != counts[chromosome]:
            failures.append(f"Figure 1 {species} {family} {chromosome}: {figure[species, family, chromosome]} != {counts[chromosome]}")

# Presence (rule b): an accepted copy, or a validated tBLASTn locus where no copy exists.
tblastn = read_tblastn_presence(args.data_root / "repository_inputs/shared/tblastn_presence.tsv")
tblastn_only = {pair for pair in tblastn if not cells.get(pair, Counter())["total"]}
present = {pair for pair, counts in cells.items() if counts["total"]} | tblastn
with s2.open(newline="") as handle:
    part_a = {r["family_id"]: r["genomes_with_occurrence"] for r in csv.DictReader(handle, delimiter="\t")}
presence_counts = {}
for family in FAMILIES:
    presence_counts[family] = sum(f == family for _, f in present)
    shown_a = part_a.get(family.replace("-", "_"))
    if shown_a != f"{presence_counts[family]}/49":
        failures.append(f"Outgroup table {family}: {shown_a} != {presence_counts[family]}/49")
# Figure 1: rows drawn as question marks with no copy are exactly the tBLASTn-only pairs.
pale = set()
with (args.build / "figures/tree_counts/full_combined_display_states.tsv").open(newline="") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        species = tip_to_paper.get(row["species"], row["species"]).replace("_", " ")
        family = "B-lectin" if row["family"] == "CBM13" else row["family"]
        if row["state"] == "unplaced" and not cells.get((species, family), Counter())["total"]:
            pale.add((species, family))
if pale != tblastn_only:
    failures.append(f"Figure 1 tBLASTn-only rows {sorted(pale)} != {sorted(tblastn_only)}")
# numbers.json carries the same counts the products show.
numbers = json.loads((args.build / "numbers.json").read_text())
if numbers["genomes_present_rule_b"] != presence_counts:
    failures.append(f"numbers.json presence {numbers['genomes_present_rule_b']} != products {presence_counts}")
expected_origin = {r["family_id"]: int(r["focal_count"]) + int(r["outgroup_count"]) for r in csv.DictReader(s2.open(newline=""), delimiter="\t") if r["family_id"] in ("GH5_2", "PL1_8", "PL4_4")}
if numbers["origin_domain_totals"] != expected_origin:
    failures.append(f"numbers.json origin totals {numbers['origin_domain_totals']} != Table S2 {expected_origin}")
with s2.open(newline="") as handle:
    expected_focal = {r["family_id"]: int(r["focal_count"]) for r in csv.DictReader(handle, delimiter="\t") if r["family_id"] in ("PL1_8", "PL4_4")}
if numbers["origin_cynipoid_domains"] != expected_focal:
    failures.append(f"numbers.json origin cynipoid counts {numbers['origin_cynipoid_domains']} != Table S2 {expected_focal}")
for family in FAMILIES:
    copies = sum(c["total"] for (_, f), c in cells.items() if f == family)
    s2_copies = sum(int(r[columns[family]]) for r in part_b)
    if not numbers["gene_copies"][family] == copies == s2_copies:
        failures.append(f"numbers.json {family} copies {numbers['gene_copies'][family]} != copy table {copies} / deposit inventory {s2_copies}")
for key, column in (("gene_copies_on_bk_chromosome", "named_chromosome_count"), ("gene_copies_unplaced", "unplaced_count")):
    if numbers[key]["all"] != sum(int(r[column]) for r in s3):
        failures.append(f"numbers.json {key} {numbers[key]['all']} != Table S3 {column}")

if failures:
    for failure in failures:
        print(f"FAIL\t{failure}")
    raise SystemExit(1)
print("PRESENCE_PASS genomes " + "/".join(str(presence_counts[f]) for f in FAMILIES) + f" tblastn_only={len(tblastn_only)}")
print(f"LOCUS_CHROMOSOME_LABELS_PASS views={len(drawn_views)} labels={drawn_labels}")
print(f"FIGURE5_PLACEMENTS_PASS tips={b_lectin_tips}")
print(f"COPY_COUNTS_PASS table_s2_rows={len(part_b)} table_s3_rows={len(s3)} locus_genome_tips={len(genome_tips)} figure1_cells={len(shown)} copy_cells={len(cells)}")
