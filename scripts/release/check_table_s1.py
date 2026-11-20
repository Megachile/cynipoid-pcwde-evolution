#!/usr/bin/env python3
"""Fail the build unless rendered Table S1 carries the built source inventory with RNA-seq shown in resource class.

Checks, independently of the producer:
  - columns, row order and every cell other than resource_class are unchanged;
  - the "Genome + RNA-seq" rows are exactly the deposit Genome rows whose species
    appears in the built expression table;
  - the "Transcriptome (RNA-seq)" rows are exactly the deposit
    "Transcriptome-derived proteome" rows;
  - every other resource_class is unchanged;
  - the caption is the deposit caption plus exactly one sentence;
  - the typeset Table S1 covered all rows.

Environment: project Python 3 (standard library only).
"""
import argparse
import csv
import json
from pathlib import Path


def read_tsv(path):
    with Path(path).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader), list(reader.fieldnames or [])


ap = argparse.ArgumentParser()
ap.add_argument("--data-root", type=Path, required=True)
ap.add_argument("--build", type=Path, required=True)
a = ap.parse_args()
accepted = a.data_root / "repository_inputs/tables/accepted"
base, base_fields = read_tsv(a.build / "table_s1_inventory/Table_S1.tsv")
built, built_fields = read_tsv(a.build / "tables/Table_S1.tsv")
with (Path(__file__).resolve().parents[2] / "config/table_numbers.tsv").open(newline="") as handle:
    labels = {r["table_id"]: r["label"] for r in csv.DictReader(handle, delimiter="\t")}
expression_label = labels["expression"]
expression_species = {r["species"] for r in read_tsv(a.build / f"tables/Table_{expression_label}b.tsv")[0]}
failures = []
if built_fields != base_fields:
    failures.append(f"Table S1 columns {built_fields} != deposit columns {base_fields}")
if len(built) != len(base):
    failures.append(f"Table S1 rows {len(built)} != deposit rows {len(base)}")
genome_rnaseq, transcriptome_rnaseq = set(), set()
for i, (old, new) in enumerate(zip(base, built), start=1):
    if any(old[f] != new[f] for f in base_fields if f != "resource_class"):
        failures.append(f"Table S1 row {i} ({old['species_display']}) changed outside resource_class")
    if new["resource_class"] == "Genome + RNA-seq":
        genome_rnaseq.add(i)
    elif new["resource_class"] == "Transcriptome (RNA-seq)":
        transcriptome_rnaseq.add(i)
    elif new["resource_class"] != old["resource_class"]:
        failures.append(f"Table S1 row {i} ({old['species_display']}): resource class changed to {new['resource_class']!r}")
want_genome = {i for i, r in enumerate(base, start=1) if r["resource_class"] == "Genome" and r["species_display"] in expression_species}
want_transcriptome = {i for i, r in enumerate(base, start=1) if r["resource_class"] == "Transcriptome-derived proteome"}
if genome_rnaseq != want_genome:
    failures.append(f"Genome + RNA-seq rows {sorted(genome_rnaseq)} != genome rows with a Table {expression_label} library {sorted(want_genome)}")
if transcriptome_rnaseq != want_transcriptome:
    failures.append(f"Transcriptome (RNA-seq) rows {sorted(transcriptome_rnaseq)} != transcriptome-derived rows {sorted(want_transcriptome)}")
missing = sorted(expression_species - {r["species_display"] for r in base})
if missing:
    failures.append(f"Table {expression_label} species without a Table S1 row: {missing}")
base_caption = (accepted / "Table_S1_caption.txt").read_text(encoding="utf-8").strip()
built_caption = (a.build / "table_s1/Table_S1_caption.txt").read_text(encoding="utf-8").strip()
if not built_caption.startswith(base_caption + " ") or built_caption[len(base_caption) + 1:].count(". ") != 0:
    failures.append("Table S1 caption is not the deposit caption plus one sentence")
render = json.loads((a.build / "tables/Table_S1.render.json").read_text())
if render["rows"] != len(base) or sum(render["rows_per_page"]) != len(base):
    failures.append("typeset Table S1 does not cover every row")
for f in failures:
    print(f"FAIL\t{f}")
if failures:
    raise SystemExit(1)
print(f"TABLE_S1_PASS rows={len(built)} genome_rnaseq={len(genome_rnaseq)} transcriptome_rnaseq={len(transcriptome_rnaseq)} "
      f"expression_species={len(expression_species)} pages={render['pages']}")
