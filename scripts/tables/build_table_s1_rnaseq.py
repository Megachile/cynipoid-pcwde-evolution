#!/usr/bin/env python3
"""Table S1 with RNA-seq availability shown in its resource class, computed from the built Table S6.

Every row and cell of the built source inventory is carried over unchanged and in order,
except `resource_class`:
  "Genome + RNA-seq"         a Genome row whose species has a library in the built Table S6
  "Transcriptome (RNA-seq)"  every Transcriptome-derived proteome row
  unchanged                  every other row
Species join by exact `species_display` string; no fuzzy matching, so sorted inquiline
bins stay distinct from their hosts. A Table S6 species with no Table S1 row fails the
build. One caption sentence defines the two labels. `join_report.tsv` lists every
Table S6 species and the Table S1 rows it joins.

Environment: project Python 3 (standard library only).
"""
import argparse
import csv
from collections import defaultdict
from pathlib import Path

GENOME, TRANSCRIPTOME = "Genome", "Transcriptome-derived proteome"
GENOME_RNASEQ, TRANSCRIPTOME_RNASEQ = "Genome + RNA-seq", "Transcriptome (RNA-seq)"
with (Path(__file__).resolve().parents[2] / 'config/table_numbers.tsv').open(newline='') as handle:
    LABELS = {r['table_id']:r['label'] for r in csv.DictReader(handle,delimiter='\t')}
CAPTION_SENTENCE = (f"In Resource class, \u201c{GENOME_RNASEQ}\u201d marks genome species with RNA-seq libraries in Table {LABELS['expression']}, "
                    f"and \u201c{TRANSCRIPTOME_RNASEQ}\u201d marks proteomes derived from RNA-seq transcriptome assemblies.")


def read_tsv(path):
    with Path(path).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader), list(reader.fieldnames or [])


def resource_class(row, s6_species):
    if row["resource_class"] == TRANSCRIPTOME:
        return TRANSCRIPTOME_RNASEQ
    if row["resource_class"] == GENOME and row["species_display"] in s6_species:
        return GENOME_RNASEQ
    return row["resource_class"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--table-s1-dir", type=Path, required=True, help="built source inventory directory")
    ap.add_argument("--caption", type=Path, required=True, help="pinned accepted caption")
    ap.add_argument("--table-s6", type=Path, required=True, help="built Table_S6.tsv")
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    if a.out.exists():
        raise SystemExit(f"refusing existing output: {a.out}")
    rows, fields = read_tsv(a.table_s1_dir / "Table_S1.tsv")
    if any(r["resource_class"] in (GENOME_RNASEQ, TRANSCRIPTOME_RNASEQ) for r in rows):
        raise SystemExit("deposit Table S1 already carries RNA-seq resource classes")
    s6_rows, _ = read_tsv(a.table_s6)
    s6_species = {r["species"] for r in s6_rows}
    s6_counts = defaultdict(int)
    for r in s6_rows:
        s6_counts[r["species"]] += 1
    joined = defaultdict(list)
    for row in rows:
        if row["species_display"] in s6_species:
            joined[row["species_display"]].append(row["resource_class"])
    a.out.mkdir(parents=True)
    with (a.out / "Table_S1.tsv").open("w", newline="", encoding="utf-8") as handle:
        w = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for row in rows:
            w.writerow({**row, "resource_class": resource_class(row, s6_species)})
    caption = a.caption.read_text(encoding="utf-8").strip()
    (a.out / "Table_S1_caption.txt").write_text(f"{caption} {CAPTION_SENTENCE}\n", encoding="utf-8")
    failures = sorted(s for s in s6_species if s not in joined)
    with (a.out / "join_report.tsv").open("w", newline="", encoding="utf-8") as handle:
        w = csv.writer(handle, delimiter="\t", lineterminator="\n")
        w.writerow(["table_s6_species", "table_s6_rows", "table_s1_rows_joined", "status"])
        for s in sorted(s6_species):
            w.writerow([s, s6_counts[s], ";".join(joined.get(s, [])), "joined" if s in joined else "NO_TABLE_S1_ROW"])
    if failures:
        raise SystemExit(f"Table S6 species without a Table S1 row: {failures}")
    classes = [resource_class(r, s6_species) for r in rows]
    print(f"TABLE_S1_RNASEQ_PASS rows={len(rows)} genome_rnaseq={classes.count(GENOME_RNASEQ)} "
          f"transcriptome_rnaseq={classes.count(TRANSCRIPTOME_RNASEQ)} table_s6_species={len(s6_species)} join_failures=0")


if __name__ == "__main__":
    main()
