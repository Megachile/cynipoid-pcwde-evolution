#!/usr/bin/env python3
"""Fail closed when current intron and GC3 paper products disagree."""
import csv
import json
import sys
from pathlib import Path


def rows(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def same_file(left, right):
    if not left.is_file() or not right.is_file() or left.read_bytes() != right.read_bytes():
        raise SystemExit(f"INTRON_PRODUCTS_FAIL different/missing: {left} {right}")


def check_projection(technical, reader, columns, technical_fields=None):
    """The reader table must be exactly the declared display projection of the technical table:
    same rows in the same order, declared headings, and each value equal to its declared rounding."""
    tech_rows, reader_rows = rows(technical), rows(reader)
    if not tech_rows or list(tech_rows[0]) != (technical_fields or [c for c, _, _ in columns]):
        raise SystemExit(f"INTRON_PRODUCTS_FAIL technical schema differs: {technical}")
    if not reader_rows or list(reader_rows[0]) != [label for _, label, _ in columns]:
        raise SystemExit(f"INTRON_PRODUCTS_FAIL reader headings differ: {reader}")
    if len(tech_rows) != len(reader_rows):
        raise SystemExit(f"INTRON_PRODUCTS_FAIL reader row count differs: {reader}")
    for i, (tech, shown) in enumerate(zip(tech_rows, reader_rows)):
        for column, label, fmt in columns:
            if shown[label] != fmt(tech[column]):
                raise SystemExit(f"INTRON_PRODUCTS_FAIL {reader.name} row {i + 1} {label!r}: "
                                 f"{shown[label]!r} is not the display of {tech[column]!r}")


def final_table_label(table_id):
    label_path = Path(__file__).resolve().parents[2] / "config/table_numbers.tsv"
    labels = {row["table_id"]: row["label"] for row in rows(label_path)}
    label = labels.get(table_id)
    if not label or not label.startswith("S") or not label[1:].isdigit():
        raise SystemExit(f"INTRON_PRODUCTS_FAIL missing final label for {table_id}")
    return label


def main(build):
    intron = build / "introns"
    numbers = json.loads((build / "numbers.json").read_text())
    gc3 = json.loads((intron / "gc3/pcwde_gc3_stats.json").read_text())
    blectin = json.loads((intron / "gc3/blectin_gc3_stats.json").read_text())
    if numbers["gc3"]["pcwde"] != gc3 or numbers["gc3"]["blectin"] != blectin:
        raise SystemExit("INTRON_PRODUCTS_FAIL GC3 numbers differ from source statistics")
    comparison = rows(intron / "gc3/blectin_gc3_comparison.tsv")
    if len(comparison) != 28 or {r["bacterial_accession"] for r in comparison} != {
            numbers["gc3"]["blectin_bacterial_comparator_accession"]}:
        raise SystemExit("INTRON_PRODUCTS_FAIL bacterial comparator differs")
    table = rows(intron / "gc3/Table_S16.tsv")
    if len(table) != 6:
        raise SystemExit("INTRON_PRODUCTS_FAIL GC3 table has wrong row count")
    for row in table:
        stat = gc3[row["cohort"]] if row["cohort"] != "B-lectin" else blectin[row["stratum"]]
        median = stat.get("median_species_host_percentile", stat.get("species_median_host_percentile"))
        if (int(row["copies"]), int(row["species"]), float(row["median_host_gc3_percentile"]),
            int(row["species_above_host_median"])) != (stat["copies"], stat["species"], median,
                                                    stat["species_above_host_median"]):
            raise SystemExit("INTRON_PRODUCTS_FAIL GC3 table differs from source statistics")
        copy = stat.get("host_distance_pp", stat.get("median_copy_distance_pp"))
        bacterial = stat.get("bacterial_distance_pp", stat.get("median_bacterial_distance_pp"))
        p = stat.get("wilcoxon_bacteria_minus_host_distance_p", stat.get("wilcoxon_distance_p_two_sided"))
        if (float(row["median_copy_host_distance_pp"]), float(row["median_bacterial_host_distance_pp"]),
            float(row["paired_wilcoxon_p"])) != (copy, bacterial, p):
            raise SystemExit("INTRON_PRODUCTS_FAIL GC3 distances differ")
        for key in ("median_copy_signed_offset_pp", "median_bacterial_signed_offset_pp"):
            if float(row[key]) != stat[key]:
                raise SystemExit("INTRON_PRODUCTS_FAIL GC3 signed offsets differ")
    one_decimal = lambda v: f"{float(v):.1f}"
    integer = lambda v: str(int(v))
    check_projection(intron / "gc3/Table_S16.tsv", build / f"tables/Table_{final_table_label('gc3')}.tsv", [
        ("cohort", "Gene family", str), ("stratum", "Copy stratum", str),
        ("copies", "Copies (n)", integer), ("species", "Species (n)", integer),
        ("median_host_gc3_percentile",
         "Median transferred-gene percentile in cynipoid genomic background", one_decimal),
        ("species_above_host_median", "Species with median family-copy percentile above 50 (n)", integer),
        ("median_copy_host_distance_pp",
         "Median absolute transferred-gene difference from genomic background (pp)", one_decimal),
        ("median_bacterial_host_distance_pp",
         "Median absolute bacterial difference from genomic background (pp)", one_decimal),
        ("paired_wilcoxon_p", "Paired Wilcoxon P (two-sided; species units)", lambda v: f"{float(v):.1e}"),
        ("median_copy_signed_offset_pp", "Median signed transferred-gene offset from genomic background (pp)", lambda v: f"{float(v):+.1f}"),
        ("median_bacterial_signed_offset_pp", "Median signed bacterial offset from genomic background (pp)", lambda v: f"{float(v):+.1f}"),
    ])
    pct = lambda v: f"{float(v):.1f}".rstrip("0").rstrip(".") if v else v
    intron_fields = list(rows(intron / "table_s5/Table_S5.tsv")[0])
    intron_display = {"observation_count": ("copies with this intron (n)", str),
                      "observation_proportion_pct": ("copies with this intron (%)", pct),
                      "protein_prevalence_pct": ("protein_prevalence_pct", pct),
                      "taxonomic_families": ("cynipoid families with this intron", str)}
    check_projection(intron / "table_s5/Table_S5.tsv", build / f"tables/Table_{final_table_label('intron')}.tsv",
                     [(f, *intron_display.get(f, (f, str))) for f in intron_fields if f != 'taxonomic_breadth'],
                     intron_fields)
    source_s5 = rows(intron / "table_s5/Table_S5.tsv")
    taxonomy = rows(Path(__file__).resolve().parents[2] / "config/cazy_lineage_groups.tsv")
    species_family = {r["species"].replace("_", " "): r["source_family"] for r in taxonomy}
    family_order = ("Diplolepididae", "Figitidae", "Cynipidae")
    for family in ("GH5_2", "PL1_8", "PL4_4"):
        members = {}
        for position in rows(intron / "summary" / f"{family}_intron_positions.tsv"):
            if position["group_id"]:
                species = position["species"]
                source_family = species_family.get(species)
                if source_family not in family_order:
                    raise SystemExit(f"INTRON_PRODUCTS_FAIL missing cynipoid family for {species!r}")
                members.setdefault(position["group_id"], set()).add(source_family)
        for site in (r for r in source_s5 if r["part"] == "B" and r["family_id"] == family):
            group_id = f"{family}_g{site['reference_position_aa']}"
            expected = ", ".join(f for f in family_order if f in members.get(group_id, set()))
            if not expected or site["taxonomic_families"] != expected:
                raise SystemExit(f"INTRON_PRODUCTS_FAIL incorrect families for {group_id}")
    part_a = {r["family_id"]: {k: int(r[k]) for k in
              ("eligible_proteins", "proteins_with_introns", "total_introns")}
              for r in source_s5 if r["part"] == "A"}
    if numbers["introns"]["table_s5_part_a"] != part_a:
        raise SystemExit("INTRON_PRODUCTS_FAIL Table S5 counts differ")
    site_rows = {(r["family_id"], r["reference_position_aa"]): r for r in source_s5 if r["part"] == "B"}
    for key, fact in numbers["introns"]["table_s5_sites"].items():
        family, site = key.split(":")
        source = site_rows[(family, site)]
        if fact != {"observations": int(source["observation_count"]), "phase": source["phase"]}:
            raise SystemExit("INTRON_PRODUCTS_FAIL Table S5 site differs")
    deposit = build / "data_deposit/pcwde_supporting_data/data/introns"
    for source in (intron / "table_s5/Table_S5.tsv", intron / "gc3/pcwde_gc3_stats.json",
                   intron / "gc3/blectin_gc3_stats.json", intron / "gc3/pcwde_species_pairs.tsv",
                   intron / "gc3/blectin_gc3_comparison.tsv",
                   intron / "summary/recurrent_position_summary.tsv"):
        same_file(source, deposit / source.name)
    for family in ("GH5_2", "PL1_8", "PL4_4"):
        for suffix in ("intron_positions", "groups", "context_rich", "context_summary"):
            same_file(intron / "summary" / f"{family}_{suffix}.tsv",
                      intron / "figure_inputs" / f"{family}_{suffix}.tsv")
    same_file(intron / "host_background/background_intron_details.tsv",
              intron / "figure_inputs/baseline_intron_details.tsv")
    observations = rows(intron / "observations/per_intron_observations.tsv")
    focal = [r for r in observations if r["family"] in ("GH5_2", "PL1_8", "PL4_4") and int(r["intron_number"]) > 0]
    if numbers["introns"]["first_cds_quarter"] != {
        "introns": sum(float(r["relative_cds_position"]) <= .25 for r in focal), "total": len(focal)}:
        raise SystemExit("INTRON_PRODUCTS_FAIL first-quarter numbers differ")
    print("INTRON_PRODUCTS_PASS table_s5 gc3_table reader_projections numbers deposit figure_inputs")


if __name__ == "__main__":
    main(Path(sys.argv[1]))
