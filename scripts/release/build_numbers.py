#!/usr/bin/env python3
"""Write numbers.json: in-text counts computed once from the shared intermediates.

Covers gene copies (total, on a BK chromosome, unplaced) and genome presence
under rule b (accepted copy or validated tBLASTn-only locus) per family.
check_copy_counts.py fails the build if any product disagrees with these values.

Environment: project Python 3 (standard library only).
"""
import argparse
import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from copy_table import FAMILIES, copy_cells, read_copies, read_tblastn_presence  # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parent))
import supplement_numbers  # noqa: E402
from blectin_domain_counts import load_counts as blectin_domain_counts  # noqa: E402

parser = argparse.ArgumentParser()
parser.add_argument("--data-root", type=Path, required=True)
parser.add_argument("--out", type=Path, required=True)
parser.add_argument("--copy-table", type=Path, required=True)
parser.add_argument("--table-s2", type=Path, required=True)
parser.add_argument("--domain-summary", type=Path, required=True)
parser.add_argument("--origin-key-nodes", type=Path, required=True)
parser.add_argument("--new-tables", type=Path, required=True)
parser.add_argument("--rendered-tables", type=Path, required=True)
parser.add_argument("--synteny-table", type=Path, required=True)
parser.add_argument("--synteny-summary", type=Path, required=True)
parser.add_argument("--chromosome-gain-summary", type=Path, required=True)
parser.add_argument("--chromosome-gain-table", type=Path, required=True)
parser.add_argument("--cazy-rank-table", type=Path, required=True)
parser.add_argument("--intron-summary", type=Path, required=True)
parser.add_argument("--intron-observations", type=Path, required=True)
parser.add_argument("--intron-composition", type=Path, required=True)
parser.add_argument("--intron-comparator", type=Path, required=True)
parser.add_argument("--table-s5", type=Path, required=True)
parser.add_argument("--gc3-stats", type=Path, required=True)
parser.add_argument("--blectin-gc3-stats", type=Path, required=True)
parser.add_argument("--blectin-gc3-comparison", type=Path, required=True)
parser.add_argument("--source-panel", type=Path, required=True)
parser.add_argument("--expression-table", type=Path, required=True)
parser.add_argument("--expression-sample-count-matrix", type=Path, required=True)
parser.add_argument("--family-membership", type=Path, required=True)
parser.add_argument("--sequence-sets", type=Path, required=True)
parser.add_argument("--b-lectin-manifest", type=Path, required=True)
parser.add_argument("--b-lectin-tree-counts", type=Path, required=True)
parser.add_argument("--catalytic-domains", type=Path, required=True)
args = parser.parse_args()
with args.expression_sample_count_matrix.open(newline="", encoding="utf-8") as handle:
    sample_count_rows = list(csv.DictReader(handle, delimiter="\t"))
with args.expression_table.open(newline="", encoding="utf-8") as handle:
    expression_rows = list(csv.DictReader(handle, delimiter="\t"))
matrix_cells = {(r["species"], field): int(value)
                for r in sample_count_rows for field, value in r.items()
                if field not in ("species", "source") and value}
if len(sample_count_rows) != 15 or len(matrix_cells) != 37 or sum(matrix_cells.values()) != 158:
    raise SystemExit("Table S2a sample-count dimensions or total changed")
context_columns = {
    ("adult", "accessory_gland"): "adult_accessory_gland",
    ("adult", "venom_gland"): "adult_venom_gland",
    ("adult", "carcass"): "adult_carcass",
    ("adult", "whole_body"): "adult_whole_body",
    ("larva", "gut"): "larval_gut",
    ("larva", "salivary_gland"): "larval_salivary_gland",
    ("larva", "malpighian_tubule"): "larval_malpighian_tubule",
    ("larva", "whole_body"): "larval_whole_body",
    ("larva", "whole_body (in gall)"): "larval_whole_body",
}
expression_cells = {}
for r in expression_rows:
    key = (r["species"], context_columns[(r["life_stage"], r["tissue"])])
    original = (key, r["tissue"])
    count = int(r["sample_count"])
    if original in expression_cells and expression_cells[original] != count:
        raise SystemExit(f"Table S2b family count mismatch: {original}")
    expression_cells[original] = count
merged_expression = {}
for (key, _), count in expression_cells.items():
    merged_expression[key] = merged_expression.get(key, 0) + count
for key, count in merged_expression.items():
    if matrix_cells.get(key) != count:
        raise SystemExit(f"Table S2a/S2b sample-count mismatch: {key}")
shared = args.data_root / "repository_inputs/shared"
cells = copy_cells(read_copies(args.copy_table))
tblastn = read_tblastn_presence(shared / "tblastn_presence.tsv")


def per_family(value):
    out = {family: sum(value(counts) for (_, fam), counts in cells.items() if fam == family) for family in FAMILIES}
    out["all"] = sum(out.values())
    return out


present = {pair for pair, counts in cells.items() if counts["total"]} | tblastn
with args.table_s2.open(newline="", encoding="utf-8") as fh:
    origin_rows = {r["family_id"]: r for r in csv.DictReader(fh, delimiter="\t")}
origin_totals = {f: int(origin_rows[f]["focal_count"]) + int(origin_rows[f]["outgroup_count"]) for f in ("GH5_2", "PL1_8", "PL4_4")}
origin_cynipoid = {f: int(origin_rows[f]["focal_count"]) for f in ("PL1_8", "PL4_4")}
numbers = {
    "blectin_domain_tree": blectin_domain_counts(args.b_lectin_tree_counts),
    "table_s2_sample_counts": {"species": len(sample_count_rows), "nonempty_cells": len(matrix_cells),
                               "profiles": sum(matrix_cells.values())},
    "origin_domain_totals": origin_totals,
    "origin_cynipoid_domains": origin_cynipoid,
    "gene_copies": per_family(lambda c: c["total"]),
    "gene_copies_on_bk_chromosome": per_family(lambda c: c["total"] - c["unplaced"]),
    "gene_copies_unplaced": per_family(lambda c: c["unplaced"]),
    "genomes_present_rule_b": {family: sum(f == family for _, f in present) for family in FAMILIES},
    "genomes_tblastn_only": {family: sum(f == family for _, f in tblastn if not cells.get((_, f), {}).get("total", 0)) for family in FAMILIES},
}
import statistics
with args.intron_observations.open(newline="") as handle:
    intron_rows = [r for r in csv.DictReader(handle, delimiter="\t")
                   if r["family"] in ("GH5_2", "PL1_8", "PL4_4") and int(r["intron_number"]) > 0]
with args.intron_summary.open(newline="") as handle:
    author_intron_summary = list(csv.DictReader(handle, delimiter="\t"))
with args.intron_composition.open(newline="") as handle:
    focal_rows = [r for r in csv.DictReader(handle, delimiter="\t")
                  if r["family"] in ("GH5_2", "PL1_8", "PL4_4")]
with args.table_s5.open(newline="") as handle:
    s5_rows = list(csv.DictReader(handle, delimiter="\t"))
sites = {(r["family_id"], r["reference_position_aa"]):
         {"observations": int(r["observation_count"]), "phase": r["phase"]}
         for r in s5_rows if r["part"] == "B"}
numbers["introns"] = {
    "first_cds_quarter": {"introns": sum(float(r["relative_cds_position"]) <= .25 for r in intron_rows),
                          "total": len(intron_rows)},
    "median_relative_cds_position_pct": 100 * statistics.median(float(r["relative_cds_position"]) for r in intron_rows),
    "median_density_per_kb": statistics.median(float(r["intron_density_per_kb"]) for r in focal_rows),
    "whole_genome_comparator": json.loads(args.intron_comparator.read_text()),
    "table_s5_part_a": {r["family_id"]: {k: int(r[k]) for k in
                          ("eligible_proteins", "proteins_with_introns", "total_introns")}
                        for r in s5_rows if r["part"] == "A"},
    "table_s5_sites": {f"{family}:{site}": sites[(family, site)] for family, site in
                       (("GH5_2", "18"), ("PL1_8", "9"), ("PL4_4", "37"), ("PL4_4", "83"))},
}
numbers["gc3"] = {"pcwde": json.loads(args.gc3_stats.read_text()),
                  "blectin": json.loads(args.blectin_gc3_stats.read_text())}
with args.blectin_gc3_comparison.open(newline="") as handle:
    blectin_comparison = list(csv.DictReader(handle, delimiter="\t"))
if len(blectin_comparison) != 28 or len({r["bacterial_accession"] for r in blectin_comparison}) != 1:
    raise SystemExit("B-lectin GC3 comparator is incomplete")
numbers["gc3"]["blectin_bacterial_comparator_accession"] = blectin_comparison[0]["bacterial_accession"]
with args.chromosome_gain_summary.open(newline="", encoding="utf-8") as handle:
    gain_rows = list(csv.DictReader(handle, delimiter="\t"))
if len(gain_rows) != 4 or {r["mode"] for r in gain_rows} != {"deltran"} or len({r["family"] for r in gain_rows}) != 4:
    raise SystemExit("incomplete eCAMI chromosome-gain summary")
with args.chromosome_gain_table.open(newline="", encoding="utf-8") as handle:
    gain_table = {(r["family_or_stratum"], r["mode"]): r for r in csv.DictReader(handle, delimiter="\t")
                  if r["section"] == "A"}
if len(gain_table) != 4 or {mode for _, mode in gain_table} != {"deltran"}:
    raise SystemExit("incomplete reader gain table lineage split")
for row in gain_rows:
    reader = gain_table[row["family"], row["mode"]]
    if any(int(reader[key]) != int(row[value]) for key, value in
           (("gains", "gains"), ("cynipini_gains", "cynipini_gains"),
            ("other_lineage_gains", "elsewhere_gains"))):
        raise SystemExit(f"gain lineage split differs between reader table and summary: {row['family']} {row['mode']}")
    if int(reader["gains"]) != int(reader["cynipini_gains"]) + int(reader["other_lineage_gains"]):
        raise SystemExit(f"reader gain split does not close: {row['family']} {row['mode']}")
numbers["ecami_chromosome_gain"] = {
    "eligible_native_groups": int(gain_rows[0]["all_native_groups"]),
    "families": {
        family: {
            "root_present": next(r for r in gain_rows if r["family"] == family)["root_present"],
            "root_ambiguous": next(r for r in gain_rows if r["family"] == family)["root_ambiguous"],
            "median_copies": float(next(r for r in gain_rows if r["family"] == family)["median_copies"]),
            "copy_stratum": next(r for r in gain_rows if r["family"] == family)["copy_stratum"],
            "own_stratum_groups": int(next(r for r in gain_rows if r["family"] == family)["own_stratum_groups"]),
            "cynipini_gains": int(gain_table[family, "deltran"]["cynipini_gains"]),
            "elsewhere_gains": int(gain_table[family, "deltran"]["other_lineage_gains"]),
            "modes": {
            mode: {
                "gains": int(r["gains"]), "losses": int(r["losses"]),
                "cynipini_gains": int(gain_table[family, mode]["cynipini_gains"]),
                "elsewhere_gains": int(gain_table[family, mode]["other_lineage_gains"]),
                "placed_tips": int(r["placed_tips"]),
                "gains_per_placed_tip": float(r["gains_per_placed_tip"]),
                "all_native_midrank_percentile": float(r["all_midrank_percentile"]) if r["all_midrank_percentile"] else None,
                "own_stratum_midrank_percentile": float(r["own_midrank_percentile"]) if r["own_midrank_percentile"] else None,
            }
            for mode in ("deltran",)
            for r in gain_rows if r["family"] == family and r["mode"] == mode
            },
        }
        for family in ("GH5_2", "PL1_8", "PL4_4", "B-lectin")
    },
}
with args.cazy_rank_table.open(newline="", encoding="utf-8") as handle:
    rank_rows = list(csv.DictReader(handle, delimiter="\t"))
expected_rank_keys = {
    (group, family)
    for group in ("Cynipini", "non-Cynipini Cynipidae", "Diplolepididae", "E. primus", "other Figitidae")
    for family in ("GH5_2", "PL1_8", "PL4_4")
}
if len(rank_rows) != 15 or {(r["lineage_group"], r["family"]) for r in rank_rows} != expected_rank_keys:
    raise SystemExit("condensed CAZy rank table is incomplete")
numbers["ecami_copy_rank_by_lineage"] = {
    r["lineage_group"]: {
        f["family"]: {"genomes": int(f["genomes_n"]),
                      "present": int(f["genomes_with_family_present"]),
                      "median_rank": float(f["median_rank"]) if f["median_rank"] else None,
                      "median_midrank_percentile": float(f["median_midrank_percentile"]) if f["median_midrank_percentile"] else None,
                      "rank_range": f["rank_range"]}
        for f in rank_rows if f["lineage_group"] == r["lineage_group"]}
    for r in rank_rows}
with args.synteny_summary.open(newline="", encoding="utf-8") as handle:
    summary_rows = list(csv.DictReader(handle, delimiter="\t"))
if not summary_rows or list(summary_rows[0]) != ["metric", "family", "block_id", "chromosome", "value"]:
    raise SystemExit("Invalid synteny summary schema")
summary = {(r["metric"], r["family"], r["block_id"]): int(r["value"]) for r in summary_rows}
if len(summary) != len(summary_rows):
    raise SystemExit("Duplicate synteny summary key")
def fact(metric):
    return summary[(metric, "all", "")]
numbers["synteny_patterns"] = fact("family_specific_blocks")
numbers["synteny_cross_genome_patterns"] = fact("cross_genome_detectable_patterns")
numbers["synteny_patterns_targets_two_or_more_genomes"] = fact("patterns_targets_two_or_more_genomes")
numbers["synteny_coverage"] = {key: fact(metric) for key, metric in (
    ("covered", "covered_copies"), ("accepted", "accepted_copies"))}
numbers["synteny_maximum_targets_per_locus"] = fact("maximum_targets_per_locus")
numbers["synteny_close_PL1_8_blocks"] = fact("close_PL1_8_blocks")
numbers["synteny_T_remus_recovered"] = {
    "recovered": fact("T_remus_recovered_cross_genome"),
    "cross_genome_total": fact("T_remus_cross_genome_total")}
modal = {}
for row in summary_rows:
    if row["metric"].startswith("modal_"):
        key = f"{row['family']}:{row['chromosome']}:{row['block_id']}"
        modal.setdefault(key, {})[row["metric"].removeprefix("modal_")] = int(row["value"])
numbers["synteny_modal_blocks"] = modal
with args.synteny_table.open(newline="", encoding="utf-8") as handle:
    table_rows = list(csv.DictReader(handle, delimiter="\t"))
part_a = [r for r in table_rows if r["section"] == "A"]
part_b = [r for r in table_rows if r["section"] == "B"]
if sum(int(r["eligible_target_genes"]) for r in part_a) != fact("accepted_copies"):
    raise SystemExit("Synteny summary differs from Table S4 roster")
table_modal = {(r["family_id"], r["modal_locus_id"], r["detected_genomes"], r["occupied_genomes"], r["empty_genomes"]) for r in part_b}
summary_modal = {(r["family"], r["chromosome"], str(modal[f"{r['family']}:{r['chromosome']}:{r['block_id']}"]["detected"]),
                  str(modal[f"{r['family']}:{r['chromosome']}:{r['block_id']}"]["occupied"]),
                  str(modal[f"{r['family']}:{r['chromosome']}:{r['block_id']}"]["empty"]))
                 for r in summary_rows if r["metric"] == "modal_detected"}
if not table_modal <= summary_modal:
    raise SystemExit("Synteny summary differs from Table S4 modal blocks")
with args.origin_key_nodes.open(newline="") as handle:
    key_rows = list(csv.DictReader(handle, delimiter="\t"))
if len(key_rows) != 12:
    raise SystemExit("expected twelve origin key-node rows")
numbers["origin_key_nodes"] = {
    family: {
        row["node"]: {"support": row["support"], "sister_taxa": row["sister_taxa"],
                      "sister_tip_count": row["sister_tip_count"]}
        for row in key_rows if row["family"] == family and row["view"] == "fig2"
    }
    for family in ("GH5_2", "PL1_8", "PL4_4")
}
if any(set(nodes) != {"cynipoid_clade", "join"} for nodes in numbers["origin_key_nodes"].values()):
    raise SystemExit("missing origin key nodes")
for n in (9,):
    with (args.new_tables / f"Table_S{n}.tsv").open(newline="") as handle:
        numbers[f"table_s{n}"] = list(csv.DictReader(handle, delimiter="\t"))
# Supplement prose facts (authority/supplement/supplement_source.md tokens).
numbers_panel = supplement_numbers.panel_facts(args.source_panel)
numbers_screen = supplement_numbers.screen_facts(args.data_root / "repository_inputs/supplement_facts", args.domain_summary)
numbers["supplement"] = {
    **numbers_panel,
    **supplement_numbers.expression_facts(args.expression_table),
    **supplement_numbers.sequence_set_facts(args.family_membership, args.sequence_sets, args.table_s2,
                                            args.b_lectin_manifest),
    **supplement_numbers.catalytic_domain_facts(args.domain_summary),
    **supplement_numbers.bassettia_attribution_facts(
        args.data_root / "repository_inputs/deposit/base/sidecars/removed_cb_origin_bp_tips.tsv"),
    **supplement_numbers.e_primus_facts(args.data_root / "repository_inputs/supplement_facts"),
    **supplement_numbers.later_proteome_facts(args.data_root / "repository_inputs/supplement_facts"),
    **supplement_numbers.cynog_catalog_facts(args.data_root / "repository_inputs/supplement_facts"),
    **numbers_screen,
    **supplement_numbers.cazy_lineage_group_facts(Path(__file__).resolve().parents[2] / "config/cazy_lineage_groups.tsv"),
    "synteny_covered_pct": round(100 * numbers["synteny_coverage"]["covered"] / numbers["synteny_coverage"]["accepted"], 1),
    "gc3": supplement_numbers.gc3_facts(numbers["gc3"]),
    "author_intron_sites": {key: {"pct": 100 * value["observations"] / numbers["introns"]["table_s5_part_a"][key.split(":")[0]]["eligible_proteins"]} for key, value in numbers["introns"]["table_s5_sites"].items()},
    "author_intron_chr": {f"{r['family']}:{r['ref_aa_pos']}:{chromosome}": count for r in author_intron_summary for chromosome, count in json.loads(r["locus_class_breakdown"]).items()},
    "author_blectin_intron_pct": 100 * numbers["introns"]["table_s5_part_a"]["B_lectin"]["proteins_with_introns"] / numbers["introns"]["table_s5_part_a"]["B_lectin"]["eligible_proteins"],
    "author_background_median_pct": 100 * numbers["introns"]["whole_genome_comparator"]["median_relative_cds_position"],
    "author_alternate_assemblies": numbers_screen["tblastn_assemblies"] - numbers_panel["genomes_considered"],
    "author_excluded_panel": numbers_panel["genomes_considered"] - numbers_panel["retained_genomes"],
    "intron_background_target_genes": sum(v["eligible_proteins"] for v in numbers["introns"]["table_s5_part_a"].values()),
    "background_first_quarter_pct": 100 * numbers["introns"]["whole_genome_comparator"]["first_cds_quarter_fraction"],
    "target_first_quarter_pct": 100 * numbers["introns"]["first_cds_quarter"]["introns"] / numbers["introns"]["first_cds_quarter"]["total"],
}
if numbers["supplement"]["PL1_8"]["displayed_domains"] != origin_totals["PL1_8"] or \
        numbers["supplement"]["PL4_4"]["displayed_domains"] != origin_totals["PL4_4"]:
    raise SystemExit("supplement displayed domain sets differ from Table S2 origin totals")
if numbers["supplement"]["e_primus_busco_complete_pct"] != numbers["supplement"]["e_primus_genome_busco_pct"]:
    raise SystemExit("E. primus BUSCO report differs from the Table S1 source panel")
args.out.parent.mkdir(parents=True, exist_ok=True)
args.out.write_text(json.dumps(numbers, indent=2) + "\n")
print("NUMBERS_PASS " + json.dumps({
    k: (v.get("all", v) if isinstance(v, dict) else len(v) if isinstance(v, list) else v)
    for k, v in numbers.items()
}))
