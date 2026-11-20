#!/usr/bin/env python3
"""Build Tables S3/S4 and public evidence from the accepted final-49 handoff."""

import csv
import hashlib
import json
import os
from collections import Counter, defaultdict
from pathlib import Path
import sys
import subprocess
import itertools
import re

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from copy_table import read_copies  # noqa: E402


ROOT = Path(__file__).resolve().parents[2]
with (ROOT / "config/table_numbers.tsv").open(newline="") as handle:
    TABLE_LABELS = {r["table_id"]:r["label"] for r in csv.DictReader(handle,delimiter="\t")}
DATA = Path(os.environ.get("PCWDE_REPOSITORY_INPUTS", ROOT / "figshare_data/repository_inputs")).resolve()
BUILD = Path(os.environ.get("PCWDE_BUILD_ROOT", ROOT / "build"))
AUTH = DATA / "tables/authority/synteny_current"
COMPONENT = BUILD / "tables_s3_s4"
OUT = COMPONENT / "candidate"
STAGING = OUT
QC = COMPONENT / "qc"
RECEIPTS = COMPONENT / "receipts"
FAMILY_ORDER = ("GH5_2", "PL1_8", "PL4_4", "B-lectin")
FAMILY_MAP = {"beta_lectin/GH30_lineage": "B-lectin"}
COPY_TABLE = Path(os.environ["PCWDE_COPY_TABLE"]).resolve()
CORRECTED = DATA / "tables/authority/synteny_corrected"
with (ROOT / "config/synteny_derived_pins.tsv").open(newline="") as handle:
    DERIVED_SHA = {r["name"]: r["sha256"] for r in csv.DictReader(handle, delimiter="\t")}
if set(DERIVED_SHA) != {"blocks_long.tsv", "block_summary.tsv", "current_region_status_complete.tsv"}:
    raise SystemExit("Synteny derived pins are incomplete")
CORRECTED_SHA = {
    "membership": DERIVED_SHA["blocks_long.tsv"],
    "blocks": DERIVED_SHA["block_summary.tsv"],
    "states": DERIVED_SHA["current_region_status_complete.tsv"],
    "order": "255046ed9b38cb038c733ffaa0194a409b274259bb5282478780a3f735ccfee8",
    "affected": "b79704a78b361221fcfbd2e9e1ca02dd3c967da1525eae997226428935e63694",
    "labels": "1a57f0122154db7ced8459658501a32c509c70444608e1dbcfabf6570bf2ac5a",
}


def read_tsv(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def pinned(path, digest):
    if sha(path) != digest:
        raise SystemExit(f"Corrected synteny input hash mismatch: {path}")
    return path


membership_path = pinned(CORRECTED / "blocks_long.tsv", CORRECTED_SHA["membership"])
blocks_path = pinned(CORRECTED / "block_summary.tsv", CORRECTED_SHA["blocks"])
states_path = pinned(CORRECTED / "current_region_status_complete.tsv", CORRECTED_SHA["states"])
order_path = pinned(CORRECTED / "gene_order_final.tsv", CORRECTED_SHA["order"])
affected_path = pinned(CORRECTED / "RECONSTRUCTED_AFFECTED_ORDER.tsv", CORRECTED_SHA["affected"])
labels_path = pinned(CORRECTED / "gca_to_species.tsv", CORRECTED_SHA["labels"])


def stable(prefix, value):
    return f"{prefix}-{hashlib.sha256(value.encode()).hexdigest()[:12].upper()}"


def public_family(value):
    return FAMILY_MAP.get(value, value)


def initialize_output_tree():
    """Create the complete component tree once, before its first write."""
    if COMPONENT.exists():
        raise SystemExit(f"refusing to overwrite existing S3/S4 component root: {COMPONENT}")
    for path in (
        STAGING / "data/synteny",
        STAGING / "supplementary_tables",
        QC,
        RECEIPTS,
    ):
        path.mkdir(parents=True, exist_ok=False)


source_panel = read_tsv(DATA / "tables/authority/tables_s3_s4_current/source_panel.tsv")
species_lookup = {row["species_display"].replace(" ", "_"): (row["species_id"], row["species_display"]) for row in source_panel}
roster = read_tsv(AUTH / "final_49_roster.tsv")
if len(roster) != 49 or len({row["builder_genome_key"] for row in roster}) != 49:
    raise SystemExit("Final synteny roster is not 49 unique genomes")
if any(row["species"].replace(" ", "_") not in species_lookup for row in roster):
    raise SystemExit("Final roster species do not close against public source-panel IDs")
roster_public = [(species_lookup[row["species"].replace(" ", "_")][0], species_lookup[row["species"].replace(" ", "_")][1]) for row in roster]
species_id_by_source = {row["species"].replace(" ", "_"): species_lookup[row["species"].replace(" ", "_")][0] for row in roster}
species_id_by_genome = {row["builder_genome_key"]: species_lookup[row["species"].replace(" ", "_")][0] for row in roster}

all_loci = read_tsv(membership_path)
if len(all_loci) != 675 or len({row["accepted_locus_id"] for row in all_loci}) != 675:
    raise SystemExit("Corrected locus authority is not 675 unique copies")

# Table S3 counts gene copies from the shared copy table: a copy on a BK
# chromosome (RagTag, JCVI or the BK reference) is named; the rest are unplaced.
placement_rows = []
for copy in read_copies(COPY_TABLE):
    placed = copy["chromosome"] != "unplaced"
    placement_rows.append({
        "placement_id": stable("PLC", copy["protein_id"]),
        "locus_id": copy["protein_id"],
        "species_id": species_id_by_source[copy["species"].replace(" ", "_")],
        "family_id": copy["family"],
        "placement_class": "direct_chromosome" if placed else "unplaced",
        "reference_block": copy["chromosome"] if placed else "",
        "source_scaffold": copy["scaffold"],
        "placement_confidence": copy["placement_route"] if placed else "unresolved",
    })
placement_rows.sort(key=lambda r: (r["species_id"], FAMILY_ORDER.index(r["family_id"]), r["locus_id"]))
synteny_ids = {row["accepted_locus_id"] for row in all_loci}
placement_ids = {row["locus_id"] for row in placement_rows}
rescue_id = "NODE_38372_length_2218_cov_4.820177:1323-346"
if placement_ids - synteny_ids != {rescue_id} or synteny_ids - placement_ids:
    raise SystemExit("Table S3 placement differs from synteny cohort beyond the unplaced rescue")
rescue_placement = [r for r in placement_rows if r["locus_id"] == rescue_id]
if len(rescue_placement) != 1 or rescue_placement[0]["placement_class"] != "unplaced":
    raise SystemExit("rescue placement changed")
identity = read_tsv(DATA / "members/synteny_existing_identity.tsv")
if len(identity) != len(all_loci) or {r["accepted_id"] for r in identity} != synteny_ids:
    raise SystemExit("Synteny identity disposition does not cover the accepted genome copies")
eligible_ids = {r["accepted_id"] for r in identity if r["order_id"]}
loci = [r for r in all_loci if r["accepted_locus_id"] in eligible_ids]
if len(loci) != len(eligible_ids) or any(r["assignment_status"] != "assigned" for r in loci):
    raise SystemExit("Eligible copies do not close against native block assignments")

placement_counts = Counter(row["placement_class"] for row in placement_rows)
expected_placement = Counter({"direct_chromosome": 543, "unplaced": 133})
if placement_counts != expected_placement:
    raise SystemExit(f"Placement closure changed: {placement_counts}")
table_s4 = []
grouped_placement = defaultdict(Counter)
for row in placement_rows:
    grouped_placement[(row["species_id"], row["family_id"])][row["placement_class"]] += 1
for species_id, display in roster_public:
    for family in FAMILY_ORDER:
        counts = grouped_placement[(species_id, family)]
        total = sum(counts.values())
        table_s4.append({
            "species_id": species_id,
            "species_display": display,
            "family_id": family,
            "accepted_target_count": total,
            "named_chromosome_count": counts["direct_chromosome"],
            "unplaced_count": counts["unplaced"],
        })
placement_total = len(placement_rows)
if placement_total != 676:
    raise SystemExit("Copy placement cohort is not 676 copies")
if len(table_s4) != 196 or sum(int(row["accepted_target_count"]) for row in table_s4) != placement_total:
    raise SystemExit("Table S4 species-by-family closure failed")

initialize_output_tree()
write_tsv(STAGING / "data/chromosome_placement_evidence.tsv", list(placement_rows[0]), placement_rows)
write_tsv(STAGING / "supplementary_tables/Table_S3.tsv", list(table_s4[0]), table_s4)
caption = (
    f"Table {TABLE_LABELS['loci']}. Chromosome placement of {placement_total} accepted genome-derived target copies in the final {len(roster)}-genome panel. "
    "Named-chromosome counts use RagTag, JCVI or the *B. kinseyi* reference chromosome; remaining copies are unplaced. "
    f"Synteny blocks are summarized separately in Table {TABLE_LABELS['synteny']}.\n"
)
(STAGING / "supplementary_tables/Table_S3_caption.txt").write_text(caption, encoding="utf-8")

# Stable public block and locus-assignment projections.
source_blocks = read_tsv(blocks_path)
if not source_blocks or len({(r["family"], r["block_id"]) for r in source_blocks}) != len(source_blocks):
    raise SystemExit("Derived synteny block keys are empty or duplicate")
block_id_map = {(row["family"], row["block_id"]): stable("SYN", row["family"] + "|" + row["block_id"]) for row in source_blocks}
source_states = read_tsv(states_path)
if len(source_states) != len(source_blocks) * 49:
    raise SystemExit("Derived block-state authority is not blocks x 49 genomes")

state_map = {"PRESENT": "occupied", "EMPTY": "empty", "N_D": "not_detected"}
public_states = []
state_counts_by_block = defaultdict(Counter)
for row in source_states:
    key = (row["family"], row["block_id"])
    state = state_map[row["status"]]
    state_counts_by_block[key][state] += 1
    public_states.append({
        "block_id": block_id_map[key],
        "species_id": species_id_by_genome[row["genome"]],
        "state": state,
        "target_locus_ids": row["target_unit_ids"],
    })
if len({(row["block_id"], row["species_id"]) for row in public_states}) != len(source_states):
    raise SystemExit("Block-by-genome keys are not unique")

public_blocks = []
locus_by_id = {row["accepted_locus_id"]: row for row in loci}
reference_name_aliases = {"Eschatocerus acacia": "Eschatocerus acaciae",
                          "Isocolus centaureae": "Eubothrus centaureae"}
for row in source_blocks:
    key = (row["family"], row["block_id"])
    counts = state_counts_by_block[key]
    loci_in_block = sorted(filter(None, row["accepted_locus_ids"].split(";")))
    reference = locus_by_id[loci_in_block[0]] if loci_in_block else None
    description = (reference_name_aliases.get(reference["species"].replace("_", " "),
                                              reference["species"].replace("_", " ")) + "; " +
                   ("BK " + reference["chromosome"] if reference["chromosome"] else "BK chromosome undetermined")) if reference else ""
    public_blocks.append({
        "block_id": block_id_map[key],
        "family_id": public_family(row["family"]),
        "reference_locus_id": loci_in_block[0] if loci_in_block else "",
        "reference_locus_description": description,
        "detected_genomes": counts["occupied"] + counts["empty"],
        "occupied_genomes": counts["occupied"],
        "empty_genomes": counts["empty"],
        "not_detected_genomes": counts["not_detected"],
    })
public_blocks.sort(key=lambda r: (FAMILY_ORDER.index(r["family_id"]), r["block_id"]))

public_assignments = []
for row in loci:
    assigned = row["assignment_status"] == "assigned"
    if not assigned:
        raise SystemExit(f"Eligible copy has no block: {row['accepted_locus_id']}")
    public_assignments.append({
        "locus_id": row["accepted_locus_id"],
        "block_id": block_id_map[(row["family"], row["block_id"])] if assigned else "",
        "assignment_status": "assigned",
        "evidence_summary": row["projection_class"].replace("_", " "),
    })
public_assignments.sort(key=lambda r: r["locus_id"])
if len(public_assignments) != len(eligible_ids):
    raise SystemExit("Public assignment coverage differs from eligible copies")

# Expand the producer's U12..U1/D1..D12 compact maps into the public long schema.
flank_path = BUILD / "corrected_synteny/flanking_gene_maps_source.tsv"
subprocess.run([sys.executable, str(ROOT / "scripts/analysis/synteny/flanking_maps.py"),
    "--order", str(order_path), "--order-sha256", CORRECTED_SHA["order"],
    "--affected-order", str(affected_path), "--affected-order-sha256", CORRECTED_SHA["affected"],
    "--membership", str(membership_path), "--membership-sha256", CORRECTED_SHA["membership"],
    "--historical", str(AUTH / "flanking_gene_maps_source.tsv"),
    "--species-labels", str(labels_path), "--out", str(flank_path),
    "--report", str(BUILD / "corrected_synteny/flank_regression.tsv")], check=True)
wide_flanks = read_tsv(flank_path)
flank_columns = [f"U{i}" for i in range(12, 0, -1)] + [f"D{i}" for i in range(1, 13)]
public_flanks = []
for row in wide_flanks:
    if row["accepted_locus_id"] not in eligible_ids:
        continue
    public_block = block_id_map[(row["family"], row["block_id"])]
    for position in flank_columns:
        if not row[position]:
            continue
        key = f"{public_block}|{row['accepted_locus_id']}|{position}|{row[position]}"
        public_flanks.append({
            "flank_map_id": stable("FLK", key),
            "block_id": public_block,
            "species_id": species_id_by_source[row["species"]],
            "anchor_gene_id": f"{row['accepted_locus_id']}:{position}",
            "matched_gene_id": row[position],
            "orientation": "unresolved",
        })
if not public_flanks or len({row["flank_map_id"] for row in public_flanks}) != len(public_flanks):
    raise SystemExit("Expanded flanking map is empty or non-unique")

write_tsv(STAGING / "data/synteny/block_summary.tsv", list(public_blocks[0]), public_blocks)
write_tsv(STAGING / "data/synteny/locus_assignments.tsv", list(public_assignments[0]), public_assignments)
write_tsv(STAGING / "data/synteny/flanking_gene_maps.tsv", list(public_flanks[0]), public_flanks)
write_tsv(STAGING / "data/synteny/block_by_genome_states.tsv", list(public_states[0]), public_states)

blocks_by_family = defaultdict(list)
assignments_by_family = defaultdict(list)
for block in public_blocks:
    blocks_by_family[block["family_id"]].append(block)
assignment_by_id = {row["locus_id"]: row for row in public_assignments}
if len(assignment_by_id) != len(eligible_ids) or set(assignment_by_id) != eligible_ids:
    raise SystemExit("Public assignments do not pair one-to-one with eligible copies")
if {row["accepted_locus_id"] for row in loci} != eligible_ids:
    raise SystemExit("Derived source has missing or duplicate eligible IDs")
for source in loci:
    assignments_by_family[public_family(source["family"])].append(assignment_by_id[source["accepted_locus_id"]])
expected_family_assigned = Counter(row["family"] for row in identity if row["order_id"])
if {family: sum(row["assignment_status"] == "assigned" for row in rows)
    for family, rows in assignments_by_family.items()} != expected_family_assigned:
    raise SystemExit("ID-joined family assignment counts changed")

sys.path.insert(0, str(ROOT / "scripts/analysis/synteny"))
from comprehensive_coverage import comprehensive_coverage
coverage_evidence, shared_totals, shared_covered, comprehensive_multi = comprehensive_coverage(
    placement_rows, all_loci, source_states, species_id_by_genome)
if shared_totals != Counter({family: sum(int(r["accepted_target_count"]) for r in table_s4
                                        if r["family_id"] == family) for family in FAMILY_ORDER}):
    raise SystemExit("Shared S4/S5 family totals disagree")
write_tsv(COMPONENT / "comprehensive_copy_coverage.tsv", list(coverage_evidence[0]), coverage_evidence)

table_s5 = []
for family in FAMILY_ORDER:
    blocks = blocks_by_family[family]
    assignments = assignments_by_family[family]
    assigned = sum(row["assignment_status"] == "assigned" for row in assignments)
    table_s5.append({
        "section": "A", "family_id": family,
        "eligible_target_genes": shared_totals[family],
        "detectable_patterns": len(blocks),
        "cross_genome_patterns": sum(int(row["detected_genomes"]) > 1 for row in blocks),
        "target_positive_patterns": sum(int(row["occupied_genomes"]) > 0 for row in blocks),
        "modal_locus_id": "", "detected_genomes": "", "occupied_genomes": "", "empty_genomes": "", "not_detected_genomes": "",
        "block_id": "", "reference_locus": "", "recovery_breadth": "", "target_positive_breadth": "",
    })

def modal_chromosome(family, block_id):
    chromosomes = Counter(row["chromosome"] for row in loci
                          if row["family"] == family and row["block_id"] == block_id
                          and row["assignment_status"] == "assigned"
                          and row["chromosome"] not in ("", "fragment"))
    if not chromosomes:
        raise SystemExit(f"No placed chromosome for modal block {family}/{block_id}")
    highest = max(chromosomes.values())
    winners = sorted(chromosome for chromosome, count in chromosomes.items() if count == highest)
    if len(winners) != 1:
        raise SystemExit(f"Tied modal chromosomes for {family}/{block_id}: {winners}")
    return winners[0]


# Apply one occupancy-first selection rule to every family. Chromosome-specific
# rows use only placed blocks on that chromosome; modal_chromosome rejects ties.
def select_modal_block(family, chromosome=None):
    candidates = [row for row in source_blocks
                  if public_family(row["family"]) == family
                  and row["block_type"] == "synteny_block"]
    if chromosome is not None:
        placed_blocks = {(row["family"], row["block_id"]) for row in loci
                         if row["assignment_status"] == "assigned"
                         and row["chromosome"] not in ("", "fragment")}
        candidates = [row for row in candidates
                      if (row["family"], row["block_id"]) in placed_blocks
                      and modal_chromosome(row["family"], row["block_id"]) == chromosome]
    if not candidates:
        raise SystemExit(f"No modal block candidates for {family}/{chromosome}")
    def rank(row):
        counts = state_counts_by_block[(row["family"], row["block_id"])]
        return (-counts["occupied"], -(counts["occupied"] + counts["empty"]), row["block_id"])
    chosen = min(candidates, key=rank)
    return chosen["family"], chosen["block_id"]

modal_blocks = [select_modal_block(family, chromosome)
                for family, chromosome in (("GH5_2", "chr10"), ("GH5_2", "chr6"),
                                           ("PL1_8", None), ("PL4_4", None),
                                           ("B-lectin", "chr3"), ("B-lectin", "chr6"))]
if {public_family(family) for family, _ in modal_blocks} != set(FAMILY_ORDER):
    raise SystemExit("Part B must represent every family")
modal_source = []
for family, block_id in modal_blocks:
    modal = modal_chromosome(family, block_id)
    for source_state, public_state in (("PRESENT", "occupied"), ("EMPTY", "empty"), ("N_D", "not_detected")):
        value = state_counts_by_block[(family, block_id)][public_state]
        modal_source.append({"modal_locus": modal, "family": family, "block_id": block_id, "statistic": source_state, "value": str(value)})
modal_grouped = defaultdict(dict)
for row in modal_source:
    modal_grouped[(row["family"], row["modal_locus"], row["block_id"])][row["statistic"]] = row["value"]
for (family, modal, source_block), values in sorted(modal_grouped.items()):
    table_s5.append({
        "section": "B", "family_id": public_family(family),
        "eligible_target_genes": "", "detectable_patterns": "", "cross_genome_patterns": "", "target_positive_patterns": "",
        "modal_locus_id": modal, "detected_genomes": int(values["PRESENT"]) + int(values["EMPTY"]), "occupied_genomes": values["PRESENT"], "empty_genomes": values["EMPTY"], "not_detected_genomes": values["N_D"],
        "block_id": "", "reference_locus": "", "recovery_breadth": "", "target_positive_breadth": "",
    })
summary_fields = ("blocks_one_genome", "blocks_2_to_24_genomes", "blocks_25_plus_genomes")
for row in table_s5:
    row.update({field: "" for field in summary_fields})
for family in FAMILY_ORDER:
    blocks = [r for r in public_blocks if r["family_id"] == family]
    counts = [int(r["detected_genomes"]) for r in blocks]
    table_s5.append({**{field: "" for field in table_s5[0]},
                     "section": "C", "family_id": family,
                     "blocks_one_genome": sum(n == 1 for n in counts),
                     "blocks_2_to_24_genomes": sum(2 <= n <= 24 for n in counts),
                     "blocks_25_plus_genomes": sum(n >= 25 for n in counts),
                     "target_positive_patterns": sum(int(r["occupied_genomes"]) > 0 for r in blocks)})
if any(not row["section"] or not row["family_id"] for row in table_s5):
    raise SystemExit("Table S5 row has an empty section or family")
modal_keys = [(row["family_id"], row["modal_locus_id"]) for row in table_s5 if row["section"] == "B"]
if len(modal_keys) != len(set(modal_keys)):
    raise SystemExit("Table S5 Part B has duplicate family/modal locus rows")
write_tsv(STAGING / "supplementary_tables/Table_S4.tsv", list(table_s5[0]), table_s5)

# The Results counts are one build intermediate, derived from the same six
# registered inputs as the table and reader workbook.
summary = []
def add(metric, value, family="all", block="", chromosome=""):
    summary.append({"metric": metric, "family": family, "block_id": block,
                    "chromosome": chromosome, "value": str(value)})

assigned = [r for r in loci if r["assignment_status"] == "assigned"]
multi = comprehensive_multi
covered = sum(shared_covered.values())
add("family_specific_blocks", len(source_blocks))
add("cross_genome_detectable_patterns", sum(int(r["detected_genomes"]) > 1 for r in public_blocks))
add("patterns_targets_two_or_more_genomes", len(multi))
add("accepted_copies", sum(shared_totals.values()))
add("covered_copies", covered)
max_tandem = max(Counter((r["family"], r["builder_locus_id"]) for r in assigned).values())
add("maximum_targets_per_locus", max_tandem)

# Distinct PL1_8 blocks on the same ordered scaffold, fewer than 20 intervening
# annotated genes. The builder locus encodes the rank interval in gene_order_final.
intervals = defaultdict(list)
for row in assigned:
    if row["family"] != "PL1_8":
        continue
    match = re.search(r"_pos(\d+)(?:_to_(\d+))?$", row["builder_locus_id"])
    if not match:
        raise SystemExit(f"Missing ordered rank: {row['builder_locus_id']}")
    intervals[(row["species"], row["coordinate_subject"], row["block_id"])].append(
        (int(match[1]), int(match[2] or match[1])))
by_scaffold = defaultdict(list)
for (species, scaffold, block), spans in intervals.items():
    by_scaffold[(species, scaffold)].append((block, min(x[0] for x in spans), max(x[1] for x in spans)))
close = set()
for (species, scaffold), spans in by_scaffold.items():
    for a, b in itertools.combinations(spans, 2):
        if a[0] == b[0]:
            continue
        left, right = sorted((a, b), key=lambda x: x[1])
        if 0 <= right[1] - left[2] - 1 < 20:
            close.add((species, scaffold, tuple(sorted((a[0], b[0])))))
add("close_PL1_8_blocks", len(close))

remus_genomes = {r["genome"] for r in source_states if r["species"] == "Telenomus_remus"}
if len(remus_genomes) != 1:
    raise SystemExit(f"Expected one Telenomus remus genome, got {remus_genomes}")
remus = next(iter(remus_genomes))
cross = {(r["family"], r["block_id"]) for r in source_blocks
         if state_counts_by_block[(r["family"], r["block_id"])]["occupied"] +
            state_counts_by_block[(r["family"], r["block_id"])]["empty"] > 1}
add("T_remus_recovered_cross_genome", sum(
    (r["family"], r["block_id"]) in cross and r["status"] != "N_D"
    for r in source_states if r["genome"] == remus))
add("T_remus_cross_genome_total", len(cross))

for family, block in modal_blocks:
    chromosome = modal_chromosome(family, block)
    counts = state_counts_by_block[(family, block)]
    for metric, value in (("modal_detected", counts["occupied"] + counts["empty"]),
                          ("modal_occupied", counts["occupied"]),
                          ("modal_empty", counts["empty"]),
                          ("modal_not_detected", counts["not_detected"])):
        add(metric, value, public_family(family), block, chromosome)
if len(assigned) != len(eligible_ids) or len(source_states) != len(source_blocks) * len(roster):
    raise SystemExit("Eligible synteny summary does not close against roster and blocks")
write_tsv(COMPONENT / "synteny_summary.tsv", list(summary[0]), summary)
block_counts = []
for block in source_blocks:
    family, block_id = block["family"], block["block_id"]
    counts = state_counts_by_block[(family, block_id)]
    detected = counts["occupied"] + counts["empty"]
    chromosomes = Counter(row["chromosome"] for row in loci
                          if row["family"] == family and row["block_id"] == block_id
                          and row["assignment_status"] == "assigned"
                          and row["chromosome"] not in ("", "fragment"))
    chromosome = sorted(chromosomes, key=lambda value: (-chromosomes[value], value))[0] if chromosomes else ""
    block_counts.append({
        "block_id": block_id, "family": family,
        "chromosome": chromosome,
        "detected_genomes": detected, "occupied_genomes": counts["occupied"],
        "empty_genomes": counts["empty"], "not_detected_genomes": counts["not_detected"],
        "cross_genome": "yes" if detected > 1 else "no",
        "targets_two_or_more_genomes": "yes" if (public_family(family), block_id) in multi else "no",
    })
write_tsv(COMPONENT / "synteny_block_counts.tsv", list(block_counts[0]), block_counts)
(STAGING / "supplementary_tables/Table_S4_caption.txt").write_text(
    f"Table {TABLE_LABELS['synteny']}. Conserved-locus recovery for {sum(shared_totals.values())} accepted target copies across {len(roster)} unique genomes. Part A summarizes family-level target assignment and block recovery; Part B reports selected reference chromosomes and counts of genomes; Part C counts blocks by recovery breadth (one genome, 2–24 genomes, or at least 25 genomes). Reference chromosome is the most frequent BK-equivalent chromosome among placed copies in a block. The full {len(public_blocks)}-block list is in the data deposit. Target absent from recovered region means the block was detected without an accepted target; region not recovered is insufficient evidence of absence.\n",
    encoding="utf-8",
)

public_outputs = [
    STAGING / "data/chromosome_placement_evidence.tsv",
    STAGING / "supplementary_tables/Table_S3.tsv",
    STAGING / "supplementary_tables/Table_S3_caption.txt",
    STAGING / "data/synteny/block_summary.tsv",
    STAGING / "data/synteny/locus_assignments.tsv",
    STAGING / "data/synteny/flanking_gene_maps.tsv",
    STAGING / "data/synteny/block_by_genome_states.tsv",
    STAGING / "supplementary_tables/Table_S4.tsv",
    STAGING / "supplementary_tables/Table_S4_caption.txt",
]
lower = b"\n".join(path.read_bytes().lower() for path in public_outputs)
for token in (b"dcan", b"dryocosmus cannoni", b"beta lectin", b"/" + b"carc/", b"/" + b"easley/"):
    if token in lower:
        raise SystemExit(f"Forbidden public token: {token.decode()}")

qc4 = [
    ("evidence_rows", len(placement_rows), len(placement_rows) == placement_total),
    ("species_family_rows", len(table_s4), len(table_s4) == 196),
    ("direct_chromosome", placement_counts["direct_chromosome"], placement_counts["direct_chromosome"] == expected_placement["direct_chromosome"]),
    ("unplaced", placement_counts["unplaced"], placement_counts["unplaced"] == expected_placement["unplaced"]),
    ("target_closure", sum(int(r["accepted_target_count"]) for r in table_s4), sum(int(r["accepted_target_count"]) for r in table_s4) == placement_total),
    ("forbidden_tokens", 0, True),
]
qc5 = [
    ("blocks", len(public_blocks), bool(public_blocks) and {r["family_id"] for r in public_blocks} == set(FAMILY_ORDER)),
    ("eligible_loci", len(public_assignments), len(public_assignments) == len(eligible_ids)),
    ("assigned_loci", sum(r["assignment_status"] == "assigned" for r in public_assignments), all(r["assignment_status"] == "assigned" for r in public_assignments)),
    ("state_cells", len(public_states), len(public_states) == len(public_blocks) * len(roster)),
    ("state_partition", len(public_states), sum(r["state"] in {"occupied", "empty", "not_detected"} for r in public_states) == len(public_states)),
    ("forbidden_tokens", 0, True),
]
for path, checks in ((QC / "Table_S3_validation.tsv", qc4), (QC / "Table_S4_validation.tsv", qc5)):
    write_tsv(path, ["check", "value", "status"], [{"check": k, "value": v, "status": "PASS" if ok else "FAIL"} for k, v, ok in checks])
    if not all(ok for _, _, ok in checks):
        raise SystemExit(f"QC failed: {path.name}")

for artifact, paths, counts in (
    ("ART-TABLE-S3", public_outputs[:3], {"accepted_targets": placement_total, "species_family_rows": 196}),
    ("ART-TABLE-S4", public_outputs[3:], {"blocks": len(public_blocks), "state_cells": len(public_states), "eligible_loci": len(public_assignments)}),
):
    receipt = {
        "artifact_id": artifact,
        "status": "qc_passed",
        "input_set": "repository_inputs/tables/authority/synteny_current",
        "inputs": {str(path.relative_to(DATA)): sha(path) for path in
                   (membership_path, blocks_path, states_path, order_path, affected_path,
                    labels_path, AUTH / "final_49_roster.tsv",
                    DATA / "tables/authority/tables_s3_s4_current/source_panel.tsv")} |
                  {"built:members/copy_placements.tsv": sha(COPY_TABLE)},
        "outputs": {str(path.relative_to(OUT)): sha(path) for path in paths},
        "counts": counts,
    }
    (RECEIPTS / f"{artifact}.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
print(json.dumps({"table_s4_rows": len(table_s4), "table_s5_rows": len(table_s5), "flank_rows": len(public_flanks)}, sort_keys=True))
