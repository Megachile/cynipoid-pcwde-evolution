#!/usr/bin/env python3
"""Validate, describe, and stage the clean final-49/no-DCAN synteny release."""
from __future__ import annotations

import argparse
import csv
import hashlib
import shutil
from collections import Counter
from pathlib import Path

import pandas as pd

def read(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def write(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    fields = fields or list(rows[0])
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reconciliation-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--build-receipt", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    WORK = args.reconciliation_dir.resolve()
    BUILD = WORK / "current_rebuild"
    STAGE = args.output_dir.resolve()
    RECEIPT = args.build_receipt.resolve()
    if not WORK.is_dir():
        raise RuntimeError(f"Reconciliation directory does not exist: {WORK}")
    if STAGE.exists():
        raise RuntimeError(f"Output directory already exists; refusing overwrite: {STAGE}")
    if RECEIPT.exists():
        raise RuntimeError(f"Build receipt already exists; refusing overwrite: {RECEIPT}")
    reconciliation_qc = read(WORK / "RECONCILIATION_QC.tsv")
    if list(reconciliation_qc.columns) != ["check", "observed", "expected", "status"] or not reconciliation_qc.status.eq("PASS").all():
        raise RuntimeError("Reconciliation QC is absent, malformed, or not PASS")
    reconciliation_receipt = read(WORK / "RECONCILIATION_RECEIPT.tsv")
    if list(reconciliation_receipt.columns) != ["kind", "role", "path", "bytes", "sha256"]:
        raise RuntimeError("Reconciliation receipt schema changed")
    for row in reconciliation_receipt.loc[reconciliation_receipt.kind.eq("output")].itertuples(index=False):
        path = WORK / row.path
        if not path.is_file() or path.stat().st_size != int(row.bytes) or digest(path) != row.sha256:
            raise RuntimeError(f"Reconciliation output receipt mismatch: {row.path}")

    roster = read(WORK / "FINAL_49_SYNTENY_ROSTER.tsv")
    roster_species = set(roster.species.str.replace(" ", "_"))
    long = read(WORK / "table_s06_synteny_blocks_long.tsv")
    detail = read(WORK / "table_s06_synteny_locus_detail.tsv")
    matrix = read(WORK / "current_region_status_complete.tsv")
    blocks = read(WORK / "table_s06_synteny_block_summary.tsv")
    stats = read(WORK / "MANUSCRIPT_SYNTENY_STATISTICS.tsv")
    excluded = read(BUILD / "excluded_non_genomic_units.tsv")
    beta_audit = read(BUILD / "beta_anchor_multiplicity_audit.tsv")

    if len(roster) != 49 or roster.species.nunique() != 49 or roster.builder_genome_key.nunique() != 49:
        raise RuntimeError("Roster is not 49 unique biological genomes")
    if len(long) != 570 or len(detail) != 570:
        raise RuntimeError(f"Expected 570 genome-derived eligible loci, observed long={len(long)} detail={len(detail)}")
    if set(long.record_class) != {"genome_locus"} or set(long.source_type) != {"genome"}:
        raise RuntimeError("Non-genomic rows entered public locus tables")
    if not set(long.species).issubset(roster_species):
        raise RuntimeError("Public locus table contains a non-roster species")
    if long.accepted_locus_id.nunique() != 570 or detail.accepted_locus_id.nunique() != 570:
        raise RuntimeError("Accepted genome-derived locus identifiers are not unique")
    status_counts = Counter(long.assignment_status)
    if status_counts != Counter({"assigned": 566, "unassigned": 4}):
        raise RuntimeError(f"Eligible locus assignment contract changed: {status_counts}")
    if len(excluded) != 22 or set(excluded.record_class) != {"historical_transcript"}:
        raise RuntimeError("Transcript-only source-edge exclusion audit is incomplete")
    beta_long = long.loc[
        long.family.eq("beta_lectin/GH30_lineage") & long.assignment_status.eq("assigned")
    ]
    if (
        len(beta_audit) != 42
        or beta_audit.accepted_biological_id.nunique() != 42
        or beta_audit.projected_anchor_rows.astype(int).sum() != 44
        or Counter(beta_audit.relationship) != Counter({"one_anchor_to_one_unit": 40, "many_anchors_to_one_unit": 2})
        or Counter(beta_audit.panel_key_rebound) != Counter({"no": 33, "yes": 9})
        or set(beta_audit.accepted_biological_id) != set(beta_long.accepted_locus_id)
    ):
        raise RuntimeError("B-lectin 44-anchor/42-accepted-unit audit failed")
    nblocks = blocks[["family", "block_id"]].drop_duplicates().shape[0]
    if nblocks != 116 or len(matrix) != 116 * 49:
        raise RuntimeError("Block-by-genome matrix is not 116 x 49")
    if matrix[["family", "block_id", "genome"]].drop_duplicates().shape[0] != len(matrix):
        raise RuntimeError("Block-by-genome matrix has duplicate cells")
    if set(matrix.genome) != set(roster.builder_genome_key):
        raise RuntimeError("Matrix genome roster differs from final-49 roster")
    if set(matrix.status) != {"PRESENT", "EMPTY", "N_D"} or len(matrix) != sum(Counter(matrix.status).values()):
        raise RuntimeError(f"Block-state matrix has invalid state coverage: {Counter(matrix.status)}")

    assigned_ids = set(long.loc[long.assignment_status.eq("assigned"), "accepted_locus_id"])
    matrix_ids = {
        target
        for serialized in matrix.target_unit_ids
        for target in serialized.split(";")
        if target
    }
    if matrix_ids != assigned_ids:
        raise RuntimeError(
            "Block-by-genome target membership differs from assigned accepted loci: "
            f"missing={len(assigned_ids - matrix_ids)} extra={len(matrix_ids - assigned_ids)}"
        )
    serialized_counts = matrix.target_unit_ids.map(lambda value: len([x for x in value.split(";") if x]))
    if not serialized_counts.eq(matrix.n_target_units.astype(int)).all():
        raise RuntimeError("Block-by-genome target counts disagree with serialized target identifiers")

    final_stats = stats.loc[stats.authority.eq("final49")].set_index("statistic").value.astype(int).to_dict()
    assigned = long.loc[long.assignment_status.eq("assigned")]
    multi_keys = set(zip(
        blocks.loc[blocks.species_count.astype(int).ge(2), "family"],
        blocks.loc[blocks.species_count.astype(int).ge(2), "block_id"],
    ))
    coverage = sum(key in multi_keys for key in zip(assigned.family, assigned.block_id))
    required_stats = {
        "total_family_specific_blocks", "patterns_detectable_more_than_one_genome",
        "patterns_targets_two_or_more_genomes", "target_coverage_numerator",
        "target_coverage_denominator", "close_block_cases",
        "T_remus_recovered_cross_genome_numerator", "T_remus_recovered_cross_genome_denominator",
    }
    if not required_stats.issubset(final_stats):
        raise RuntimeError(f"Missing recalculated manuscript statistics: {sorted(required_stats - set(final_stats))}")
    if final_stats["target_coverage_numerator"] != coverage or final_stats["target_coverage_denominator"] != len(assigned):
        raise RuntimeError("Recalculated target coverage disagrees with accepted locus rows")
    directed_stats = {
        "total_family_specific_blocks": 116,
        "patterns_detectable_more_than_one_genome": 84,
        "patterns_targets_two_or_more_genomes": 27,
        "target_coverage_numerator": 455,
        "target_coverage_denominator": 566,
        "close_block_cases": 6,
        "T_remus_recovered_all_blocks_numerator": 27,
        "T_remus_recovered_all_blocks_denominator": 116,
        "T_remus_recovered_cross_genome_numerator": 27,
        "T_remus_recovered_cross_genome_denominator": 84,
        "modal_matrix_PRESENT": 309,
        "modal_matrix_EMPTY": 2223,
        "modal_matrix_N_D": 3152,
    }
    mismatches = {
        key: (final_stats.get(key), expected)
        for key, expected in directed_stats.items()
        if final_stats.get(key) != expected
    }
    if mismatches:
        raise RuntimeError(
            f"Recalculated final-49 statistics differ from the authorized package contract: {mismatches}"
        )
    if Counter(matrix.status) != Counter({"PRESENT": 309, "EMPTY": 2223, "N_D": 3152}):
        raise RuntimeError(f"Modal matrix totals changed: {Counter(matrix.status)}")

    family_rows = []
    for family, group in blocks.groupby("family", sort=True):
        family_rows.append({
            "family": family,
            "family_specific_blocks": group.block_id.nunique(),
            "accepted_assigned_loci": int(assigned.family.eq(family).sum()),
            "multi_species_target_blocks": int(group.species_count.astype(int).ge(2).sum()),
        })
    STAGE.mkdir(parents=True)
    write(STAGE / "family_synteny_statistics.tsv", family_rows)

    modal_rows = []
    for block_id, chromosome, held_species in [("BLOCK_0001", "chr6", 22), ("BLOCK_0000", "chr10", 26)]:
        members = assigned.loc[assigned.family.eq("GH5_2") & assigned.block_id.eq(block_id)]
        placed_chromosomes = members.loc[members.chromosome.ne("fragment"), "chromosome"]
        if placed_chromosomes.empty or placed_chromosomes.mode().tolist() != [chromosome]:
            raise RuntimeError(f"{block_id} does not retain expected modal {chromosome} identity")
        if members.species.nunique() != held_species:
            raise RuntimeError(f"{block_id} species count changed")
        states = Counter(matrix.loc[matrix.family.eq("GH5_2") & matrix.block_id.eq(block_id), "status"])
        for state in ["PRESENT", "EMPTY", "N_D"]:
            modal_rows.append({"modal_locus": chromosome, "family": "GH5_2", "block_id": block_id,
                               "statistic": state, "value": states[state], "unit": "retained genomes"})
        modal_rows.append({"modal_locus": chromosome, "family": "GH5_2", "block_id": block_id,
                           "statistic": "target_positive_species", "value": members.species.nunique(),
                           "unit": "retained species with accepted targets"})
    write(STAGE / "modal_GH5_2_locus_statistics.tsv", modal_rows)

    held = [
        ("total_family_specific_blocks", 117, final_stats["total_family_specific_blocks"], "family-specific blocks"),
        ("cross_genome_detectable_patterns", 83, final_stats["patterns_detectable_more_than_one_genome"], "blocks"),
        ("target_positive_multi_genome_patterns", 27, final_stats["patterns_targets_two_or_more_genomes"], "blocks"),
        ("eligible_target_coverage_numerator", 455, final_stats["target_coverage_numerator"], "genome-derived accepted loci in multi-genome blocks"),
        ("eligible_target_coverage_denominator", 570, final_stats["target_coverage_denominator"] + int(long.assignment_status.eq("unassigned").sum()), "all genome-derived accepted loci"),
        ("close_block_cases", 6, final_stats["close_block_cases"], "unique species/family/block-pair cases"),
        ("T_remus_cross_genome_recovery_numerator", 28, final_stats["T_remus_recovered_cross_genome_numerator"], "blocks"),
        ("T_remus_cross_genome_recovery_denominator", 83, final_stats["T_remus_recovered_cross_genome_denominator"], "blocks"),
    ]
    comparison = [
        {"statistic": key, "held_historical_value": old, "clean_final49_value": new,
         "delta": new - old, "status": "unchanged" if old == new else "changed", "counting_unit": unit}
        for key, old, new, unit in held
    ]
    write(STAGE / "held_constant_comparison.tsv", comparison)

    eligible_loci = len(long)
    assigned_loci = len(assigned)
    conditional_percent = 100 * coverage / assigned_loci
    eligible_percent = 100 * coverage / eligible_loci
    states = Counter(matrix.status)
    brief = f"""# Manuscript synteny brief — exact final 49 genomes

The rebuilt authority contains **{nblocks} family-specific blocks across exactly {len(roster)} unique current genomes**. Of these, {final_stats['patterns_detectable_more_than_one_genome']} patterns are recoverable as PRESENT or EMPTY in more than one genome, while {final_stats['patterns_targets_two_or_more_genomes']} contain accepted targets in at least two genomes. Those multi-genome target-bearing patterns account for {coverage}/{eligible_loci} eligible accepted biological target units ({eligible_percent:.1f}%). Among assigned targets, the corresponding conditional fraction is {coverage}/{assigned_loci} ({conditional_percent:.1f}%).

The complete modal-locus state matrix contains {states['PRESENT']:,} PRESENT, {states['EMPTY']:,} EMPTY, and {states['N_D']:,} N_D cells. The maximum is {final_stats['maximum_targets_per_locus']} accepted targets in one maintained builder locus. Under the contracted definition (distinct blocks, same family/species/exact ordered scaffold, fewer than 20 intervening annotated genes), there are {final_stats['close_block_cases']} unique species/family/block-pair cases.

For *Telenomus remus*, recovery is {final_stats['T_remus_recovered_all_blocks_numerator']}/{final_stats['T_remus_recovered_all_blocks_denominator']} across all blocks and {final_stats['T_remus_recovered_cross_genome_numerator']}/{final_stats['T_remus_recovered_cross_genome_denominator']} across the independently reproduced cross-genome-detectable subset.

The final-49 dataset contains 570 eligible genomic target units (566 assigned and four unassigned), the 42-biological-unit/44-anchor B-lectin relationship, and nine current-panel key rebindings.
"""
    (STAGE / "MANUSCRIPT_SYNTENY_BRIEF.md").write_text(brief)

    public = {
        "FINAL_49_SYNTENY_ROSTER.tsv": WORK / "FINAL_49_SYNTENY_ROSTER.tsv",
        "table_s05_synteny_block_summary.tsv": WORK / "table_s06_synteny_block_summary.tsv",
        "table_s05_synteny_blocks_long.tsv": WORK / "table_s06_synteny_blocks_long.tsv",
        "table_s05_synteny_locus_detail.tsv": WORK / "table_s06_synteny_locus_detail.tsv",
        "table_s05d_synteny_flanking_gene_maps.tsv": WORK / "table_s06d_synteny_flanking_gene_maps.tsv",
        "block_by_genome_state.tsv": WORK / "current_region_status_complete.tsv",
        "manuscript_synteny_statistics.tsv": WORK / "MANUSCRIPT_SYNTENY_STATISTICS.tsv",
        "modal_GH5_2_locus_statistics.tsv": STAGE / "modal_GH5_2_locus_statistics.tsv",
        "held_constant_comparison.tsv": STAGE / "held_constant_comparison.tsv",
        "family_synteny_statistics.tsv": STAGE / "family_synteny_statistics.tsv",
    }
    for name, source in public.items():
        destination = STAGE / name
        if source.resolve() != destination.resolve():
            shutil.copyfile(source, destination)
    staged_roster = read(STAGE / "FINAL_49_SYNTENY_ROSTER.tsv")
    if "roster_source" in staged_roster.columns:
        staged_roster["roster_source"] = "final-49 retained manifest authority"
        staged_roster.to_csv(
            STAGE / "FINAL_49_SYNTENY_ROSTER.tsv", sep="\t", index=False, lineterminator="\n"
        )

    definitions = {
        "family": "Target gene family or lineage.", "block_id": "Maintained family-specific synteny block identifier.",
        "accepted_locus_id": "Unique accepted genome-derived biological target locus.", "species": "Current biological species label.",
        "source_type": "Sequence evidence source; public locus rows are genome-derived.", "record_class": "Accepted record unit; public locus rows are genome_locus.",
        "builder_target_id": "Target identifier supplied to the maintained synteny builder.", "builder_locus_id": "Physical target-locus identifier produced by the builder.",
        "assignment_status": "Whether the eligible locus could be assigned to a maintained synteny block.", "block_status": "Builder classification of the assigned locus/block.",
        "chromosome": "Chromosome label where available.", "coordinate_subject": "Assembly sequence carrying the locus.",
        "coordinate_start": "One-based accepted locus start when available.", "coordinate_end": "One-based accepted locus end when available.",
        "placement_state": "Coordinate placement status from the accepted locus authority.", "projection_class": "Provenance class for projection into maintained gene order.",
        "unassigned_reason": "Reason an eligible genome-derived locus has no block assignment.", "accepted_aliases": "Semicolon-separated accepted identifiers for the same locus.",
        "block_type": "Synteny block or true singleton.", "loci": "Number of accepted locus rows in the block.",
        "species_count": "Number of distinct retained species with accepted targets in the block.", "multi_species": "Whether targets occur in at least two retained species.",
        "accepted_locus_ids": "Semicolon-separated accepted loci in the block.", "old_packaged_relation": "Relation to the prior packaged block authority.",
        "position_class": "Maintained locus position classification.", "copy_status": "Maintained ancestral/extra-copy classification.",
        "genome_has_core": "Whether the genome contains the family core locus.", "j_vs_core": "Jaccard similarity to the family core signature.",
        "n_conserved_cynogs": "Count of conserved flanking CynOGs.", "n_cynogs": "Count of CynOGs in the locus signature.",
        "block_size": "Number of loci assigned to the builder block.", "scaffold": "Ordered assembly sequence used for the flanking map.",
        "genome": "Unique retained builder genome key.", "status": "Modal-locus state: PRESENT, EMPTY, or N_D.",
        "recovered": "Whether the modal locus was recovered as PRESENT or EMPTY.", "region_start_rank": "First ordered-gene rank in the recovered region.",
        "region_end_rank": "Last ordered-gene rank in the recovered region.", "n_target_units": "Accepted target loci in the block/genome cell.",
        "target_unit_ids": "Semicolon-separated target loci in the block/genome cell.", "n_conserved_flanking_cynogs": "Conserved flank count supporting recovery.",
        "search_route": "Recovery route used for the block/genome cell.", "inapplicable_reason": "Reason a state could not be assigned.",
        "accession_or_project": "Current assembly accession or project identifier.", "builder_genome_key": "Unique key used by the maintained synteny builder.",
        "roster_source": "Path to the retained-49 authority used for this build.", "authority": "Statistic authority/version.",
        "statistic": "Named statistic with its counting unit defined separately.", "value": "Integer statistic value.",
        "unit_definition": "Exact counting-unit definition.", "modal_locus": "Modal chromosome identity for the GH5_2 block.",
        "unit": "Counting unit for the statistic.", "held_historical_value": "Previously held manuscript/package value.",
        "clean_final49_value": "Value recalculated from this clean final-49 release.", "delta": "Clean value minus held value.",
        "counting_unit": "Exact unit counted in the comparison.",
        "accepted_biological_id": "Accepted biological B-lectin locus identifier represented by projected gene-model anchors.",
        "internal_target_id": "Temporary internal projection identifier retained only for multiplicity provenance and never serialized as public block membership.",
        "accepted_species": "Current accepted biological species for the projected B-lectin unit.",
        "projected_anchor_rows": "Number of retained gene-order anchor rows representing this accepted biological unit.",
        "declared_anchor_rows": "Anchor-row count recorded by the deterministic projection source.",
        "projected_genome_key": "Genome key originally emitted by the accepted projection authority.",
        "current_panel_genome_key": "Current final-panel key resolved through the accepted species/accession authority.",
        "panel_key_rebound": "Whether the projected key required rebinding to the current final-panel key.",
        "relationship": "Explicit one- or many-anchor relationship to one accepted biological unit.",
        "family_specific_blocks": "Recalculated number of distinct family-specific synteny blocks.",
        "accepted_assigned_loci": "Recalculated number of assigned accepted biological loci in the family.",
        "multi_species_target_blocks": "Recalculated family blocks with accepted targets in at least two species.",
    }
    dictionary_rows = []
    for name in public:
        columns = list(read(STAGE / name).columns)
        for column in columns:
            if column.startswith("U") or column.startswith("D"):
                definition = "Ordered upstream/downstream flanking protein and CynOG cell."
            else:
                definition = definitions.get(column)
            if not definition:
                raise RuntimeError(f"Missing public-field definition for {name}:{column}")
            dictionary_rows.append({"file": name, "field": column, "definition": definition})
    write(STAGE / "PUBLIC_FIELD_DATA_DICTIONARY.tsv", dictionary_rows)
    for path in STAGE.iterdir():
        if path.is_file() and (b"DCAN" in path.read_bytes() or b"Dryocosmus_cannoni" in path.read_bytes()):
            raise RuntimeError(f"Forbidden DCAN token in public candidate {path.name}")
        if path.is_file() and b"/" + b"carc/" in path.read_bytes():
            raise RuntimeError(f"Private CARC path in public candidate {path.name}")
    checksum_rows = [
        {"file": path.name, "bytes": path.stat().st_size, "sha256": digest(path)}
        for path in sorted(STAGE.iterdir()) if path.is_file()
    ]
    write(STAGE / "SHA256SUMS.tsv", checksum_rows)
    receipt = []
    for path in sorted(STAGE.iterdir()):
        if path.is_file():
            source = "deterministic input copy" if path.name in public else "generated metadata"
            receipt.append({"file": path.name, "source": source, "bytes": path.stat().st_size,
                            "sha256": digest(path), "build_scope": "49 unique retained genomes; genome-derived loci only"})
    write(RECEIPT, receipt)
    print("PASS final49 dataset: 116x49; eligible loci=570")


if __name__ == "__main__":
    main()
