#!/usr/bin/env python3
"""Validate and materialize the checksum-pinned final-49 synteny authority.

The public producer starts at the accepted published-intermediate boundary; it
does not rerun the upstream private genome/GFF/CynOG construction route.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import shutil
from collections import Counter
from pathlib import Path

import pandas as pd


MANIFEST_COLUMNS = [
    "role", "path", "bytes", "sha256", "schema", "source_authority", "recovery_source",
]
ROLE_TO_OUTPUT = {
    "final_49_roster": "FINAL_49_SYNTENY_ROSTER.tsv",
    "block_summary": "table_s06_synteny_block_summary.tsv",
    "blocks_long": "table_s06_synteny_blocks_long.tsv",
    "locus_detail": "table_s06_synteny_locus_detail.tsv",
    "flanking_maps": "table_s06d_synteny_flanking_gene_maps.tsv",
    "region_state": "current_region_status_complete.tsv",
    "region_state_qc": "current_region_status_qc.tsv",
    "manuscript_statistics": "MANUSCRIPT_SYNTENY_STATISTICS.tsv",
    "tandem_summary": "current_tandem_locus_summary.tsv",
    "close_block_sensitivity": "close_block_definition_sensitivity.tsv",
}
ROLE_TO_REBUILD_OUTPUT = {
    "excluded_non_genomic_units": "excluded_non_genomic_units.tsv",
    "beta_builder_mapping": "beta_builder_target_mapping.tsv",
    "beta_anchor_audit": "beta_anchor_multiplicity_audit.tsv",
}
EXPECTED_ROLES = set(ROLE_TO_OUTPUT) | set(ROLE_TO_REBUILD_OUTPUT) | {
    "ordered_cynog_dictionary", "target_members", "genome_panel",
    "final49_validators", "deterministic_rerun_qc",
}
SCHEMAS = {
    "ordered_cynog_dictionary": ["genome", "chrom", "position", "protein_id", "cynog_id", "start_bp", "end_bp", "strand"],
    "target_members": ["protein_id", "family"],
    "genome_panel": ["accession", "species", "family", "tribe", "superfamily", "phylo_order", "tree_order", "status", "busco_complete", "busco_single", "busco_dup", "busco_frag", "busco_missing", "busco_status", "busco_reason", "n50_mb", "genome_mb", "repeat_families"],
    "final_49_roster": ["accession_or_project", "species", "family", "builder_genome_key", "roster_source"],
    "block_summary": ["family", "block_id", "block_type", "loci", "species_count", "multi_species", "species", "accepted_locus_ids", "old_packaged_relation"],
    "blocks_long": ["family", "accepted_locus_id", "species", "source_type", "record_class", "builder_target_id", "builder_locus_id", "block_id", "assignment_status", "block_status", "chromosome", "coordinate_subject", "coordinate_start", "coordinate_end", "placement_state", "projection_class", "unassigned_reason", "accepted_aliases"],
    "locus_detail": ["family", "accepted_locus_id", "species", "source_type", "record_class", "builder_target_id", "builder_locus_id", "block_id", "assignment_status", "block_status", "chromosome", "coordinate_subject", "coordinate_start", "coordinate_end", "placement_state", "projection_class", "unassigned_reason", "accepted_aliases", "position_class", "copy_status", "genome_has_core", "j_vs_core", "n_conserved_cynogs", "n_cynogs", "block_size"],
    "flanking_maps": ["family", "block_id", "accepted_locus_id", "species", "scaffold"] + [f"U{i}" for i in range(12, 0, -1)] + [f"D{i}" for i in range(1, 13)],
    "region_state": ["family", "block_id", "genome", "species", "status", "recovered", "scaffold", "region_start_rank", "region_end_rank", "n_target_units", "target_unit_ids", "n_conserved_flanking_cynogs", "search_route", "inapplicable_reason"],
    "region_state_qc": ["check", "observed", "expected", "status"],
    "manuscript_statistics": ["authority", "statistic", "value", "unit_definition"],
    "tandem_summary": ["builder_locus_id", "family", "species", "genome", "scaffold", "chromosome", "rank_start", "rank_end", "block_id", "n_builder_targets", "n_unique_builder_targets", "collapsed_duplicate_builder_rows", "builder_target_ids", "unique_builder_target_ids", "accepted_target_count", "accepted_target_ids", "multiple_targets"],
    "close_block_sensitivity": ["distance_rule", "counting_unit", "count", "pair_scope", "distance_definition"],
    "excluded_non_genomic_units": ["accepted_locus_id", "species", "source_type", "record_class", "exclusion_reason"],
    "beta_builder_mapping": ["internal_target_id", "biological_unit_id", "species", "species_code", "coordinate_subject", "coordinate_start", "coordinate_end", "projection_class", "anchor_type", "anchor_rows", "genome_key"],
    "beta_anchor_audit": ["internal_target_id", "accepted_biological_id", "accepted_species", "projected_anchor_rows", "declared_anchor_rows", "projected_genome_key", "current_panel_genome_key", "panel_key_rebound", "relationship"],
    "final49_validators": ["check", "observed", "expected", "status"],
    "deterministic_rerun_qc": ["check", "status"],
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-manifest", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def write_tsv(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    fields = fields or list(rows[0])
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def load_manifest(path: Path) -> dict[str, dict]:
    if not path.is_file():
        raise RuntimeError(f"Input manifest does not exist: {path}")
    rows = read_tsv(path)
    if list(rows.columns) != MANIFEST_COLUMNS:
        raise RuntimeError(f"Input manifest schema changed: {list(rows.columns)}")
    if rows.role.duplicated().any() or set(rows.role) != EXPECTED_ROLES:
        raise RuntimeError(
            f"Input role set changed: missing={sorted(EXPECTED_ROLES - set(rows.role))} "
            f"extra={sorted(set(rows.role) - EXPECTED_ROLES)}"
        )
    root = path.resolve().parent
    result: dict[str, dict] = {}
    for row in rows.to_dict("records"):
        relative = Path(row["path"])
        if relative.is_absolute() or ".." in relative.parts:
            raise RuntimeError(f"Manifest path must be relative and traversal-free: {relative}")
        source = (root / relative).resolve()
        if root not in source.parents or not source.is_file():
            raise RuntimeError(f"Manifest input is absent or outside its component: {relative}")
        if source.stat().st_size != int(row["bytes"]) or digest(source) != row["sha256"]:
            raise RuntimeError(f"Checksum/size mismatch for input role {row['role']}")
        frame = read_tsv(source)
        if row["schema"] != row["role"] + "_v1" or list(frame.columns) != SCHEMAS[row["role"]]:
            raise RuntimeError(f"Schema mismatch for input role {row['role']}")
        result[row["role"]] = {**row, "source": source, "frame": frame}
    return result


def validate(inputs: dict[str, dict]) -> list[dict]:
    frame = lambda role: inputs[role]["frame"]
    panel, roster = frame("genome_panel"), frame("final_49_roster")
    targets, gene_order = frame("target_members"), frame("ordered_cynog_dictionary")
    blocks, long, detail = frame("block_summary"), frame("blocks_long"), frame("locus_detail")
    flanks, matrix, stats = frame("flanking_maps"), frame("region_state"), frame("manuscript_statistics")
    excluded = frame("excluded_non_genomic_units")
    beta_mapping, beta_audit = frame("beta_builder_mapping"), frame("beta_anchor_audit")

    def require(condition: bool, message: str) -> None:
        if not condition:
            raise RuntimeError(message)

    require(len(panel) == 49 and panel.accession.nunique() == panel.species.nunique() == 49, "Panel is not 49 unique genomes")
    require(set(panel.status) == {"active"}, "Panel contains a non-active genome")
    require(len(roster) == 49 and roster.builder_genome_key.nunique() == roster.species.nunique() == 49, "Roster is not 49 unique genomes")
    require(set(roster.builder_genome_key) == set(panel.accession) and set(roster.species) == set(panel.species), "Panel and roster disagree")
    expected_targets = Counter({"GH5_2": 148, "PL1_8": 196, "PL4_4": 180, "CBM13": 42})
    require(len(targets) == targets.protein_id.nunique() == 566 and Counter(targets.family) == expected_targets, "Target membership changed")
    target_order = gene_order.loc[gene_order.protein_id.isin(set(targets.protein_id))]
    require(set(target_order.protein_id) == set(targets.protein_id), "Ordered CynOG dictionary does not cover all 566 targets")
    require(len(target_order) == 583, "Ordered CynOG target-row multiplicity changed")
    require(target_order.protein_id.str.startswith("BETAUNIT").sum() == 44, "B-lectin ordered-anchor count changed")
    require(set(panel.accession).issubset(set(gene_order.genome)), "A panel genome is absent from the ordered CynOG dictionary")
    require(len(long) == len(detail) == 570 and long.accepted_locus_id.nunique() == detail.accepted_locus_id.nunique() == 570, "Eligible locus contract changed")
    require(Counter(long.assignment_status) == Counter({"assigned": 566, "unassigned": 4}), "Assignment totals changed")
    require(set(long.record_class) == {"genome_locus"} and set(long.source_type) == {"genome"}, "Non-genomic locus entered public rows")
    require(set(long.accepted_locus_id) == set(detail.accepted_locus_id), "Long/detail locus identities disagree")
    assigned = long.loc[long.assignment_status.eq("assigned")]
    require(len(flanks) == 566 and set(flanks.accepted_locus_id) == set(assigned.accepted_locus_id), "Flanking maps do not cover assigned loci")
    require(blocks[["family", "block_id"]].drop_duplicates().shape[0] == 116, "Block count changed")
    require(len(matrix) == 5684 and matrix[["family", "block_id", "genome"]].drop_duplicates().shape[0] == len(matrix), "State matrix is not 116 x 49 unique cells")
    require(set(matrix.genome) == set(panel.accession), "State matrix panel changed")
    require(Counter(matrix.status) == Counter({"PRESENT": 309, "EMPTY": 2223, "N_D": 3152}), "State totals changed")
    serialized = {token for value in matrix.target_unit_ids for token in value.split(";") if token}
    require(serialized == set(assigned.accepted_locus_id), "Matrix membership differs from assigned loci")
    require(not any(value.startswith("BETAUNIT") for value in serialized), "Temporary BETAUNIT entered serialized public membership")
    final_stats = stats.loc[stats.authority.eq("final49")].set_index("statistic").value.astype(int).to_dict()
    expected_stats = {
        "total_family_specific_blocks": 116, "patterns_detectable_more_than_one_genome": 84,
        "patterns_targets_two_or_more_genomes": 27, "target_coverage_numerator": 455,
        "target_coverage_denominator": 566, "close_block_cases": 6,
        "T_remus_recovered_all_blocks_numerator": 27, "T_remus_recovered_all_blocks_denominator": 116,
        "T_remus_recovered_cross_genome_numerator": 27, "T_remus_recovered_cross_genome_denominator": 84,
        "modal_matrix_PRESENT": 309, "modal_matrix_EMPTY": 2223, "modal_matrix_N_D": 3152,
    }
    require(all(final_stats.get(key) == value for key, value in expected_stats.items()), "Final-49 statistics changed")
    require(len(excluded) == 22 and set(excluded.record_class) == {"historical_transcript"}, "Transcript-only exclusion audit changed")
    require(len(beta_mapping) == beta_mapping.biological_unit_id.nunique() == 46, "B-lectin unit mapping changed")
    require(len(beta_audit) == beta_audit.accepted_biological_id.nunique() == 42, "B-lectin accepted-unit audit changed")
    require(beta_audit.projected_anchor_rows.astype(int).sum() == 44, "B-lectin anchor count changed")
    require(Counter(beta_audit.panel_key_rebound) == Counter({"no": 33, "yes": 9}), "B-lectin rebinding count changed")
    require(Counter(beta_audit.relationship) == Counter({"one_anchor_to_one_unit": 40, "many_anchors_to_one_unit": 2}), "B-lectin multiplicity changed")
    require(set(beta_audit.accepted_biological_id) == set(assigned.loc[assigned.family.eq("beta_lectin/GH30_lineage"), "accepted_locus_id"]), "B-lectin accepted IDs disagree")
    require(frame("final49_validators").status.eq("PASS").all(), "Accepted final-49 validator is not PASS")
    require(frame("deterministic_rerun_qc").status.eq("PASS").all(), "Accepted deterministic replay receipt is not PASS")
    require(frame("region_state_qc").status.eq("PASS").all(), "Accepted region-state QC is not PASS")
    displays = pd.concat([roster.species, panel.species, long.species, matrix.species], ignore_index=True).astype(str)
    forbidden = {"Dryocosmus_cannoni", "Ceroptres_masudai", "Eschatocerus_acacia", "Isocolus_centaureae"}
    require(not set(displays).intersection(forbidden), "Deprecated/forbidden display taxon entered the final interface")
    values = [("final_genomes", 49), ("eligible_loci", 570), ("assigned_loci", 566),
              ("unassigned_loci", 4), ("blocks", 116), ("state_cells", 5684),
              ("beta_units", 42), ("beta_anchors", 44), ("beta_rebindings", 9)]
    return [{"check": key, "observed": value, "expected": value, "status": "PASS"} for key, value in values]


def main() -> None:
    args = parse_args()
    output = args.output_dir.resolve()
    if output.exists():
        raise RuntimeError(f"Output directory already exists; refusing overwrite: {output}")
    inputs = load_manifest(args.input_manifest)
    qc = validate(inputs)
    output.mkdir(parents=True)
    rebuild = output / "current_rebuild"
    rebuild.mkdir()
    for role, name in ROLE_TO_OUTPUT.items():
        shutil.copyfile(inputs[role]["source"], output / name)
    for role, name in ROLE_TO_REBUILD_OUTPUT.items():
        shutil.copyfile(inputs[role]["source"], rebuild / name)
    shutil.copyfile(inputs["final49_validators"]["source"], output / "FINAL49_VALIDATORS.tsv")
    shutil.copyfile(inputs["deterministic_rerun_qc"]["source"], output / "DETERMINISTIC_RERUN_QC.tsv")
    write_tsv(output / "RECONCILIATION_QC.tsv", qc)
    receipt = [
        {"kind": "input", "role": role, "path": row["path"], "bytes": row["bytes"], "sha256": row["sha256"]}
        for role, row in sorted(inputs.items())
    ]
    for path in sorted(output.rglob("*")):
        if path.is_file() and path.name != "RECONCILIATION_RECEIPT.tsv":
            receipt.append({"kind": "output", "role": path.stem, "path": path.relative_to(output).as_posix(),
                            "bytes": path.stat().st_size, "sha256": digest(path)})
    write_tsv(output / "RECONCILIATION_RECEIPT.tsv", receipt)
    print("PASS final49 reconciliation: 49 genomes; 570 loci; 116x49 states; 42 units/44 anchors")


if __name__ == "__main__":
    main()
