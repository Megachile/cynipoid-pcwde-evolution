#!/usr/bin/env python3
import csv
import gzip
import hashlib
import os
import re
import shutil
import sys
import tarfile
import subprocess
from pathlib import Path
from typing import Optional

sys.path.insert(0, str(Path(__file__).resolve().parent))
from data_root import resolve_data_root, validate_reader_inputs, validate_projected_inputs, write_projected_inputs, validate_logical_inputs
from synteny_description import workbook_description


# These files contain private provenance omitted from the supporting archive.
# Scientific rows and non-hash scientific cells retain their original order.
OMIT_PRIVATE_PROVENANCE = {
    "analysis/synteny/BUILD_RECEIPT.tsv",
    "analysis/synteny/INPUT_MANIFEST.tsv",
    "analysis/synteny/RESULT_MANIFEST.tsv",
    "analysis/synteny/results/SHA256SUMS.tsv",
    "tables/authority/synteny_current/producer_build_receipt.tsv",
    "tables/introns_current/source_manifest.tsv",
}
ROW_HASH_TSVS = {
    "INPUT_MANIFEST.tsv",
    "figures/INPUT_MANIFEST.tsv",
    "tables/authority/table_s2/b_lectin_accepted_members.tsv",
    "tables/authority/table_s2/beta_units.tsv",
    "tables/authority/table_s2/b_lectin_domain_manifest.tsv",
    "tables/authority/selection_current/alignment_qc.tsv",
    "tables/authority/synteny_current/producer_build_receipt.tsv",
    "deposit/base/sidecars/B_lectin_domain_manifest.tsv",
    "deposit/base/sidecars/B_lectin_architecture.tsv",
    "deposit/base/data/selection/selection_alignment_qc.tsv",
    "analysis/synteny/BUILD_RECEIPT.tsv",
    "analysis/synteny/RESULT_MANIFEST.tsv",
    "analysis/synteny/INPUT_MANIFEST.tsv",
    "analysis/synteny/results/SHA256SUMS.tsv",
    "figures/presence_counts/beta_units.tsv",
    "figures/b_lectin_tree/interface.tsv",
    "figures/b_lectin_domain_tree/architecture.tsv",
    "figures/b_lectin_domain_tree/domain_manifest.tsv",
    "figures/accessory_domain_trees/label_inputs/tip_crosswalk.tsv",
}
FIGURE_MANIFEST = "figures/INPUT_MANIFEST.tsv"
CDS_MAPPING = "tables/authority/selection_current/current_cds_mapping.tsv"
S4_SEQUENCE_DIGEST_TSVS = {"analysis/s4_architecture/full_sequence_mapping.tsv",
    "figures/accessory_domain_trees/architecture.tsv",
}
ACCEPTED_PRIOR_PROJECTION_SHA256 = "8c8ed54598751f6893b74a49f4612171f6b2915fb90f9265980e4bd6f68477e6"
PRIVATE_PATTERN = re.compile(r"/(?:carc|easley|users|root)/")
REMOVED_PATTERN = re.compile(r"dryocosmus[ _]cannoni|\bDCAN\b", re.IGNORECASE)
TEXT_SUFFIXES = {".tsv", ".txt", ".md", ".faa", ".fa", ".fasta", ".treefile", ".nwk", ".tree", ".csv", ".r", ".py", ".bash", ".sbatch", ".fna", ".gff", ".gff3"}
INTRON_GFF = re.compile(r"analysis/introns/sources/[^/]+/gff\.gff$")


def read_tsv(path: Path) -> tuple[list[str], list[list[str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t", strict=True)
        header = next(reader, None)
        rows = list(reader)
    if not header or len(header) != len(set(header)):
        raise SystemExit(f"unsafe TSV header: {path}")
    return header, rows


def write_tsv(path: Path, header: list[str], rows: list[list[str]]) -> None:
    if any(len(row) != len(header) for row in rows):
        raise SystemExit(f"refusing nonrectangular projected TSV: {path}")
    temporary = path.with_name(f".{path.name}.public_projection.tmp")
    if temporary.exists() or temporary.is_symlink():
        raise SystemExit(f"stale public-projection temporary file: {temporary}")
    original_mode = path.stat().st_mode & 0o777
    try:
        with temporary.open("x", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(header)
            writer.writerows(rows)
        os.chmod(temporary, original_mode)
        os.replace(temporary, path)
    finally:
        if temporary.exists() and not temporary.is_symlink():
            temporary.unlink()


def write_and_verify_tsv(path: Path, header: list[str], rows: list[list[str]]) -> None:
    """Serialize once and prove no row/cell loss, duplication, or reordering."""
    expected = (header, rows)
    write_tsv(path, header, rows)
    if read_tsv(path) != expected:
        raise SystemExit(f"TSV projection identity failure: {path}")


def hash_fields(header: list[str], relative: str = "") -> list[str]:
    # S4 sequence_sha256 is the scientific identity of a full-length protein,
    # not a row-level provenance checksum. Keep it in the public architecture.
    return [field for field in header if "sha256" in field.casefold()
            and not (relative in S4_SEQUENCE_DIGEST_TSVS and field == "sequence_sha256")]


def public_problem_set(inputs: Path) -> set[tuple[str, str]]:
    problems: set[tuple[str, str]] = set()
    for path in sorted(p for p in inputs.rglob("*") if p.is_file()):
        relative = path.relative_to(inputs).as_posix()
        if path.is_symlink() or Path(relative).is_absolute() or ".." in Path(relative).parts:
            problems.add(("unsafe_path", relative))
            continue
        if path.suffix.casefold() == ".tsv":
            header, rows = read_tsv(path)
            if any(len(row) != len(header) for row in rows):
                problems.add(("nonrectangular_tsv", relative))
            if hash_fields(header, relative):
                problems.add(("row_hash_header", relative))
        if path.suffix.casefold() in TEXT_SUFFIXES:
            with path.open(encoding="utf-8", errors="replace") as handle:
                for line in handle:
                    if PRIVATE_PATTERN.search(line):
                        problems.add(("private_absolute_path", relative))
                    if path.suffix.casefold() not in {".fna", ".gff", ".gff3"} and not relative.startswith("analysis/synteny/") and REMOVED_PATTERN.search(line):
                        problems.add(("removed_taxon_token", relative))
    return problems


def expected_preprojection_problems(inputs: Path) -> set[tuple[str, str]]:
    expected = {("removed_taxon_token", FIGURE_MANIFEST),
                ("private_absolute_path", CDS_MAPPING),
                ("removed_taxon_token", "tables/authority/synteny_current/producer_build_receipt.tsv"),
                ("private_absolute_path", "tables/introns_current/source_manifest.tsv")}
    expected.update(("row_hash_header", path) for path in ROW_HASH_TSVS)
    declared_private_gffs = {("private_absolute_path", relative)
                             for relative in declared_gff_paths(inputs)}
    expected.update(intron_expected_preprojection_problems(inputs, ROOT) & declared_private_gffs)
    return expected


def declared_gff_paths(inputs: Path) -> list[str]:
    header, rows = read_tsv(ROOT / "config/intron_sources.tsv")
    source_index = header.index("source_file")
    return [f"analysis/introns/sources/{row[source_index]}/gff.gff" for row in rows
            if (inputs / f"analysis/introns/sources/{row[source_index]}/gff.gff").is_file()]


def rebuild_input_inventory(inputs: Path) -> int:
    """Write the projected path/size inventory; file hashes live in Git/package manifests."""
    manifest = inputs / "INPUT_MANIFEST.tsv"
    rows = []
    for path in sorted(p for p in inputs.rglob("*") if p.is_file() and p != manifest):
        rows.append([
            path.relative_to(inputs).as_posix(),
            str(path.stat().st_size),
            "accepted_public_intermediate",
        ])
    write_tsv(manifest, ["relative_path", "size_bytes", "source_class"], rows)
    return len(rows) + 1


def refresh_input_inventory_sizes(inputs: Path) -> int:
    """Update byte sizes in an already-projected inventory, keeping row order and source classes."""
    manifest = inputs / "INPUT_MANIFEST.tsv"
    header, rows = read_tsv(manifest)
    if header != ["relative_path", "size_bytes", "source_class"]:
        raise SystemExit("projected repository inventory schema mismatch")
    refreshed = [[relative, str((inputs / relative).stat().st_size), source_class]
                 for relative, _size, source_class in rows]
    if refreshed != rows:
        write_tsv(manifest, header, refreshed)
    return len(rows) + 1


def validate_projected_inventory(inputs: Path) -> int:
    manifest = inputs / "INPUT_MANIFEST.tsv"
    header, rows = read_tsv(manifest)
    if header != ["relative_path", "size_bytes", "source_class"]:
        raise SystemExit("projected repository inventory schema mismatch")
    expected = {row[0] for row in rows}
    actual = {
        path.relative_to(inputs).as_posix()
        for path in inputs.rglob("*") if path.is_file() and path != manifest
    }
    if expected != actual or len(expected) != len(rows):
        raise SystemExit(
            f"projected repository inventory closure mismatch missing={sorted(actual-expected)[:5]} "
            f"extra={sorted(expected-actual)[:5]}"
        )
    for relative, size, source_class in rows:
        path = inputs / relative
        corrected = (relative.startswith("tables/authority/synteny_corrected/")
                     and source_class == "corrected_synteny_authority")
        frozen_cazy = (relative.startswith("cazy_inputs/")
                       and source_class == "frozen_cazy_scan_input")
        if path.stat().st_size != int(size) or not (source_class == "accepted_public_intermediate" or corrected or frozen_cazy):
            raise SystemExit(f"projected repository inventory mismatch: {relative}")
    return len(rows) + 1


def refresh_figure_inventory(inputs: Path) -> None:
    """Bind the public figure inventory to the fully projected file sizes.

    Public projection deliberately removes row-level SHA fields and repairs
    trailing empty TSV cells.  Those transformations change seven figure-input
    byte sizes after the source manifest has already admitted their original
    hashes.  Refresh every sibling size only after all member transformations,
    while proving exact inventory closure and preserving manifest row order and
    all non-size fields.
    """
    figure_root = inputs / "figures"
    manifest = inputs / FIGURE_MANIFEST
    header, rows = read_tsv(manifest)
    if header != ["component", "relative_path", "role", "size_bytes"]:
        raise SystemExit("projected figure inventory schema mismatch before size refresh")
    path_index = header.index("relative_path")
    size_index = header.index("size_bytes")
    declared = [row[path_index] for row in rows]
    if len(declared) != len(set(declared)):
        raise SystemExit("projected figure inventory contains duplicate paths")
    actual = {
        path.relative_to(figure_root).as_posix()
        for path in figure_root.rglob("*")
        if path.is_file() and path != manifest
    }
    if set(declared) != actual:
        raise SystemExit(
            f"projected figure inventory closure mismatch "
            f"missing={sorted(actual-set(declared))[:5]} "
            f"extra={sorted(set(declared)-actual)[:5]}"
        )
    refreshed = []
    for row in rows:
        path = figure_root / row[path_index]
        if path.is_symlink() or not path.is_file():
            raise SystemExit(f"projected figure inventory member is nonregular: {row[path_index]}")
        updated = list(row)
        updated[size_index] = str(path.stat().st_size)
        refreshed.append(updated)
    write_and_verify_tsv(manifest, header, refreshed)


def projection_spec(inputs: Path, logical_count: int) -> list[dict[str, str]]:
    rows = []
    for operation, relative, detail in [
        ("neutralize_excluded_source_metadata", FIGURE_MANIFEST, "replace eight exclusion-label cells; identifiers unchanged"),
        ("omit_private_provenance_field", CDS_MAPPING, "omit cds_path; retain ten nonprivate fields"),
    ]:
        _, data = read_tsv(inputs / relative)
        rows.append({"operation": operation, "target": relative, "rows": str(len(data)), "detail": detail})
    for relative in sorted(ROW_HASH_TSVS - OMIT_PRIVATE_PROVENANCE):
        _, data = read_tsv(inputs / relative)
        rows.append({"operation": "omit_row_hash_fields", "target": relative, "rows": str(len(data)), "detail": "scientific rows and non-hash fields retained"})
    for relative in sorted(OMIT_PRIVATE_PROVENANCE):
        rows.append({"operation": "omit_private_provenance_file", "target": relative, "rows": "not_applicable", "detail": "receipt or manifest only; no scientific payload rows"})
    for path in (inputs / relative for relative in declared_gff_paths(inputs)):
        relative = path.relative_to(inputs).as_posix()
        digest_state = hashlib.sha256()
        count = 0
        with path.open("rb") as handle:
            for line in handle:
                digest_state.update(line)
                count += line == b"# private source path omitted\n"
        rows.append({"operation": "neutralize_private_gff_header", "target": relative, "rows": str(count), "detail": f"feature lines unchanged; projected_sha256={digest_state.hexdigest()}"})
    rows.append({"operation": "bind_file_level_logical_manifest", "target": "config/logical_inputs.tsv", "rows": str(logical_count), "detail": "assembler emits one immutable file-level digest and byte count per public logical input"})
    return rows


def project_repository_inputs(inputs: Path, prior_projection: Optional[Path] = None) -> list[dict[str, str]]:
    problems = public_problem_set(inputs)
    omitted_present = {relative for relative in OMIT_PRIVATE_PROVENANCE if (inputs / relative).is_file()}
    # stage.py starts from the accepted public deposit, not from the historical
    # private preprojection tree.  Its pinned projection receipt certifies the
    # old removals; the current S4 architecture input still needs its scientific digest retained.
    present_s4 = {relative for relative in S4_SEQUENCE_DIGEST_TSVS if (inputs / relative).is_file()}
    if present_s4 and present_s4 != S4_SEQUENCE_DIGEST_TSVS:
        raise SystemExit(f"incomplete staged S4 architecture input: {sorted(present_s4)}")
    gff_paths = [inputs / relative for relative in declared_gff_paths(inputs)]
    projected_gffs = [path for path in gff_paths if any(line == b"# private source path omitted\n" for line in path.open("rb"))]
    if projected_gffs and problems:
        raise SystemExit("partially projected intron GFF sources")
    gff_problems = {("private_absolute_path", path.relative_to(inputs).as_posix())
                    for path in gff_paths if gff_private_header_count(path)}
    s4_problems = gff_problems
    if prior_projection is not None and prior_projection.is_file():
        if hashlib.sha256(prior_projection.read_bytes()).hexdigest() != ACCEPTED_PRIOR_PROJECTION_SHA256:
            raise SystemExit("staged public projection receipt content drift")
        prior_header, prior_rows = read_tsv(prior_projection)
        if prior_header != ["operation", "target", "rows", "detail"]:
            raise SystemExit("staged public projection receipt schema drift")
        prior = [dict(zip(prior_header, row)) for row in prior_rows]
        declared = {(row["operation"], row["target"]) for row in prior}
        if len(declared) != len(prior):
            raise SystemExit("staged public projection receipt duplicate identity")
        allowed_operations = {"pad_trailing_empty_cells", "neutralize_excluded_source_metadata",
                              "omit_private_provenance_field", "omit_row_hash_fields",
                              "omit_private_provenance_file", "bind_file_level_logical_manifest"}
        if any(row["operation"] not in allowed_operations for row in prior):
            raise SystemExit("staged public projection receipt unknown operation")
        required = {("omit_private_provenance_field", CDS_MAPPING),
                    ("neutralize_excluded_source_metadata", FIGURE_MANIFEST)}
        # The accepted public stage renamed several historical figure and
        # intron paths.  Their original targets remain in its sealed receipt.
        old_hashes = [row for row in prior if row["operation"] == "omit_row_hash_fields"]
        old_omissions = [row for row in prior if row["operation"] == "omit_private_provenance_file"]
        if (not required <= declared or len(old_hashes) != 15 or len(old_omissions) != 6
                or omitted_present or problems != s4_problems
                or (present_s4 and any("sequence_sha256" not in read_tsv(inputs / relative)[0]
                                       for relative in present_s4))):
            raise SystemExit(f"unsafe/stale staged public projection: missing_receipts={sorted(required-declared)[:8]} "
                             f"problems={sorted(problems)} omissions={sorted(omitted_present)}")
        validate_projected_inventory(inputs)
        for path in gff_paths:
            project_gff_header(path)
        if public_problem_set(inputs):
            raise SystemExit("staged GFF projection left forbidden content")
        # The staged inventory is already projected and carries curated
        # source classes; refresh sizes only, so an unchanged staged input
        # layer stays a verbatim copy for the closure check.
        logical_count = refresh_input_inventory_sizes(inputs)
        validate_projected_inventory(inputs)
        prior = [row for row in prior if row["operation"] != "bind_file_level_logical_manifest"]
        for row in projection_spec(inputs, logical_count):
            if row["operation"] == "neutralize_private_gff_header":
                prior.append(row)
        for relative in sorted(present_s4):
            prior.append({"operation": "retain_scientific_sequence_digest", "target": relative,
                          "rows": str(len(read_tsv(inputs / relative)[1])),
                          "detail": "sequence_sha256 identifies the mapped full-length protein"})
        prior.append({"operation": "bind_file_level_logical_manifest", "target": "config/logical_inputs.tsv",
                      "rows": str(logical_count), "detail": "current projected inventory"})
        return prior
    if not problems and not omitted_present:
        logical_count = validate_projected_inventory(inputs)
        if gff_paths:
            return projection_spec(inputs, logical_count)
        return [{"operation": "bind_file_level_logical_manifest", "target": "config/logical_inputs.tsv",
                 "rows": str(logical_count), "detail": "already projected inventory"}]
    expected = expected_preprojection_problems(inputs)
    if problems != expected or omitted_present != OMIT_PRIVATE_PROVENANCE:
        missing = sorted(expected - problems)
        extra = sorted(problems - expected)
        raise SystemExit(f"unsafe/stale public projection state: missing={missing[:8]} extra={extra[:8]} omissions={sorted(OMIT_PRIVATE_PROVENANCE-omitted_present)}")

    figure_manifest = inputs / FIGURE_MANIFEST
    header, rows = read_tsv(figure_manifest)
    role_index = header.index("role")
    replacements = 0
    for row in rows:
        if "no-DCAN" in row[role_index]:
            row[role_index] = row[role_index].replace("no-DCAN", "source-excluded")
            replacements += 1
    if replacements != 8:
        raise SystemExit(f"figure-manifest excluded-source metadata drift: {replacements} != 8")
    write_and_verify_tsv(figure_manifest, header, rows)

    cds = inputs / CDS_MAPPING
    header, rows = read_tsv(cds)
    if "cds_path" not in header:
        raise SystemExit("CDS private-provenance field shape drift")
    cds_path_index = header.index("cds_path")
    cds_paths = [row[cds_path_index] for row in rows]
    if not any(PRIVATE_PATTERN.search(value) for value in cds_paths) or any(value and not PRIVATE_PATTERN.search(value) for value in cds_paths):
        raise SystemExit("CDS private-provenance values contain an undeclared public/mixed form")
    keep = [index for index, field in enumerate(header) if field != "cds_path"]
    projected_header = [header[index] for index in keep]
    projected_rows = [[row[index] for index in keep] for row in rows]
    write_and_verify_tsv(cds, projected_header, projected_rows)

    for relative in sorted(ROW_HASH_TSVS - OMIT_PRIVATE_PROVENANCE):
        path = inputs / relative
        header, rows = read_tsv(path)
        removed = hash_fields(header)
        if not removed:
            raise SystemExit(f"declared row-hash projection field absent: {relative}")
        keep = [index for index, field in enumerate(header) if field not in removed]
        projected_header = [header[index] for index in keep]
        projected_rows = [[row[index] for index in keep] for row in rows]
        write_and_verify_tsv(path, projected_header, projected_rows)

    for relative in sorted(OMIT_PRIVATE_PROVENANCE):
        path = inputs / relative
        path.unlink()

    for path in gff_paths:
        project_gff_header(path)

    # Member transformations above intentionally change byte sizes.  The
    # public four-field figure manifest is their inventory consumer, so bind it
    # to the final projected bytes only after every member rewrite is complete.
    refresh_figure_inventory(inputs)

    remaining = public_problem_set(inputs)
    if remaining:
        raise SystemExit(f"public projection left forbidden content: {sorted(remaining)[:12]}")
    logical_count = rebuild_input_inventory(inputs)
    validate_projected_inventory(inputs)
    return projection_spec(inputs, logical_count)
from package_archive import write_archive
from public_inputs import project_repository_inputs as project_intron_inputs, project_gff_header, gff_private_header_count, expected_preprojection_problems as intron_expected_preprojection_problems

ROOT = Path(os.environ.get("PCWDE_REPO_ROOT", Path(__file__).resolve().parents[2])).resolve()

if len(sys.argv) == 2 and sys.argv[1] == "--validate-inputs-only":
    selected_root = resolve_data_root(ROOT)
    validate_reader_inputs(ROOT, selected_root)
    print("DEPOSIT_INPUT_FORM_PASS")
    raise SystemExit(0)

if len(sys.argv) == 3 and sys.argv[1] == "--project-inputs-only":
    selected_inputs = Path(sys.argv[2]).resolve()
    prior = selected_inputs.parent / "PUBLIC_PROJECTION.tsv"
    if prior.is_file() and all((selected_inputs / relative).is_file() for relative in S4_SEQUENCE_DIGEST_TSVS):
        project_repository_inputs(selected_inputs, prior)
    else:
        project_intron_inputs(selected_inputs, ROOT)
    if selected_inputs.name == "repository_inputs":
        if not (selected_inputs.parent / "PROJECTED_INPUTS.tsv").exists():
            write_projected_inputs(selected_inputs.parent)
        validate_projected_inputs(selected_inputs.parent)
    print("PUBLIC_PROJECTION_PASS")
    raise SystemExit(0)

DATA_ROOT = resolve_data_root(ROOT)
validate_reader_inputs(ROOT, DATA_ROOT)
SOURCE_INPUTS = DATA_ROOT / "repository_inputs"
LOGICAL_MANIFEST = ROOT / "config/logical_inputs.tsv"
logical_manifest_before = hashlib.sha256(LOGICAL_MANIFEST.read_bytes()).hexdigest()
OUTROOT = Path(os.environ.get("PCWDE_DEPOSIT_OUT", ROOT / "build/data_deposit")).resolve()
OUT = OUTROOT / "pcwde_supporting_data"

if OUTROOT.exists():
    raise SystemExit(f"refusing to overwrite existing deposit output: {OUTROOT}")
OUT.mkdir(parents=True)
shutil.copytree(SOURCE_INPUTS, OUT / "repository_inputs")
INPUTS = OUT / "repository_inputs"
# The chromosome-gain evidence is already sealed as repository logical inputs.
# Keep it in the public input layer so a clean checkout can rebuild the report
# from the extracted archive without resolving upstream workspace paths.
with (ROOT / "config/chromosome_gain_sources.tsv").open(newline="", encoding="utf-8") as handle:
    chromosome_gain_sources = list(csv.DictReader(handle, delimiter="\t"))
if len(chromosome_gain_sources) != 59 or len({r["role"] for r in chromosome_gain_sources}) != 59:
    raise SystemExit("incomplete chromosome-gain input manifest")
for record in chromosome_gain_sources:
    packaged = INPUTS / "chromosome_gain" / record["role"]
    if not packaged.is_file() or packaged.stat().st_size != int(record["deposit_bytes"]) or hashlib.sha256(packaged.read_bytes()).hexdigest() != record["deposit_sha256"]:
        raise SystemExit(f"chromosome-gain input mismatch: {record['role']}")
projection_rows = project_repository_inputs(INPUTS, SOURCE_INPUTS.parent / "PUBLIC_PROJECTION.tsv")
write_projected_inputs(OUT)
validate_projected_inputs(OUT)
validate_reader_inputs(ROOT, OUT)
SRC = INPUTS / "deposit/base"
READER = INPUTS / "deposit/reader"
BUILD = OUTROOT.parent

for name in ("alignments", "data", "sequences", "sidecars", "trees"):
    shutil.copytree(SRC / name, OUT / name,
                    ignore=(lambda directory, names: ["source_panel.tsv"]
                            if Path(directory) == SRC / "data" else []))

# Use the selection membership table produced by this manuscript build.
selection_member_output = BUILD / "table_s7/selection_membership_and_public_source.tsv"
if not selection_member_output.is_file():
    raise SystemExit(f"missing current selection membership table: {selection_member_output}")
shutil.copy2(selection_member_output, OUT / "data/selection/selection_membership_and_public_source.tsv")
sequence_set_dir = OUT / "data/sequence_sets"
sequence_set_dir.mkdir()
for name in ("pcwde_sequence_set_membership.tsv", "pcwde_sequence_set_summary.tsv"):
    source = BUILD / "sequence_sets" / name
    if not source.is_file():
        raise SystemExit(f"missing same-build sequence set output: {source}")
    shutil.copy2(source, sequence_set_dir / name)
rank_dir = OUT / "data/cazy_subfamily"
rank_dir.mkdir()
for built, name in (("cazy_deposit/all_unit_counts.tsv", "all_unit_counts.tsv"),
                    ("cazy_ranks/paired_ranks.tsv", "paired_ranks.tsv"),
                    ("cazy_ranks/primary_ecami_ranks.tsv", "primary_ecami_ranks.tsv"),
                    ("cazy_rank_table/per_genome_cazy_ranks.tsv", "per_genome_cazy_ranks.tsv"),
                    ("cazy_ranks/native_component_counts.tsv", "native_component_counts.tsv"),
                    ("cazy_summary/lineage_summary.tsv", "lineage_summary.tsv"),
                    ("cazy_summary/top_native.tsv", "top_native.tsv")):
    source = BUILD / built
    if not source.is_file():
        raise SystemExit(f"missing same-build CAZy rank output: {source}")
    shutil.copy2(source, rank_dir / name)
gain_dir = OUT / "data/chromosome_gain"
gain_dir.mkdir()
for built, name in (("chromosome_gain/ecami_scores.tsv", "ecami_scores.tsv"),
                    ("chromosome_gain/ecami_tip_states.tsv", "ecami_tip_states.tsv"),
                    ("chromosome_gain/ecami_node_changes.tsv", "ecami_node_changes.tsv"),
                    ("chromosome_gain/background_eligibility.tsv", "background_eligibility.tsv"),
                    ("chromosome_gain/hand_checks.tsv", "hand_checks.tsv"),
                    ("chromosome_gain/report/focal_summary.tsv", "focal_summary.tsv"),
                    ("chromosome_gain/report/copy_strata.tsv", "copy_strata.tsv"),
                    ("chromosome_gain/report/focal_scan_sensitivity.tsv", "focal_scan_sensitivity.tsv"),
                    ("chromosome_gain/focal_component_sizes.tsv", "focal_component_sizes.tsv"),
                    ("chromosome_gain/focal_scan_eligibility.tsv", "focal_scan_eligibility.tsv"),
                    ("chromosome_gain/report/focal_scan_max_diagnostic.tsv", "focal_scan_max_diagnostic.tsv")):
    source = BUILD / built
    if not source.is_file():
        raise SystemExit(f"missing same-build eCAMI chromosome-gain output: {source}")
    if name in {"ecami_scores.tsv", "ecami_node_changes.tsv", "focal_scan_sensitivity.tsv", "focal_scan_max_diagnostic.tsv"}:
        with source.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields = [field for field in reader.fieldnames if not field.startswith(("acctran_", "sankoff_"))]
            records = [row for row in reader if name != "ecami_node_changes.tsv" or row["mode"] == "deltran"]
        if not records:
            raise SystemExit(f"empty DELTRAN projection: {source}")
        with (gain_dir / name).open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader(); writer.writerows({field: row[field] for field in fields} for row in records)
    else:
        shutil.copy2(source, gain_dir / name)

# Catalytic reader outputs come from this build's current origin-domain set.
for name in ("catalytic_reader_table.tsv", "catalytic_residue_summary.tsv", "catalytic_residue_by_sequence.tsv"):
    built = BUILD / "table_s7" / name
    if not built.is_file():
        raise SystemExit(f"missing current catalytic product: {built}")
    shutil.copyfile(built, OUT / "data/selection" / name)

# The reader FASTAs carry the exact sequence/ID pairs used by the accepted
# complete origin trees; the saved subset exports paired some IDs incorrectly.
for family in ("GH5_2", "PL1_8", "PL4_4"):
    shutil.copy2(SRC / f"sequences/origin_context_{family}_full.faa",
                 OUT / f"sequences/origin_context_{family}.faa")

# Reader membership and placement tables come from this build (gene-copy unit), not from saved copies.
# Reader labels: the family is written "B-lectin" in both tables, and the placement column
# holding RagTag/JCVI/NCBI-BK/unresolved is named placement_route. Values are otherwise unchanged.
READER_HEADERS = {"placement_confidence": "placement_route"}


def write_reader_table(built, target):
    lines = built.read_text(encoding="utf-8").splitlines(keepends=True)
    header = [READER_HEADERS.get(c, c) for c in lines[0].rstrip("\r\n").split("\t")]
    family = header.index("family_id")
    out = ["\t".join(header) + "\n"]
    for line in lines[1:]:
        cells = line.rstrip("\r\n").split("\t")
        if cells[family] == "B_lectin":
            cells[family] = "B-lectin"
        out.append("\t".join(cells) + "\n")
    target.write_text("".join(out), encoding="utf-8")


for rel, built in (("data/family_membership.tsv", BUILD / "table_s2/candidate/data/family_membership.tsv"),
                   ("data/chromosome_placement_evidence.tsv", BUILD / "tables_s3_s4/candidate/data/chromosome_placement_evidence.tsv")):
    if not built.is_file():
        raise SystemExit(f"missing build output for the reader deposit: {built}")
    if (OUT / rel).exists():
        raise SystemExit(f"a saved copy of {rel} is still in deposit/base; the build's version must be the only source")
    write_reader_table(built, OUT / rel)
for rel, built in (("data/accepted_copy_counts.tsv", BUILD / "table_s2/candidate/data/accepted_copy_counts.tsv"),
                   ("data/domain_summary.tsv", BUILD / "table_s2/candidate/data/domain_summary.tsv")):
    if not built.is_file() or (OUT / rel).exists():
        raise SystemExit(f"missing or conflicting new deposit source: {built}")
    shutil.copy2(built, OUT / rel)

# Publish the intron observations and summaries produced in this build.
intron_out = OUT / "data/introns"
if intron_out.exists():
    if any(intron_out.iterdir()):
        raise SystemExit("saved intron reader data remain in deposit/base")
else:
    intron_out.mkdir(parents=True)
for source, names in (
    (BUILD / "introns/observations", ("GH5_2_intron_metadata.tsv", "PL1_8_intron_metadata.tsv", "PL4_4_intron_metadata.tsv", "B_lectin_intron_metadata.tsv", "per_intron_observations.tsv")),
    (BUILD / "introns/summary", ("recurrent_position_summary.tsv",)),
    (BUILD / "introns/table_s5", ("Table_S5.tsv",)),
    (BUILD / "introns/gc3", ("pcwde_gc3_per_copy.tsv", "blectin_gc3_per_copy.tsv", "pcwde_species_pairs.tsv", "blectin_gc3_comparison.tsv", "pcwde_gc3_stats.json", "blectin_gc3_stats.json")),
):
    for name in names:
        built = source / name
        if not built.is_file():
            raise SystemExit(f"missing built intron data: {built}")
        shutil.copy2(built, intron_out / name)

# Replace analysis-internal synteny/expression projections with the two
# author-approved reader-facing objects.
shutil.rmtree(OUT / "data/synteny")
(OUT / "data/synteny").mkdir()
subprocess.run([sys.executable, str(ROOT / "scripts/analysis/synteny/build_reader_workbook.py"),
    "--blocks", str(INPUTS / "tables/authority/synteny_corrected/block_summary.tsv"),
    "--states", str(INPUTS / "tables/authority/synteny_corrected/current_region_status_complete.tsv"),
    "--flanks", str(BUILD / "corrected_synteny/flanking_gene_maps_source.tsv"),
    "--roster", str(INPUTS / "tables/authority/synteny_current/final_49_roster.tsv"),
    "--counts", str(BUILD / "tables_s3_s4/synteny_block_counts.tsv"),
    "--out", str(OUT / "data/synteny/conserved_synteny_blocks.xlsx")], check=True)
shutil.copy2(BUILD / "tables_s3_s4/candidate/data/synteny/block_summary.tsv",
             OUT / "data/synteny/block_summary.tsv")
shutil.rmtree(OUT / "data/expression")
(OUT / "data/expression").mkdir()
sys.path.insert(0, str(ROOT / "scripts/tables"))
from expression_projection import project_reader_expression
project_reader_expression(
    BUILD / "table_s6/data/expression/retained_sample_crosswalk.tsv",
    READER / "gene_by_sample_tpm.tsv",
    OUT / "data/expression/gene_by_sample_tpm.tsv",
)

# Public executable inputs were copied and validated before any package
# projection; all subsequent reads use that output-local immutable subtree.

with (OUT / "PUBLIC_PROJECTION.tsv").open("w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=list(projection_rows[0]), delimiter="\t", lineterminator="\n")
    w.writeheader(); w.writerows(projection_rows)

# Public metadata projection consumes the complete table from this build.
from source_panel import build_source_panel
build_source_panel(BUILD / "table_s2/candidate/data/source_panel.tsv", OUT / "data/source_panel.tsv")

selection_qc = OUT / "data/selection/selection_alignment_qc.tsv"
with selection_qc.open(newline="", encoding="utf-8") as fh:
    selection_rows = list(csv.DictReader(fh, delimiter="\t"))
selection_fields = ["target_full_protein_count" if f == "locked_full_protein_targets" else f
                    for f in selection_rows[0] if not f.endswith("_sha256")]
for row in selection_rows:
    row["target_full_protein_count"] = row.pop("locked_full_protein_targets")
with selection_qc.open("w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=selection_fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
    w.writeheader(); w.writerows(selection_rows)

# Row-level sequence hashes were production/QC fields. Package-level checksums
# remain in MANIFEST.tsv and SHA256SUMS.txt.
for sidecar in sorted((OUT / "sidecars").glob("*.tsv")):
    with sidecar.open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        continue
    fields = [f for f in rows[0] if not f.endswith("_sha256")]
    with sidecar.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader(); w.writerows(rows)

membership_file = OUT / "data/family_membership.tsv"
with membership_file.open(newline="", encoding="utf-8") as fh:
    membership_rows = list(csv.DictReader(fh, delimiter="\t"))
membership_count = len(membership_rows)
genome_copy_count = sum(r["unit_type"] == "gene_copy" for r in membership_rows)

descriptions = {
    "data/chromosome_gain/ecami_scores.tsv": "Retained BK-chromosome gain and loss counts for accepted PCWDE families, their largest-component eCAMI sensitivities and eligible native eCAMI groups.",
    "data/chromosome_gain/ecami_tip_states.tsv": "Placed BK-chromosome states and gene-copy counts per family/group and cynipoid genome.",
    "data/chromosome_gain/ecami_node_changes.tsv": "Focal retained gains and losses on the pruned Figure 1 cynipoid species tree.",
    "data/chromosome_gain/background_eligibility.tsv": "Native eCAMI crown-presence and five-placed-tip eligibility under the shared placement route.",
    "data/chromosome_gain/hand_checks.tsv": "Five focal and five native gain events re-scored against explicit tip chromosome states.",
    "data/chromosome_gain/focal_summary.tsv": "DELTRAN focal roots, inferred chromosome-presence gains and losses, rates, copy strata, native percentiles and clade-localized gains.",
    "data/chromosome_gain/copy_strata.tsv": "Native eCAMI copy-stratum group counts, median retained gains and zero-gain shares.",
    "data/chromosome_gain/focal_scan_sensitivity.tsv": "Fixed largest-component eCAMI scan route for each focal PCWDE family across all 48 genomes.",
    "data/chromosome_gain/focal_component_sizes.tsv": "Gene-collapsed focal component sizes and selected fixed eCAMI group IDs.",
    "data/chromosome_gain/focal_scan_eligibility.tsv": "Fixed focal eCAMI component coverage and unchanged crown eligibility.",
    "data/chromosome_gain/focal_scan_max_diagnostic.tsv": "Historical genome-wise maximum scan trajectory retained as diagnostic only.",
    "data/cazy_subfamily/all_unit_counts.tsv": "Gene counts for every genome and CAZy unit in the competitive scan.",
    "data/cazy_subfamily/paired_ranks.tsv": "Accepted-member and scan-derived PCWDE counts and ranks against the same native background.",
    "data/cazy_subfamily/primary_ecami_ranks.tsv": "Largest focal eCAMI group, pooled scan and accepted-member counts ranked against native eCAMI groups.",
    "data/cazy_subfamily/per_genome_cazy_ranks.tsv": "Full per-genome primary, pooled, accepted-member and mixed-unit ranks underlying the condensed lineage reader table.",
    "data/cazy_subfamily/native_component_counts.tsv": "Per-genome native eCAMI group counts used as the primary rank background.",
    "data/cazy_subfamily/lineage_summary.tsv": "Lineage summaries of paired PCWDE subfamily ranks.",
    "data/cazy_subfamily/top_native.tsv": "Most abundant native CAZy units by lineage.",
    "data/source_panel.tsv": "Sequence-resource roster, public identifiers, annotation basis, BUSCO values where applicable, and analysis-panel status.",
    "data/family_membership.tsv": f"Accepted and excluded family members ({membership_count:,} rows). unit_type: gene_copy for the {genome_copy_count:,} genome gene copies (protein models after fragment merging), domain, protein or locus for other units. Families: GH5_2, PL1_8, PL4_4, B-lectin.",
    "data/accepted_copy_counts.tsv": "Accepted gene-copy counts for each of 49 genomes and four target families; the detailed companion to the condensed reader tables.",
    "data/domain_summary.tsv": "Shared Table S3 domain and screen-hit summary from pinned origin and accessory domain sets.",
    "data/chromosome_placement_evidence.tsv": f"Placement of the {genome_copy_count:,} genome gene copies on the Belonocnema kinseyi reference chromosomes (NCBI assembly 31): placement_class direct_chromosome or unplaced; placement_route is RagTag, JCVI, NCBI-BK or unresolved. Families: GH5_2, PL1_8, PL4_4, B-lectin.",
    "data/synteny/conserved_synteny_blocks.xlsx": workbook_description(BUILD / "tables_s3_s4/synteny_summary.tsv"),
    "data/synteny/block_summary.tsv": "Full 201-block synteny list with a readable reference-locus description (species and BK chromosome, or undetermined) and per-block recovery counts.",
    "data/introns/per_intron_observations.tsv": "Individual observed introns mapped to family reference coordinates.",
    "data/introns/recurrent_position_summary.tsv": "Summary of recurrent intron positions and eligible-protein denominators.",
    "data/expression/gene_by_sample_tpm.tsv": "Tidy gene-by-sample TPM values with family, species, sample, study, life-stage, tissue, and available mapping-QC context.",
    "data/selection/selection_reader_table.tsv": "Family-level selection summary corresponding to Supplementary Table S7.",
    "data/selection/selection_family_summary.tsv": "Family-level FUBAR and aBSREL result counts.",
    "data/selection/selection_alignment_qc.tsv": "Protein/CDS recovery and codon-alignment filtering statistics.",
    "data/selection/selection_membership_and_public_source.tsv": "Selection-analysis membership and public source information.",
    "data/selection/fubar_site_evidence.tsv": "Per-site FUBAR parameter estimates and posterior classifications.",
    "data/selection/absrel_branch_evidence.tsv": "Per-branch aBSREL statistics, including corrected significance.",
    "data/selection/catalytic_reader_table.tsv": "Reader-level catalytic-residue conservation summary.",
    "data/selection/catalytic_residue_summary.tsv": "Catalytic-residue conservation by family and reference site.",
    "data/selection/catalytic_residue_by_sequence.tsv": "Sequence-level observations at each evaluated catalytic residue.",
    "data/sequence_sets/pcwde_sequence_set_membership.tsv": "Per-sequence membership in locus-tree, selection, de novo catalytic alignment, and origin-tree sets, with source IDs and link status.",
    "data/sequence_sets/pcwde_sequence_set_summary.tsv": "Per-family set totals, linked and unlinked sequence counts, origin outgroups, and pairwise differences among linked rows.",
}

def description(rel: str) -> str:
    if rel in descriptions:
        return descriptions[rel]
    if rel.startswith("alignments/origin_context_") and "_full.aln.faa" in rel:
        return "Complete transfer-context inference alignment: 427 GH5_2, 364 PL1_8 or 444 PL4_4 sequences after post-inference removal of CB-derived Bassettia pallida rows."
    if rel.startswith("sequences/origin_context_") and rel.endswith("_full.faa"):
        return "Original unaligned MAFFT input for the complete transfer-context inference, with domain trimming already applied and CB-derived Bassettia pallida records removed."
    if rel.startswith("trees/origin_context_") and rel.endswith("_full.treefile"):
        return "Complete transfer-context inference tree, before representative Figure 3 or Figure S1 display pruning."
    if rel.startswith("sidecars/origin_context_") and rel.endswith("_method.txt"):
        return "Original IQ-TREE version, model, seed, bootstrap method and input relationship."
    if rel == "sidecars/removed_cb_origin_bp_tips.tsv":
        return "Bassettia pallida tips removed after inference as Callirhytis quercusbatatoides-derived, with the tree distances and protein identity that met the removal rule."
    if rel == "sidecars/origin_context_full_tip_mapping.tsv":
        return "One-to-one mapping of original alignment identifiers to IQ-TREE tip spelling."
    if "analysis/domain_preparation/" in rel:
        return "Frozen family proteins, admitted-protein hmmsearch hits and historical trimmed FASTAs for preprocessing comparisons."
    if rel.startswith("alignments/origin_context_"):
        return "Retained origin-context subset alignment; current executable Figure S1 inputs are under repository_inputs/figures/figure3."
    if rel.startswith("alignments/complete_locus_"):
        return "Protein alignment used to infer the complete locus tree corresponding to Supplementary Figure S3, including the single bacterial rooting sequence."
    if rel == "alignments/B_lectin_domains.aln.faa":
        return "B-lectin domain-copy alignment corresponding to Supplementary Figure S4."
    if rel.startswith("alignments/catalytic_domains/"):
        return "Profile-guided catalytic-domain alignment corresponding to Supplementary Figure S2."
    if rel.startswith("sequences/catalytic_domains/"):
        return "Unaligned catalytic-domain sequences corresponding to Supplementary Figure S2."
    if rel.startswith("sequences/origin_context_display_domains/"):
        return "Domain-trimmed sequences used in the representative Figure 3 display panels (29/28/29 tips)."
    if rel.startswith("sequences/complete_locus_"):
        return "Unaligned proteins used to infer the complete locus tree for Supplementary Figure S3, including the single bacterial rooting sequence."
    if rel.startswith("sequences/"):
        return "Unaligned protein sequences corresponding to the named supplementary tree."
    if rel.startswith("sidecars/"):
        return "Public identifier, label, or domain-architecture sidecar for the named tree."
    if rel.endswith(".display_rooted_ladderized.treefile"):
        return "Rooted and ladderized catalytic-domain display tree corresponding to Supplementary Figure S2."
    if rel.startswith("trees/complete_locus_") and rel.endswith(".inference_with_bacterial_root.treefile"):
        return "Complete locus inference tree retaining the single bacterial rooting sequence used in the phylogenetic analysis for Supplementary Figure S3."
    if rel.startswith("trees/complete_locus_") and rel.endswith(".display_diplolepididae_rooted_ladderized.treefile"):
        return "Supplementary Figure S3 display tree after removal of the bacterial tip, rerooting on the largest exclusively Diplolepididae branch, and ladderization."
    if rel.endswith(".unrooted.treefile"):
        return "Unrooted maximum-likelihood catalytic-domain inference tree corresponding to Supplementary Figure S2."
    if rel.startswith("trees/"):
        return "Newick tree used for the named supplementary figure."
    return "Supporting data file."

def dictionary_format(path: Path) -> str:
    suffix = path.suffix.lower()
    if suffix == ".domtblout": return "HMMER domain table"
    if suffix == ".tsv": return "TSV"
    if suffix in {".faa", ".fa", ".fasta", ".fna"}: return "FASTA"
    if suffix == ".xlsx": return "Excel workbook"
    if suffix in {".txt", ".md", ".r", ".log"}: return "Plain text"
    if suffix == ".fai": return "FASTA index"
    if suffix == ".pbi" and path.name.endswith(".bam.pbi"): return "PacBio BAM index"
    if suffix in {".treefile", ".nwk", ".tree"}: return "Newick"
    if suffix in {".gff3", ".gff"} or path.name.endswith((".gff.gz", ".gff3.gz")): return "GFF3"
    if suffix == ".bed": return "BED"
    if suffix == ".agp" or path.name.startswith("agp_"): return "AGP"
    if suffix == ".simple" and path.name.endswith(".anchors.simple"): return "JCVI anchor table"
    if suffix == ".svg": return "SVG"
    if suffix == ".json": return "JSON"
    raise SystemExit(f"unknown deposit format: {path}")


def count_records(path: Path, fmt: str) -> int:
    if fmt == "HMMER domain table":
        with path.open(encoding="utf-8") as fh:
            return sum(bool(line.strip()) and not line.startswith("#") for line in fh)
    if fmt == "FASTA index":
        with path.open(encoding="utf-8") as fh:
            return sum(bool(line.strip()) for line in fh)
    if fmt == "PacBio BAM index":
        return 1  # Binary index file; its read count is computed by supplement_numbers.
    if path.suffix == ".tsv":
        with path.open(encoding="utf-8") as fh:
            return max(sum(1 for _ in fh) - 1, 0)
    if path.suffix in {".faa", ".fa", ".fasta", ".fna"}:
        with path.open(encoding="utf-8") as fh:
            return sum(line.startswith(">") for line in fh)
    if path.suffix == ".xlsx":
        from openpyxl import load_workbook
        workbook = load_workbook(path, read_only=True)
        try:
            return len(workbook.sheetnames) - 1
        finally:
            workbook.close()
    if fmt in {"GFF3", "AGP", "BED", "JCVI anchor table"}:
        count = 0
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path, "rt", encoding="utf-8") as fh:
            for line in fh:
                stripped = line.strip()
                if fmt == "GFF3" and stripped == "##FASTA":
                    break
                if not stripped or stripped.startswith("#"):
                    continue
                if fmt == "BED" and (stripped.startswith("track ") or stripped.startswith("browser ")):
                    continue
                if fmt == "JCVI anchor table" and stripped.startswith("###"):
                    continue
                count += 1
        return count
    return 1

dictionary_rows = []
for path in sorted(p for p in OUT.rglob("*") if p.is_file()):
    rel = path.relative_to(OUT).as_posix()
    fmt = dictionary_format(path)
    columns = ""
    if fmt == "TSV":
        with path.open(encoding="utf-8") as fh:
            columns = fh.readline().rstrip("\r\n").replace("\t", "; ")
    desc = description(rel)
    if fmt in {"GFF3", "AGP", "BED", "JCVI anchor table"}:
        unit = {"GFF3": "feature rows", "AGP": "assembly component/gap rows", "BED": "interval rows", "JCVI anchor table": "paired syntenic span rows"}[fmt]
        desc = f"{desc} Records count nonblank {unit}, excluding comments and separators."
    dictionary_rows.append({"file": rel, "format": fmt, "records": count_records(path, fmt), "columns_if_tsv": columns, "description": desc})

membership_dictionary = next(r for r in dictionary_rows if r["file"] == "data/family_membership.tsv")
if membership_dictionary["records"] != membership_count or f"({membership_count:,} rows)" not in membership_dictionary["description"] or f"{genome_copy_count:,} genome gene copies" not in membership_dictionary["description"]:
    raise SystemExit("family membership data dictionary disagrees with shipped table")
with (OUT / "DATA_DICTIONARY.tsv").open("w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=list(dictionary_rows[0]), delimiter="\t", lineterminator="\n")
    w.writeheader(); w.writerows(dictionary_rows)

(OUT / "README.md").write_text("""# Supporting data for *Multiple ancient gene transfers from bacteria facilitated herbivory in cynipoid wasps*

This package contains the machine-readable evidence supporting the comparative genomic, phylogenetic, synteny, intron, expression, and molecular-evolution analyses in the manuscript. Reader-formatted figures and supplementary tables are distributed with the manuscript; this archive provides their underlying sequences, alignments, trees, and tabular evidence.

This local pre-submission package has not been published to the reserved Figshare DOI.

## Contents

- `data/`: source roster, family membership, chromosome placement, retained chromosome-gain comparisons against native eCAMI groups, conserved synteny, introns, expression, selection, and catalytic-residue conservation. Conserved synteny is presented in one workbook with direct `Block N` names; expression is one tidy gene-by-sample table.
- `data/sequence_sets/`: one row per sequence used in the locus tree, frozen selection analysis, de novo catalytic alignment, or transfer-context origin tree. The membership file records the source ID for each set. `link_status=not_linked` denotes a model-version difference whose locus equivalence is not established; it is not an exclusion. Pairwise membership differences in the summary use linked rows; origin-tree outgroups are counted separately.
- The sequence-set membership columns give family, canonical sequence ID, public species name, source type and accession where recorded, BK chromosome where available, four `in_*` booleans, each set's source IDs, and `exclusion_note`. A blank metadata cell or exclusion note means the pinned producer did not supply it. Links use identical IDs or pinned crosswalks first, then unique reciprocal CDS overlap on pinned same-contig, same-strand coordinates if available. Protein identity alone does not link rows. Unlinked IDs retain their own rows. The summary reports each set's total, link tier and unlinked count; origin outgroups and post-inference pruning are explicit.
- `sequences/`: unaligned amino-acid sequences used in supplementary phylogenies.
- `alignments/`: aligned amino-acid sequences used in supplementary phylogenies.
- `trees/`: Newick trees. Transfer-context files ending in `_full.treefile` contain 427/364/444 of the 427/365/446 inference tips (see below). Matching `_full.aln.faa` alignments and `_full.faa` original unaligned inputs are distinct from the smaller subset exports. The retained `origin_context` exports have 186/95/125 tips; the current executable Figure S1 trees have 170/67/82 tips and Figure 3 has 29/28/29, all subsets of the complete trees. Current figure inputs are under `repository_inputs/figures/origin_trees/`. The full tip mapping and method records are in `sidecars/`. Complete-locus files distinguish the inference tree retaining its bacterial rooting sequence from the bacterial-tip-free Diplolepididae-rooted, ladderized Figure S3 display tree. Catalytic-domain files named `display_rooted_ladderized` record display rooting/order; files named `unrooted` preserve maximum-likelihood inference topology.
- `sidecars/`: public label, identifier, domain-copy, and domain-architecture crosswalks.
- `repository_inputs/`: checksum-pinned scientific inputs consumed by the public Git code through `config/logical_inputs.tsv` and `PCWDE_DATA_ROOT`.
- `DATA_DICTIONARY.tsv`: one row per scientific file, with record counts, TSV columns, and descriptions.
- `MANIFEST.tsv` and `SHA256SUMS.txt`: byte sizes and SHA-256 checksums.

The `repository_inputs/analysis/domain_preparation/` files support standalone preprocessing. HMMER hit rows are preserved for proteins in the paired frozen FASTAs; metadata uses portable names. Historical trimmed comparison files remain on CARC and are omitted from the current public inputs. The corrected best-E-value, inclusive-overlap rule changes four proteins and yields 426/363/443 entries instead of the historical 427/364/444. These diagnostic outputs have not replaced the manuscript alignments or trees. Full alignment comparisons remove gaps and stop symbols; the tip mapping records IQ-TREE punctuation changes.

Four *Bassettia pallida* transcript-derived sequences were removed after inference because they are near-identical to *Callirhytis quercusbatatoides* sequences, consistent with *C. quercusbatatoides*-like sequence detected in the *B. pallida* RNA-seq libraries: in every tree containing both species, each lies within 0.03 substitutions/site of a *C. quercusbatatoides* tip (at most 0.25 times its distance to any other species) with at least 98% aligned protein identity. They are listed with these values in `sidecars/removed_cb_origin_bp_tips.tsv` and were removed from every deposited tree, from the sequences, alignments (rows only; no realignment or column changes), sidecars and figure inputs that feed those trees, and from the domain-preparation inputs. Deposited trees are the accepted inferences pruned of these tips; all other topology, branch lengths and support values are unchanged. Rerunning inference from the deposited alignments will therefore not reproduce the deposited trees byte-for-byte. Catalytic-residue tables use the current origin-domain membership and exclude three of these recorded tips. Expression, selection and family-membership tables are unchanged.

Blank BUSCO fields in `data/source_panel.tsv` are expected for individual outgroup protein hits and other resources for which assembly/proteome BUSCO was not evaluated. B-lectin records explicitly identify whether the counting unit is a protein, locus, or individual domain copy. Expression values are TPM rounded to the nearest 0.1; mapping QC is included in the same gene-by-sample table where available. Selection classifications use FUBAR posterior probability > 0.9, and aBSREL significance uses Holm-corrected P < 0.05. Catalytic-residue percentages use observable residues as the denominator.

""", encoding="utf-8")

(OUT / "FIGSHARE_METADATA.md").write_text("""# Figshare metadata

Title: Supporting data for “Multiple ancient gene transfers from bacteria facilitated herbivory in cynipoid wasps”

Reserved DOI (not yet published): 10.6084/m9.figshare.33265659

Version: 1.0.0

Archive: `pcwde_supporting_data.tar.gz`

Archive checksum: recorded in the adjacent `ARCHIVE.sha256` sidecar generated with this archive

Authors: Adam J. Kranz; Y. Miles Zhang; Maryam R. Faghihi; Sheina B. Sim; Renee L. Corpuz; Scott M. Geib; Guerin E. Brown; Andrew A. Forbes; Jack Hearn; Scott P. Egan; Ellen O. Martinson

Description: Machine-readable source rosters, family membership, chromosome placement, synteny, intron, expression, selection, catalytic-residue, sequence, alignment, and tree data supporting the associated manuscript.

Keywords: Cynipoidea; gall wasps; lateral gene transfer; plant cell wall-degrading enzymes; comparative genomics; synteny; phylogenetics; gene expression

License: CC BY 4.0


Extraction: verify `ARCHIVE.sha256`, then run `tar -xzf pcwde_supporting_data.tar.gz`. The archive creates exactly one `pcwde_supporting_data/` root. Set `PCWDE_DATA_ROOT` to that extracted directory.
""", encoding="utf-8")

scientific_files = [p for p in OUT.rglob("*") if p.is_file() and p.name not in {"MANIFEST.tsv", "SHA256SUMS.txt"}]
def sha256_file(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()
manifest_rows = []
for path in sorted(scientific_files):
    rel = path.relative_to(OUT).as_posix()
    digest = sha256_file(path)
    manifest_rows.append((rel, path.stat().st_size, digest))
with (OUT / "MANIFEST.tsv").open("w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["relative_path", "bytes", "sha256"]); w.writerows(manifest_rows)
with (OUT / "SHA256SUMS.txt").open("w", encoding="utf-8") as fh:
    for rel, _, digest in manifest_rows:
        fh.write(f"{digest}  {rel}\n")

# Seal the byte-identical public package through the shared archive boundary.
write_archive(OUTROOT)

validate_reader_inputs(ROOT, DATA_ROOT)
if hashlib.sha256(LOGICAL_MANIFEST.read_bytes()).hexdigest() != logical_manifest_before:
    raise SystemExit("tracked logical-input manifest mutated during deposit build")

import check_closure  # noqa: E402  (same directory)
check_closure.write_receipt(OUTROOT.parent / "closure", "scripts/release/build_data_deposit.py", OUTROOT)

print(f"scientific_files={len(dictionary_rows)}")
print(f"package_files={sum(1 for p in OUT.rglob('*') if p.is_file())}")
