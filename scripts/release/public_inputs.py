#!/usr/bin/env python3
"""Import-safe public-input projection shared by the producer and fixtures."""
import csv
import hashlib
import os
import re
from pathlib import Path
from data_root import SOURCE_MANIFEST_SHA256, validate_logical_inputs

# These files contain private provenance omitted from the supporting archive.
# Scientific rows and non-hash scientific cells retain their original order.
OMIT_PRIVATE_PROVENANCE = {
    "analysis/synteny/BUILD_RECEIPT.tsv",
    "analysis/synteny/INPUT_MANIFEST.tsv",
    "analysis/synteny/RESULT_MANIFEST.tsv",
    "analysis/synteny/results/SHA256SUMS.tsv",
    "tables/authority/synteny_current/producer_build_receipt.tsv",
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
S4_SEQUENCE_DIGEST_TSVS = {"analysis/s4_architecture/full_sequence_mapping.tsv","figures/accessory_domain_trees/architecture.tsv"}
FIGURE_MANIFEST = "figures/INPUT_MANIFEST.tsv"
CDS_MAPPING = "tables/authority/selection_current/current_cds_mapping.tsv"
PRIVATE_PATTERN = re.compile(r"/(?:carc|easley|users|root)/")
PRIVATE_BYTES = re.compile(rb"/(?:carc|easley|users|root)/")
INTRON_GFF = re.compile(r"analysis/introns/sources/[^/]+/gff\.gff$")
REMOVED_PATTERN = re.compile(r"dryocosmus[ _]cannoni|\bDCAN\b", re.IGNORECASE)
REMOVED_TOKEN_SUFFIXES = {".tsv", ".txt", ".md", ".faa", ".fa", ".fasta", ".treefile", ".nwk", ".tree", ".csv", ".r", ".py", ".bash", ".sbatch"}
TEXT_SUFFIXES = REMOVED_TOKEN_SUFFIXES | {".fna", ".gff", ".gff3"}
HYBRID_SOURCE_MANIFEST_SHA256 = SOURCE_MANIFEST_SHA256


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
                    if path.suffix.casefold() in REMOVED_TOKEN_SUFFIXES and not relative.startswith("analysis/synteny/") and REMOVED_PATTERN.search(line):
                        problems.add(("removed_taxon_token", relative))
    return problems


def expected_preprojection_problems(inputs: Path, repo: Path) -> set[tuple[str, str]]:
    expected = {("removed_taxon_token", FIGURE_MANIFEST),
                ("private_absolute_path", CDS_MAPPING),
                ("removed_taxon_token", "tables/authority/synteny_current/producer_build_receipt.tsv"),
                }
    expected.update(("row_hash_header", path) for path in ROW_HASH_TSVS)
    # Only declared annotation sources may carry private GFF header provenance.
    source_ids = {row[0] for row in read_tsv(repo / "config/intron_sources.tsv")[1]}
    expected.update(("private_absolute_path", f"analysis/introns/sources/{source}/gff.gff")
                    for source in source_ids if (inputs / f"analysis/introns/sources/{source}/gff.gff").is_file()
                    and gff_private_header_count(inputs / f"analysis/introns/sources/{source}/gff.gff"))
    return expected


def gff_private_header_count(path: Path) -> int:
    count = 0
    with path.open("rb") as handle:
        for line in handle:
            if PRIVATE_BYTES.search(line):
                if not line.startswith(b"#"):
                    raise SystemExit(f"private path in GFF feature line: {path}")
                count += 1
    return count


def project_gff_header(path: Path) -> int:
    count = gff_private_header_count(path)
    if not count:
        return 0
    temporary = path.with_name(f".{path.name}.public_projection.tmp")
    if temporary.exists() or temporary.is_symlink():
        raise SystemExit(f"stale public-projection temporary file: {temporary}")
    mode = path.stat().st_mode & 0o777
    try:
        with path.open("rb") as source, temporary.open("xb") as output:
            for line in source:
                if PRIVATE_BYTES.search(line):
                    output.write(b"# private source path omitted\n")
                else:
                    output.write(line)
        os.chmod(temporary, mode)
        os.replace(temporary, path)
    finally:
        if temporary.exists() and not temporary.is_symlink():
            temporary.unlink()
    return count


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
    for path in sorted((inputs / "analysis/introns/sources").glob("*/gff.gff")):
        relative = path.relative_to(inputs).as_posix()
        digest_state = hashlib.sha256()
        count = 0
        with path.open("rb") as handle:
            for line in handle:
                digest_state.update(line)
                count += line == b"# private source path omitted\n"
        digest = digest_state.hexdigest()
        rows.append({"operation": "neutralize_private_gff_header", "target": relative, "rows": str(count), "detail": f"feature lines unchanged; projected_sha256={digest}"})
    rows.append({"operation": "bind_file_level_logical_manifest", "target": "config/logical_inputs.tsv", "rows": str(logical_count), "detail": "assembler emits one immutable file-level digest and byte count per public logical input"})
    return rows


def project_repository_inputs(inputs: Path, repo: Path, source_manifest_digest: str = HYBRID_SOURCE_MANIFEST_SHA256) -> list[dict[str, str]]:
    problems = public_problem_set(inputs)
    gffs = sorted((inputs / "analysis/introns/sources").glob("*/gff.gff"))
    projected_gffs = []
    for path in gffs:
        with path.open("rb") as handle:
            if any(line == b"# private source path omitted\n" for line in handle):
                projected_gffs.append(path)
    if projected_gffs and problems:
        raise SystemExit("partially projected intron GFF sources")
    omitted_present = {relative for relative in OMIT_PRIVATE_PROVENANCE if (inputs / relative).is_file()}
    if not problems and not omitted_present:
        logical_count = validate_projected_inventory(inputs)
        return projection_spec(inputs, logical_count)
    # The known staged handoff has already undergone the TSV and omission
    # transforms, but its pinned GFF source headers still need neutralization.
    # Admit this state only when every byte matches the committed source pins.
    if problems and not omitted_present and all(kind == "private_absolute_path" and INTRON_GFF.fullmatch(relative)
                                                for kind, relative in problems):
        return migrate_pinned_hybrid(inputs, problems, repo, source_manifest_digest)
    expected = expected_preprojection_problems(inputs, repo)
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

    for path in sorted((inputs / "analysis/introns/sources").glob("*/gff.gff")):
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


def migrate_pinned_hybrid(inputs: Path, problems: set[tuple[str, str]], repo: Path, source_manifest_digest: str = HYBRID_SOURCE_MANIFEST_SHA256) -> list[dict[str, str]]:
    manifest = repo / "config/logical_inputs.tsv"
    if hashlib.sha256(manifest.read_bytes()).hexdigest() != source_manifest_digest:
        raise SystemExit("hybrid source manifest digest mismatch")
    rows = validate_logical_inputs(repo, inputs.parent)
    if len(rows) != len(read_tsv(manifest)[1]):
        raise SystemExit("hybrid source manifest count mismatch")
    expected = expected_preprojection_problems(inputs, repo) & {
        ("private_absolute_path", f"analysis/introns/sources/{row[0]}/gff.gff")
        for row in read_tsv(repo / "config/intron_sources.tsv")[1]
    }
    if problems != expected or not expected:
        raise SystemExit(f"hybrid GFF hazard mismatch: {sorted(problems ^ expected)[:8]}")
    if any((inputs / path).exists() for path in OMIT_PRIVATE_PROVENANCE | {"tables/figure5_current/source_manifest.tsv"}):
        raise SystemExit("hybrid omitted provenance state mismatch")
    if read_tsv(inputs / CDS_MAPPING)[0] != ["family", "protein_id", "species", "source", "source_file", "species_code", "aa_length", "mapping_status", "mapping_reason", "n_distinct_cds"]:
        raise SystemExit("hybrid CDS header mismatch")
    if read_tsv(inputs / FIGURE_MANIFEST)[0] != ["component", "relative_path", "role", "size_bytes"]:
        raise SystemExit("hybrid figure header mismatch")
    for relative in ROW_HASH_TSVS - OMIT_PRIVATE_PROVENANCE:
        if hash_fields(read_tsv(inputs / relative)[0], relative):
            raise SystemExit(f"hybrid row-hash state mismatch: {relative}")
    figure_rows = read_tsv(inputs / FIGURE_MANIFEST)[1]
    # validate_logical_inputs bound this entire figure manifest to its pinned
    # source row, including every role.  Hybrid inputs have already had the
    # no-DCAN projection applied; do not require a role absent from the pin.
    if any("no-DCAN" in row[2] for row in figure_rows):
        raise SystemExit("hybrid figure role state mismatch")
    expected_inventory_count = sum(
        row["figshare_relative_path"] != "repository_inputs/INPUT_MANIFEST.tsv"
        for row in rows
    )
    if len(read_tsv(inputs / "INPUT_MANIFEST.tsv")[1]) != expected_inventory_count:
        raise SystemExit("hybrid inventory cardinality mismatch")
    for _, relative in sorted(problems):
        project_gff_header(inputs / relative)
    refresh_figure_inventory(inputs)
    if public_problem_set(inputs):
        raise SystemExit("hybrid projection left public hazards")
    count = rebuild_input_inventory(inputs)
    validate_projected_inventory(inputs)
    if count != len(rows):
        raise SystemExit("hybrid projected count mismatch")
    return projection_spec(inputs, count)
