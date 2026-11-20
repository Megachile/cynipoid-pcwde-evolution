#!/usr/bin/env python3
"""Fail-closed verification of the public repository and pinned inputs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from data_root import resolve_data_root, validate_reader_inputs


TEXT_SUFFIXES = {".py", ".r", ".R", ".tsv", ".csv", ".txt", ".md", ".faa", ".fa", ".fasta", ".nwk", ".tree", ".treefile", ".json", ".bash", ".sbatch", ".cff", ".tex", ""}
CURRENT_TABLE_S5 = "76bad2a9d744bdc049741b1ace495abc6f371bedb614e7f7ad743f1211af4071"
STALE_TABLE_S5 = "c13c01157e0275c767db415db6837178ac972783cac44fc758b13a8a466ec56d"
TINY_GIT_FILE_LIMIT = 262_144
TOTAL_GIT_BYTE_LIMIT = 8_388_608
MANUSCRIPT_RELATIVE_PATH = "authority/manuscript/manuscript_current.md"
MANUSCRIPT_BYTES = 129011
MANUSCRIPT_SHA256 = "5e6b453cbbde5f5cf4dc98d16057aff9286cdf7b619bb96e9fda83aa39394355"
DRAFT_RELATIVE_PATH = "authority/manuscript/manuscript_draft_20261109T042633Z.md"
DRAFT_BYTES = 114616
DRAFT_SHA256 = "93c0dc55d8fa69e75084c6942de261e426f6a3adfac5b3f83fba107e7acdaaac"
STALE_MANUSCRIPT_SHA256 = "22f987bcdd57e23ecea951f92773a7c92bd0497c52b5d0e9972e34d067c82df9"
SOURCE_AUTHORITY_RELATIVE_PATH = "authority/manuscript/SOURCE_AUTHORITY.md"
SOURCE_AUTHORITY_BYTES = 3370
SOURCE_AUTHORITY_SHA256 = "57392131905c5d2573247ac2d43eb25ce45c44201edca5ec5a255c9eba0eb31c"
FIG8_CAPTION_RELATIVE_PATH = "authority/manuscript/figure8_caption.md"
FIG8_CAPTION_BYTES = 793
FIG8_CAPTION_SHA256 = "8497c24020123c7313af69e77258398922049c17d251afe84e8479e80830b506"
RENDER_ALL_SHA256 = "bfc431c4d96d4d04fd90cdcadbb7d981532decf5d5f8851a42d90c31357f6c44"
CURRENT_DRAFT_IDENTITIES = {'authority/manuscript/manuscript_draft_20261118T0422Z.md': [114713, '138e0fb97788587d3fafb77e8df4b6a314c2a38c39325310603b49a8777efc25'], 'authority/manuscript/manuscript_draft_20261115T1730Z.md': [114601, 'c791e1cb6dbcc4bc2b4d21364ff7e0436f768215ba286cc6d3704e17b8650502'], 'authority/manuscript/manuscript_draft_20261116T0230Z.md': [114616, 'bf2fbf3d411df38edbf09d7cf4136395e6a2e67c64c4336dab7bc5f8bd4e9f2c'], 'authority/manuscript/manuscript_draft_20261116T0330Z.md': [114579, '2d9e8e5ea8cd30ac1fef828bdaa669d3aa93ab8834e0d930620162e096a51a7b'], 'authority/manuscript/manuscript_draft_20261116T0420Z.md': [114592, 'b6975af9d695932b653d9ef5156420d91552badd62dfc0c7b197a88c0f77b723'], 'authority/manuscript/manuscript_draft_20261116T1531Z.md': [114630, '0b1a5886f00215b323db8c24d1cc34ea39f6ba5b50eeaf08c21bc58428ea9d7b'], 'authority/manuscript/manuscript_draft_20261117T1715Z.md': [114628, '12942768a47d3caad15dfd338442148d7eb7495748a8bae379c48b6bfadd432c'], 'authority/manuscript/manuscript_source_pending.md': [114713, '138e0fb97788587d3fafb77e8df4b6a314c2a38c39325310603b49a8777efc25']}
PUBLIC_DOCUMENTATION = {
    "authority/supplement/supplement_draft_20261118T0358Z.md",
    "docs/PUBLIC_SCOPE_PROJECTION.md", "docs/CAZY_PRECISION.md",
    "authority/supplement/supplement_draft_20261117T1715Z.md",
    "authority/supplement/supplement_draft_20261116T1531Z.md",
    "authority/supplement/supplement_draft_20261116T0420Z.md",
    *CURRENT_DRAFT_IDENTITIES, "authority/supplement/supplement_draft_20261115T1730Z.md",
    "scripts/analysis/synteny/modelless_recovery/METHOD.md",    "README.md", "LICENSE", "CITATION.cff", "docs/DATA_ACCESS.md",
    "docs/SUPPLEMENT_COPY_RANK_METHOD.md", "docs/SUPPLEMENT_GAIN_LINEAGE_METHOD.md",
    "docs/SUPPLEMENT_SYNTENY_SUMMARY_METHOD.md", "docs/SUPPLEMENT_TABLE_DOMAIN_SUMMARY_METHOD.md",
    "docs/READER_DAG.md", "docs/ORIGIN_DISPLAY_TREES.md", "docs/SUPPLEMENT_DRAFT_REFRESH.md", "docs/SUPPLEMENT_DOC_EDIT_REPORT_20261110.md", "docs/CATALYTIC_SET_METHOD.md", "docs/CATALYTIC_CLUSTERING_METHOD.md", "docs/figure_s4_architecture_method.md", "docs/figure_s6_rendering.md", "docs/pinned_hybrid_transition.md", "scripts/release/PACKAGING_BOUNDARY.md", "scripts/figures/README.md", "scripts/deposit_inputs/README.md", "scripts/tables/README.md", MANUSCRIPT_RELATIVE_PATH, SOURCE_AUTHORITY_RELATIVE_PATH, FIG8_CAPTION_RELATIVE_PATH, "authority/supplement/revision_draft_text.md", "authority/supplement/supplement_draft_20261102T234911Z.md", "upstream/supplement_table_sources/README.md", "docs/upstream_dbcan_sub.md", DRAFT_RELATIVE_PATH, "authority/supplement/supplement_draft_20261108T023458Z.md", "authority/supplement/supplement_draft_20261109T054239Z.md", "authority/supplement/supplement_draft_20261110T144642Z.md", "authority/supplement/supplement_source.md", "authority/captions/author_captions.json",
}
PUBLIC_MANIFESTS = {
    "config/intron_manifest_bridge.json",
    "config/source_aliases.tsv",
    "config/author_result_routes.json",
    "authority/source_refresh_pending.json", "config/manuscript_current_pointers.json",
    "scripts/analysis/synteny/modelless_recovery/RECOVERY_PINS.tsv",    "config/logical_inputs.tsv", "config/build_inputs.tsv", "config/deposit_sources.tsv", "config/figure_numbers.tsv", "config/closure_declared.tsv", "config/deposit_layout.tsv", "config/upstream_locked_inputs.tsv", "config/chromosome_gain_sources.tsv", "config/deposit_size_baseline.tsv", "config/manuscript_template.tex", "config/supplement_assembly.tex", "config/blectin_models.tsv", "config/intron_genome_sources.tsv", "config/intron_seqid_aliases.tsv", "config/intron_snapshot_sources.tsv", "config/intron_sources.tsv", "config/intron_frozen_sources.tsv", "config/table_numbers.tsv", "config/supplement_sections.tsv", "config/supplement_paragraphs.tsv", "config/supplement_seed_checks.tsv", "authority/supplement/SOURCE_MANIFEST.tsv", "config/supplement_literals.tsv",
    "config/cazy_stage_sources.tsv", "config/cazy_lineage_groups.tsv", "config/cazy_processed_proteomes.tsv", "config/cazy_paper_proteomes.tsv", "config/cazy_scans.tsv", "config/ecami_focal_parity.tsv",
    "config/synteny_derived_pins.tsv", "config/synteny_existing_sources.tsv",
    "config/supplement_accessory_domain_trees.tex", "config/supplement_b_lectin_domain_tree.tex", "config/supplement_complete_locus_trees.tex",
    "config/supplement_euceroptres_genome.tex", "config/supplement_origin_trees_full.tex", "config/supplement_tree_counts_full.tex",
    "docs/methods_code_crosswalk.tsv",
    "environment/python-lock.txt", "environment/r-lock.txt",
}
PUBLIC_CONFIGURATION = {
    "config/public_retired_inputs.txt", "config/public_scope_pins.tsv",".gitignore", "Makefile", "config/excluded_expression_samples.tsv"}
PUBLIC_INTERMEDIATES = {"scripts/analysis/catalytic/clustering/LOCKED_ANALYSIS.json"}
PUBLIC_FIXTURES = {"config/bp_cb_prune_plan.tsv", "config/ecology_states.tsv", "config/tree_counts_remnant.tsv", "config/tree_counts_outgroup_ecology.tsv", "config/origin_display_tips.tsv"}
PREVIEW_DRAFT = "config/fig2_preview_draft.md"
PREVIEW_DRAFT_SHA256 = "5b015b83fa9fc7762a7f87d48936d02f9c6e70ec63f8831982ad70880342bf94"
SYNTENY_PROVENANCE = "scripts/analysis/synteny/PROVENANCE.tsv"
REQUIRED_CLAIMS = {
    "README.md": (
        'Multiple ancient gene transfers from bacteria facilitated herbivory in cynipoid wasps',
        "Some large or externally sourced upstream inputs are not redistributed",
        "Expensive MAFFT, IQ-TREE, and HyPhy results are not rerun by the manuscript target",
        "synteny release code consumes the deposited final-49 locus and block tables",
    ),
    "docs/READER_DAG.md": (
        "does not rerun MAFFT/IQ-TREE inference",
        "or read quantification",
    ),
}


def literal_claim_present(text: str, claim: str) -> bool:
    """Match an exact semantic phrase after collapsing ordinary whitespace."""
    normalized_text = re.sub(r"\s+", " ", text).strip()
    normalized_claim = re.sub(r"\s+", " ", claim).strip()
    return normalized_claim in normalized_text


def admitted_git_class(root: Path, path: Path) -> str | None:
    relative = path.relative_to(root).as_posix()
    matches = []
    if relative.startswith("scripts/") and path.suffix in {".py", ".r", ".R", ".bash", ".sbatch"}:
        matches.append("code")
    if relative.startswith("upstream/introns/") and path.suffix == ".py":
        matches.append("code")
    if relative in PUBLIC_DOCUMENTATION:
        matches.append("documentation")
    if relative in PUBLIC_MANIFESTS or relative == SYNTENY_PROVENANCE:
        matches.append("manifest_or_environment")
    if relative in PUBLIC_CONFIGURATION:
        matches.append("configuration")
    if relative in PUBLIC_FIXTURES:
        matches.append("behavioral_fixture")
    if relative in PUBLIC_INTERMEDIATES:
        matches.append("documented_intermediate")
    if relative == PREVIEW_DRAFT and digest(path) == PREVIEW_DRAFT_SHA256:
        matches.append("presentation_input")
    return matches[0] if len(matches) == 1 else None


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def manuscript_authority_failures(root: Path, files: list[Path]) -> list[str]:
    """Require every exact manifest-bound historical and current source file."""
    failures: list[str] = []
    authority_root = root / "authority/manuscript"
    if authority_root.is_symlink() or not authority_root.is_dir():
        return ["manuscript authority directory is missing or nonregular"]
    objects = sorted(path.relative_to(root).as_posix() for path in authority_root.rglob("*"))
    expected = sorted((MANUSCRIPT_RELATIVE_PATH, SOURCE_AUTHORITY_RELATIVE_PATH, FIG8_CAPTION_RELATIVE_PATH, DRAFT_RELATIVE_PATH, *CURRENT_DRAFT_IDENTITIES))
    if objects != expected:
        failures.append(f"manuscript authority set mismatch expected={expected} actual={objects}")
    identities = {
        MANUSCRIPT_RELATIVE_PATH: (MANUSCRIPT_BYTES, MANUSCRIPT_SHA256),
        SOURCE_AUTHORITY_RELATIVE_PATH: (SOURCE_AUTHORITY_BYTES, SOURCE_AUTHORITY_SHA256),
        FIG8_CAPTION_RELATIVE_PATH: (FIG8_CAPTION_BYTES, FIG8_CAPTION_SHA256),
        DRAFT_RELATIVE_PATH: (DRAFT_BYTES, DRAFT_SHA256),
        **CURRENT_DRAFT_IDENTITIES,
    }
    for relative, (expected_bytes, expected_sha256) in identities.items():
        path = root / relative
        if (
            path.is_symlink()
            or not path.is_file()
            or path.stat().st_size != expected_bytes
            or digest(path) != expected_sha256
        ):
            failures.append(f"manuscript authority identity mismatch: {relative}")
    for path in files:
        relative = path.relative_to(root).as_posix()
        file_digest = digest(path)
        if file_digest == MANUSCRIPT_SHA256 and relative != MANUSCRIPT_RELATIVE_PATH:
            failures.append(f"checksum-pinned manuscript was relabelled: {relative}")
        if file_digest == STALE_MANUSCRIPT_SHA256:
            failures.append(f"forbidden stale manuscript present: {relative}")
        if file_digest == FIG8_CAPTION_SHA256 and relative != FIG8_CAPTION_RELATIVE_PATH:
            failures.append(f"checksum-pinned Figure 8 caption was relabelled: {relative}")
        if file_digest == SOURCE_AUTHORITY_SHA256 and relative != SOURCE_AUTHORITY_RELATIVE_PATH:
            failures.append(f"checksum-pinned source authority was relabelled: {relative}")
    return failures


def lightweight_boundary_failures(root: Path, files: list[Path]) -> list[str]:
    """Apply the universal per-file and total reader-facing Git boundaries."""
    failures = manuscript_authority_failures(root, files)
    oversized = [
        path.relative_to(root).as_posix()
        for path in files
        if path.stat().st_size > TINY_GIT_FILE_LIMIT
    ]
    if oversized:
        failures.append(
            f"Git file exceeds {TINY_GIT_FILE_LIMIT}-byte lightweight boundary: {oversized[:5]}"
        )
    total_bytes = sum(path.stat().st_size for path in files)
    if total_bytes > TOTAL_GIT_BYTE_LIMIT:
        failures.append(
            f"Git candidate exceeds {TOTAL_GIT_BYTE_LIMIT}-byte total boundary: {total_bytes}"
        )
    return failures


def verify_input_binding(root: Path, data_root: Path) -> list[dict[str, str]]:
    """The raw/projected input gate used by the combined reader verifier."""
    return validate_reader_inputs(root, data_root)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--data-root", type=Path)
    parser.add_argument("--repository-only", action="store_true")
    args = parser.parse_args()
    root = args.repo.resolve()
    from portable_sources import validate as validate_portable_sources
    validate_portable_sources(root)
    data_root = None if args.repository_only else resolve_data_root(root, args.data_root)
    logical_path = root / "config/logical_inputs.tsv"
    logical_digest_before = digest(logical_path)
    logical_manifest = [] if args.repository_only else verify_input_binding(root, data_root)
    inputs = None if args.repository_only else data_root / "repository_inputs"
    failures: list[str] = []
    database_objects = [path.relative_to(inputs).as_posix() for path in (inputs.rglob("*") if inputs else [])
                        if path.is_file() and path.suffix.casefold() in {".db", ".sqlite", ".sqlite3"}]
    if database_objects:
        failures.append(f"database snapshots are not public inputs: {database_objects[:5]}")
    forbidden_git_roots = [name for name in ("data", "outputs", "figshare_data") if (root / name).exists()]
    if forbidden_git_roots:
        failures.append(f"scientific/generated roots present in Git: {forbidden_git_roots}")
    files = [path for path in root.rglob("*") if path.is_file()
             and not ({".git", "build", "figshare_data"} & set(path.relative_to(root).parts))]
    symlinks = [path for path in root.rglob("*") if path.is_symlink()
                and not ({".git", "build", "figshare_data"} & set(path.relative_to(root).parts))]
    if symlinks:
        failures.append(f"symlinks present: {symlinks[:5]}")
    forbidden_names = [
        path for path in files
        if (
            (any(part.startswith(".") for part in path.relative_to(root).parts) and path.name != ".gitignore")
            or any("cache" in part.casefold() for part in path.relative_to(root).parts)
            or path.name.endswith(("~", ".bak", ".orig", ".rej", ".pyc", ".pyo", ".tmp", ".temp", ".swp", ".swo"))
            or path.name in {".DS_Store", "Thumbs.db"}
        )
    ]
    if forbidden_names:
        failures.append(f"hidden/cache/backup files: {forbidden_names[:5]}")
    failures.extend(lightweight_boundary_failures(root, files))
    unclassified = [path.relative_to(root).as_posix() for path in files if admitted_git_class(root, path) is None]
    if unclassified:
        failures.append(f"unclassified Git objects: {unclassified[:8]}")

    renderer = root / "scripts/release/render_all.bash"
    if renderer.is_symlink() or not renderer.is_file() or digest(renderer) != RENDER_ALL_SHA256:
        failures.append("accepted render_all.bash entry point missing or mutated")

    if not args.repository_only:
        intron_inputs = inputs / "analysis/introns"
        if (intron_inputs / "sources").exists():
            failures.append("whole intron annotation sources entered repository inputs")
        frozen_inputs = intron_inputs / "frozen"
        if not frozen_inputs.is_dir() or len([p for p in frozen_inputs.rglob("*") if p.is_file()]) != 94:
            failures.append("frozen intron layer-one inventory mismatch")
        manifest_path = inputs / "INPUT_MANIFEST.tsv"
        if not manifest_path.is_file():
            failures.append("missing data/INPUT_MANIFEST.tsv")
        else:
            manifest = read_tsv(manifest_path)
            if list(manifest[0]) != ["relative_path", "size_bytes", "source_class"] if manifest else True:
                failures.append("projected repository inventory schema mismatch")
            expected = {row["relative_path"] for row in manifest}
            actual = {path.relative_to(inputs).as_posix() for path in inputs.rglob("*") if path.is_file() and path != manifest_path}
            actual.discard("figures/figure_s5/Figure_S5.pdf")
            if expected != actual:
                failures.append(f"data manifest inventory mismatch missing={sorted(actual-expected)[:5]} extra={sorted(expected-actual)[:5]}")
            for row in manifest:
                path = inputs / row["relative_path"]
                if not path.is_file() or path.stat().st_size != int(row["size_bytes"]):
                    failures.append(f"data manifest mismatch: {row['relative_path']}")

        figure_manifest = read_tsv(inputs / "figures/INPUT_MANIFEST.tsv")
        if figure_manifest and list(figure_manifest[0]) != ["component", "relative_path", "role", "size_bytes"]:
            failures.append("projected figure inventory schema mismatch")
        # Public scope retires five unused figure inputs from the reviewed 58.
        expected_figure_count = 71 + 1 - 14 - 5
        if len(figure_manifest) != expected_figure_count:
            failures.append(f"figure manifest cardinality: {len(figure_manifest)} != {expected_figure_count}")
        for row in figure_manifest:
            path = inputs / "figures" / row["relative_path"]
            if not path.is_file() or path.stat().st_size != int(row["size_bytes"]):
                failures.append(f"figure manifest mismatch: {row['relative_path']}")
        figure_actual = {
            path.relative_to(inputs / "figures").as_posix()
            for path in (inputs / "figures").rglob("*")
            if path.is_file() and path.name != "INPUT_MANIFEST.tsv"
        }
        figure_actual.discard("figure_s5/Figure_S5.pdf")
        if {row["relative_path"] for row in figure_manifest} != figure_actual:
            failures.append("figure inventory closure mismatch")
        # Analytical inputs and deposit seeds must be free of source-12 identities,
        # retired topology-test products, private paths, and Git LFS placeholders.
        public_roots = [inputs]
        checks = [
            (re.compile(rb"dryocosmus[ _]cannoni|\bDCAN(?:_|\b)", re.I), "source-12 identity"),
            (re.compile(rb"\bAU[- _]?(?:test|topology)\b", re.I), "retired AU route"),
            (re.compile(rb"/(?:carc|easley|users|root)/"), "private absolute path"),
            (re.compile(rb"version https://git-lfs.github.com/spec/v1"), "Git LFS pointer"),
        ]
        for base in public_roots:
            for path in base.rglob("*"):
                if not path.is_file() or path.suffix not in TEXT_SUFFIXES or path.name.endswith("MANIFEST.tsv"):
                    continue
                seen = set()
                with path.open('rb') as handle:
                    for line in handle:
                        for pattern, label in checks:
                            if label not in seen and pattern.search(line):
                                failures.append(f"{label}: {path.relative_to(inputs)}")
                                seen.add(label)

        table_s5 = inputs / "tables/accepted/Table_S5.tsv"
        if table_s5.exists():
            failures.append("saved Table S5 must be generated by make")

    # No canonical private filesystem path may survive in Git. Split-token
    # negative gates are allowed, but literal path prefixes are not.
    for path in files:
        if path.suffix in TEXT_SUFFIXES:
            payload = path.read_bytes()
            if any(token in payload for token in (b"/" + b"carc/", b"/" + b"easley/", b"/" + b"users/", b"/" + b"root/")):
                failures.append(f"private absolute path in Git: {path.relative_to(root)}")

    forbidden_public_objects = {
        "scripts/release/build_reader_tables.py",
        "scripts/release/test_release_boundary.py",
        "scripts/analysis/catalytic/validate_catalytic_sites.py",
        "fixtures/release_boundary_cases.json",
        "docs/artifact_producer_crosswalk.tsv",
        "docs/authority_inventory.tsv",
        "docs/portability_adaptations.tsv",
        "docs/source_admission.tsv",
        "docs/stale_source_exclusions.tsv",
        "docs/table_producer_crosswalk.tsv",
        "docs/PHASE1_REPORT.md", "docs/PHASE2_REPORT.md", "docs/authority_inventory.tsv",
        "docs/source_admission.tsv", "docs/stale_source_exclusions.tsv",
        "scripts/release/build_coauthor_bundle.py", "authority/manuscript/frozen_draft.md",
    }
    present_forbidden = sorted(relative for relative in forbidden_public_objects if (root / relative).exists())
    if present_forbidden:
        failures.append(f"internal/stale Git objects present: {present_forbidden}")
    if digest(logical_path) != logical_digest_before:
        failures.append("tracked logical-input manifest mutated during verification")

    for relative, claims in REQUIRED_CLAIMS.items():
        path = root / relative
        text = path.read_text(encoding="utf-8") if path.is_file() else ""
        for claim in claims:
            if not literal_claim_present(text, claim):
                failures.append(f"missing bounded reproducibility claim in {path.relative_to(root)}: {claim}")

    if failures:
        for failure in failures:
            print(f"FAIL\t{failure}")
        raise SystemExit(1)
    if args.repository_only:
        print(f"REPOSITORY_POLICY_PASS git_files={len(files)}")
    else:
        print(f"VERIFY_PASS git_files={len(files)} logical_inputs={len(logical_manifest)} data_manifest={len(read_tsv(manifest_path))} figure_inputs={len(figure_manifest)} authorities=20 Figure_S5=provenance_only")


if __name__ == "__main__":
    main()
