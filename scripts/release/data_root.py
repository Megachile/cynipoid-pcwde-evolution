#!/usr/bin/env python3
"""Resolve and fail-closed validate the single Figshare data-root interface."""

from __future__ import annotations

import csv
import hashlib
import os
from pathlib import Path
from typing import Optional


MANIFEST_FIELDS = (
    "logical_id", "figshare_relative_path", "sha256", "bytes", "role",
    "producer_authority", "consuming_commands",
)
SOURCE_MANIFEST_SHA256 = "62b286b2b8002af0e14977e72bbb0e608989999754f1a0ffa7bd254f119b3ad3"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def resolve_data_root(repo: Path, explicit: Optional[Path] = None) -> Path:
    configured = explicit or (Path(os.environ["PCWDE_DATA_ROOT"])
                              if os.environ.get("PCWDE_DATA_ROOT") else repo / "figshare_data")
    root = configured.expanduser().resolve()
    if not root.is_dir() or root.is_symlink():
        raise SystemExit(f"missing or unsafe Figshare data root: {root}")
    return root


def validate_logical_inputs(repo: Path, data_root: Path) -> list[dict[str, str]]:
    manifest_path = repo / "config/logical_inputs.tsv"
    if manifest_path.is_symlink() or not manifest_path.is_file():
        raise SystemExit("missing or unsafe logical-input manifest")
    manifest_digest_before = sha256(manifest_path)
    with manifest_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != MANIFEST_FIELDS:
            raise SystemExit("logical-input manifest schema mismatch")
        rows = list(reader)
    if (
        not rows
        or len({row["logical_id"] for row in rows}) != len(rows)
        or len({row["figshare_relative_path"] for row in rows}) != len(rows)
    ):
        raise SystemExit("logical-input IDs or paths are empty/duplicated")
    inputs_root = data_root / "repository_inputs"
    if inputs_root.is_symlink() or not inputs_root.is_dir():
        raise SystemExit("missing or unsafe repository_inputs root")
    expected_paths = set()
    for row in rows:
        relative = Path(row["figshare_relative_path"])
        if relative.is_absolute() or ".." in relative.parts or relative.parts[0] != "repository_inputs":
            raise SystemExit(f"unsafe logical-input path: {relative}")
        path = data_root / relative
        if path.is_symlink() or not path.is_file():
            raise SystemExit(f"missing logical input: {row['logical_id']}")
        if path.stat().st_size != int(row["bytes"]) or sha256(path) != row["sha256"]:
            raise SystemExit(f"logical-input hash mismatch: {row['logical_id']}")
        expected_paths.add(relative.as_posix())
    actual_paths = {
        path.relative_to(data_root).as_posix()
        for path in inputs_root.rglob("*") if path.is_file()
    }
    # Collaborator Figure S5 is separately pinned by the manuscript build.
    actual_paths.discard("repository_inputs/figures/figure_s5/Figure_S5.pdf")
    if actual_paths != expected_paths:
        raise SystemExit(
            f"logical-input closure mismatch missing={sorted(expected_paths-actual_paths)[:5]} "
            f"extra={sorted(actual_paths-expected_paths)[:5]}"
        )
    if sha256(manifest_path) != manifest_digest_before:
        raise SystemExit("logical-input manifest mutated during validation")
    from public_scope_policy import validate
    validate(inputs_root, repo)
    return rows


def write_projected_inputs(data_root: Path) -> None:
    """Bind the transformed public input tree without altering source pins."""
    root = data_root / "repository_inputs"
    target = data_root / "PROJECTED_INPUTS.tsv"
    if target.exists() or target.is_symlink():
        raise SystemExit("refusing to overwrite projected-input manifest")
    paths = sorted(p for p in root.rglob("*") if p.is_file())
    if any(p.is_symlink() for p in root.rglob("*")):
        raise SystemExit("projected inputs contain symlinks")
    with target.open("x", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("relative_path", "bytes", "sha256"))
        for path in paths:
            writer.writerow((path.relative_to(data_root).as_posix(), path.stat().st_size, sha256(path)))


def validate_projected_inputs(data_root: Path) -> int:
    root = data_root / "repository_inputs"
    manifest = data_root / "PROJECTED_INPUTS.tsv"
    if manifest.is_symlink() or not manifest.is_file() or root.is_symlink() or not root.is_dir():
        raise SystemExit("missing or unsafe projected-input manifest")
    with manifest.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["relative_path", "bytes", "sha256"]:
            raise SystemExit("projected-input manifest schema mismatch")
        rows = list(reader)
    expected = set()
    for row in rows:
        relative = Path(row["relative_path"])
        if relative.is_absolute() or ".." in relative.parts or relative.parts[0] != "repository_inputs":
            raise SystemExit("unsafe projected-input path")
        path = data_root / relative
        if path.is_symlink() or not path.is_file() or path.stat().st_size != int(row["bytes"]) or sha256(path) != row["sha256"]:
            raise SystemExit(f"projected-input mismatch: {relative}")
        expected.add(relative.as_posix())
    actual = {p.relative_to(data_root).as_posix() for p in root.rglob("*") if p.is_file()}
    if len(rows) != len(expected) or actual != expected or any(p.is_symlink() for p in root.rglob("*")):
        raise SystemExit("projected-input closure mismatch")
    return len(rows)


def validate_reader_inputs(repo: Path, data_root: Path) -> list[dict[str, str]]:
    """Validate either raw pinned inputs or a complete projected reader tree."""
    manifest = repo / "config/logical_inputs.tsv"
    if manifest.is_symlink() or not manifest.is_file() or sha256(manifest) != SOURCE_MANIFEST_SHA256:
        raise SystemExit("source-pin identity mismatch")
    with manifest.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != MANIFEST_FIELDS:
            raise SystemExit("logical-input manifest schema mismatch")
        rows = list(reader)
    expected = {row["figshare_relative_path"] for row in rows}
    if len(rows) != 480 or len(expected) != len(rows):
        raise SystemExit("source-pin path inventory mismatch")
    projected = data_root / "PROJECTED_INPUTS.tsv"
    if not projected.exists() and not projected.is_symlink():
        return validate_logical_inputs(repo, data_root)
    validate_projected_inputs(data_root)
    with projected.open(newline="", encoding="utf-8") as handle:
        projected_rows = list(csv.DictReader(handle, delimiter="\t"))
    if {row["relative_path"] for row in projected_rows} != expected:
        raise SystemExit("projected/source path inventory mismatch")
    from public_inputs import public_problem_set
    problems = public_problem_set(data_root / "repository_inputs")
    if problems:
        raise SystemExit(f"mixed or unsafe projected input form: {sorted(problems)[:5]}")
    from public_scope_policy import validate
    validate(data_root / "repository_inputs", repo)
    return rows
