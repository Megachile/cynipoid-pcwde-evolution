"""Validate the two frozen raw dbCAN-sub FASTAs with duplicate source IDs."""

import csv
import hashlib
from collections import Counter
from pathlib import Path


DECLARED = {
    "08c7e24f015f77e65ad98abf8bfc9c68387d08782f46e4d9c739d13eada1296f",
    "e863b115092468f73aebfbb0e40971654d3a4b445e29cff7e940701463c0aca1",
}


def source_ids(path):
    with Path(path).open(encoding="utf-8") as handle:
        return [line[1:].split()[0] for line in handle if line.startswith(">")]


def scanner_ids(ids):
    """dbCAN keeps the first ID and numbers later occurrences from _1."""
    counts = Counter()
    result = []
    for source_id in ids:
        if not source_id:
            raise ValueError("empty FASTA ID")
        n = counts[source_id]
        result.append(source_id if n == 0 else f"{source_id}_{n}")
        counts[source_id] += 1
    if len(result) != len(set(result)):
        raise ValueError("processed FASTA ID collision")
    return result


def declared_raw_sources(repo):
    repo = Path(repo)
    with (repo / "config/cazy_scans.tsv").open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    selected = {}
    for row in rows:
        digest = row["proteome_sha256"]
        if digest not in DECLARED:
            continue
        if row["mapping_method"] != "dbcan_duplicate_suffix":
            raise ValueError(f"undeclared scanner ID transform: {digest}")
        expected = row["proteome_path"]
        if not expected.startswith("repository_inputs/cazy_inputs/proteome/") or ".." in Path(expected).parts:
            raise ValueError("unsafe frozen raw FASTA registration")
        if row["proteome_path"] != expected or expected in selected:
            raise ValueError(f"unexpected frozen raw FASTA registration: {digest}")
        scan_path = row["scan_path"]
        if not scan_path.startswith("repository_inputs/cazy_inputs/scan/"):
            raise ValueError(f"unexpected frozen scan registration: {digest}")
        selected[expected] = (digest, scan_path, row["scan_sha256"])
    if {item[0] for item in selected.values()} != DECLARED:
        raise ValueError("incomplete frozen raw FASTA registrations")
    return selected


def validate_raw(path, expected_digest, scan_path=None, scan_digest=None):
    path = Path(path)
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != expected_digest:
        raise ValueError("frozen raw FASTA hash mismatch")
    ids = source_ids(path)
    if not ids:
        raise ValueError("empty frozen raw FASTA")
    first = scanner_ids(ids)
    if scan_path is not None:
        scan_path = Path(scan_path)
        if hashlib.sha256(scan_path.read_bytes()).hexdigest() != scan_digest:
            raise ValueError("frozen scan hash mismatch")
        with scan_path.open(newline="", encoding="utf-8") as handle:
            hits = {row["Target Name"] for row in csv.DictReader(handle, delimiter="\t")}
        missing = hits - set(first)
        if missing:
            raise ValueError(f"scan targets outside processed FASTA: {len(missing)}")
    return len(ids), len(set(ids)), len(first)
