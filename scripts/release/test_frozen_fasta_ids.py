#!/usr/bin/env python3
"""Small regression fixture for the frozen raw FASTA exception."""

import hashlib
import tempfile
from pathlib import Path

from frozen_fasta_ids import scanner_ids, source_ids, validate_raw


def expect_failure(call):
    try:
        call()
    except ValueError:
        return
    raise AssertionError("invalid FASTA passed")


with tempfile.TemporaryDirectory() as directory:
    root = Path(directory)
    unique = root / "unique.faa"
    unique.write_text(">A.1\nMK\n>B.1\nMK\n")
    assert scanner_ids(source_ids(unique)) == ["A.1", "B.1"]
    raw = root / "raw.faa"
    raw.write_text(">A.1\nMK\n>A.1\nMK\n>B.1\nMK\n")
    digest = hashlib.sha256(raw.read_bytes()).hexdigest()
    assert scanner_ids(source_ids(raw)) == ["A.1", "A.1_1", "B.1"]
    assert validate_raw(raw, digest) == (3, 2, 3)
    expect_failure(lambda: validate_raw(raw, "0" * 64))
    expect_failure(lambda: scanner_ids(["A.1", "A.1", "A.1_1"]))
    scan = root / "scan.tsv"
    scan.write_text("Target Name\tSubfam Name\nA.1_1\tGH1_e72\n")
    scan_digest = hashlib.sha256(scan.read_bytes()).hexdigest()
    assert validate_raw(raw, digest, scan, scan_digest) == (3, 2, 3)
    scan.write_text("Target Name\tSubfam Name\nC.1\tGH1_e72\n")
    expect_failure(lambda: validate_raw(raw, digest, scan,
                                        hashlib.sha256(scan.read_bytes()).hexdigest()))
print("FROZEN_FASTA_ID_FIXTURE_PASS")
