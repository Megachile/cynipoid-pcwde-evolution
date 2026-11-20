#!/usr/bin/env python3
"""Fail if retired saved intron inputs reenter the build contracts."""
from pathlib import Path
import csv

repo = Path(__file__).resolve().parents[2]
retired = "tables/introns_current"
for name, column in (
    ("config/closure_declared.tsv", "pattern"),
    ("config/deposit_layout.tsv", "source"),
    ("config/deposit_sources.tsv", "upstream_path"),
    ("config/build_inputs.tsv", "path"),
    ("config/logical_inputs.tsv", "figshare_relative_path"),
):
    path = repo / name
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if column not in rows[0]:
        raise SystemExit(f"missing {column} in {name}: {list(rows[0])}")
    bad = [r[column] for r in rows if retired in r[column]]
    if bad:
        raise SystemExit(f"retired intron input in {name}: {bad}")
stage = (repo / "scripts/deposit_inputs/stage.py").read_text()
if '"tables/introns_current"' not in stage or "retire_built_introns(inputs)" not in stage:
    raise SystemExit("saved intron retirement missing from stage.py")
print("INTRON_RETIREMENT_PASS closure/deposit/build contracts and stage retirement")
