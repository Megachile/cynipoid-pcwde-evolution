#!/usr/bin/env python3
"""Probe source-derived target admission and rejected family collisions."""
import csv
import tempfile
from pathlib import Path

from existing_hit_targets import accepted_order_edges, accepted_without_order, extend_native_targets

with tempfile.TemporaryDirectory() as directory:
    source = Path(directory) / "identity.tsv"
    with source.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["genome", "accepted_id", "family", "order_id"],
                                delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows([
            {"genome": "G", "accepted_id": "A", "family": "B-lectin", "order_id": "M"},
            {"genome": "G", "accepted_id": "B", "family": "B-lectin", "order_id": ""},
            {"genome": "G", "accepted_id": "C", "family": "GH5_2", "order_id": "C"},
        ])
    edges = accepted_order_edges(source)
    missing = accepted_without_order(source)
    def loader():
        return {"CBM13": {"B"}, "GH5_2": {"C"}}, {"B": "CBM13", "C": "GH5_2"}
    families, by_id = extend_native_targets(loader, edges, missing)
    assert families == {"CBM13": {"M"}, "GH5_2": {"C"}}
    assert by_id == {"M": "CBM13", "C": "GH5_2"}
    try:
        extend_native_targets(lambda: ({"CBM13": set(), "GH5_2": {"M"}}, {"M": "CBM13"}),
                              edges, missing)
    except ValueError as exc:
        assert "another native family" in str(exc)
    else:
        raise AssertionError("cross-family model was accepted")
print("EXISTING_HIT_TARGETS_PASS")
