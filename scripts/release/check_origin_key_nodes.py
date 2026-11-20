#!/usr/bin/env python3
"""Check that plotted Figure 2 key labels and published numbers use producer values."""
import argparse
import csv
import json
from pathlib import Path

p = argparse.ArgumentParser()
p.add_argument("--build", type=Path, required=True)
a = p.parse_args()

def rows(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))

key = {(r["family"], r["node"]): r for r in rows(a.build / "trees/origin_key_nodes.tsv") if r["view"] == "fig2"}
drawn = rows(a.build / "figures/origin_trees/evidence/figure2_key_node_labels.tsv")
numbers = json.loads((a.build / "numbers.json").read_text())["origin_key_nodes"]
if len(key) != 6 or len(drawn) != 6 or len({(r["family"], r["node"]) for r in drawn}) != 6:
    raise SystemExit("origin key node coverage differs from six Figure 2 labels")
for row in drawn:
    family, role = row["family"], row["node"]
    expected = key[(family, role)]["support"]
    if row["drawn_support"] != expected or numbers[family][role]["support"] != expected:
        raise SystemExit(f"origin support mismatch at {family} {role}")
print("ORIGIN_KEY_CONSISTENCY_PASS labels=6")
