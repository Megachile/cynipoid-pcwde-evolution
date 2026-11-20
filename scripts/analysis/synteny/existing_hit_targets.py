"""Add accepted B-lectin identities with existing gene-order positions to native hits."""
import csv
from pathlib import Path


def accepted_order_edges(identity_path: Path):
    with identity_path.open(newline="") as handle:
        records = list(csv.DictReader(handle, delimiter="\t"))
    if len({(r["genome"], r["accepted_id"]) for r in records}) != len(records):
        raise ValueError("duplicate accepted identity")
    edges = {}
    for r in records:
        if r["family"] != "B-lectin" or not r["order_id"]:
            continue
        if r["order_id"] in edges and edges[r["order_id"]] != r["accepted_id"]:
            raise ValueError(f"two accepted units map to one ordered model: {r['order_id']}")
        edges[r["order_id"]] = r["accepted_id"]
    return edges


def accepted_without_order(identity_path: Path):
    with identity_path.open(newline="") as handle:
        records = list(csv.DictReader(handle, delimiter="\t"))
    return {r["accepted_id"] for r in records if not r["order_id"]}


def extend_native_targets(original_loader, edges, without_order=()):
    by_family, family_by_id = original_loader()
    conflicts = sorted(pid for pid in edges
                       if family_by_id.get(pid, "CBM13") != "CBM13"
                       or any(pid in ids for family, ids in by_family.items() if family != "CBM13"))
    if conflicts:
        raise ValueError(f"accepted B-lectin model has another native family: {conflicts}")
    if "CBM13" not in by_family:
        raise ValueError("native B-lectin family missing")
    # A native hit, or a supplemental backend target, is not by itself an
    # existing gene-order position. Apply this rule to every accepted family.
    for pid in without_order:
        present = {family for family, ids in by_family.items() if pid in ids}
        if len(present) > 1:
            raise ValueError(f"unpositioned accepted identity crosses native families: {pid}")
        for family in present:
            by_family[family].remove(pid)
        family_by_id.pop(pid, None)
    for pid in edges:
        by_family["CBM13"].add(pid)
        family_by_id[pid] = "CBM13"
    return by_family, family_by_id
