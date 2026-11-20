#!/usr/bin/env python3
"""Run the native synteny backend from pinned frozen inputs and existing IDs."""

import sys as _source_sys
from pathlib import Path as _SourcePath
_source_repo = next(p for p in _SourcePath(__file__).resolve().parents if (p / 'config/logical_inputs.tsv').is_file())
_source_sys.path.insert(0, str(_source_repo / 'scripts/release'))
from source_aliases import source_path
import argparse
import csv
import hashlib
import importlib.util
import os
from pathlib import Path

import pandas as pd

from existing_hit_targets import accepted_order_edges, accepted_without_order, extend_native_targets


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1048576), b""):
            h.update(chunk)
    return h.hexdigest()


def pinned_sources(project, manifest):
    with manifest.open(newline="") as handle:
        records = list(csv.DictReader(handle, delimiter="\t"))
    if not records or set(records[0]) != {"role", "project_relative_path", "sha256"}:
        raise ValueError("source manifest schema changed")
    if len({r["role"] for r in records}) != len(records):
        raise ValueError("duplicate source role")
    found = {}
    project = project.resolve(strict=True)
    for r in records:
        rel = Path(r["project_relative_path"])
        if rel.is_absolute() or ".." in rel.parts:
            raise ValueError(f"escaping source path: {rel}")
        path = (source_path(project, rel)).resolve(strict=True)
        if not path.is_file() or (not str(rel).startswith("source_inputs/") and not path.is_relative_to(project)):
            raise ValueError(f"invalid source path: {rel}")
        if sha(path) != r["sha256"]:
            raise ValueError(f"source pin mismatch: {r['role']}")
        found[r["role"]] = path
    needed = {"final_gene_order", "genome_species", "missing_order", "reconstructed_order",
              "protein_to_chromosome", "backend_baseline_loci"}
    needed |= {"corrected_hit_" + family for family in ("GH5_2", "PL1_8", "PL4_4", "CBM13")}
    if needed - set(found):
        raise ValueError(f"missing native sources: {sorted(needed - set(found))}")
    return found


def proteins(path):
    frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if not {"protein_ids", "family"}.issubset(frame):
        raise ValueError(f"native loci schema changed: {path}")
    members = {}
    for row in frame.itertuples(index=False):
        for protein in row.protein_ids.split(","):
            if protein:
                if protein in members:
                    raise ValueError(f"duplicate native target: {protein}")
                members[protein] = row.family
    return members


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--project", type=Path, required=True)
    p.add_argument("--sources", type=Path, required=True)
    p.add_argument("--identity", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    inputs = pinned_sources(args.project, args.sources)
    with args.identity.open(newline="") as handle:
        identity = list(csv.DictReader(handle, delimiter="\t"))
    if not identity or len({(r["genome"], r["accepted_id"]) for r in identity}) != len(identity):
        raise ValueError("invalid accepted identity disposition")
    expected = {r["order_id"]: ("CBM13" if r["family"] == "B-lectin" else r["family"])
                for r in identity if r["order_id"]}
    if len(expected) != sum(bool(r["order_id"]) for r in identity):
        raise ValueError("two accepted copies map to one backend identity")
    hits = [inputs["corrected_hit_" + fam] for fam in ("GH5_2", "PL1_8", "PL4_4", "CBM13")]
    if len({x.parent for x in hits}) != 1:
        raise ValueError("native hit tables are not in one pinned directory")
    os.environ["PCWDE_PROJECT_ROOT"] = str(args.project.resolve())
    source = Path(__file__).with_name("01_build_synteny.py")
    spec = importlib.util.spec_from_file_location("synteny_native_existing", source)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    mod.BASE = args.project.resolve()
    mod.PROTEIN_TO_CHR_FILE = inputs["protein_to_chromosome"]
    mod.CHR_CYNOG_DICT = inputs["final_gene_order"]
    mod.DBCAN_HITS_DIR = hits[0].parent
    mod.GCA_SPECIES = inputs["genome_species"]
    mod.MISSING_ORDER_FILE = str(inputs["missing_order"])
    mod.RECONSTRUCTED_ORDER_FILE = str(inputs["reconstructed_order"])
    mod.EXCLUDED_GENOMES = {"Pediaspis_aceris"}  # existing native superseded-key rule
    mod.OUTPUT_DIR = args.out.resolve()
    mod.OUTPUT_DIR.parent.mkdir(parents=True, exist_ok=True)
    mod.OUTPUT_DIR.mkdir(parents=True, exist_ok=False)
    loader = mod.load_dbcan_target_proteins
    edges = accepted_order_edges(args.identity)
    without_order = accepted_without_order(args.identity)
    mod.INELIGIBLE_ACCEPTED_IDS = set(without_order)
    mod.load_dbcan_target_proteins = lambda: extend_native_targets(loader, edges, without_order)
    mod.main()
    observed = proteins(mod.OUTPUT_DIR / "pcwde_loci.tsv")
    if observed != expected:
        missing = sorted(set(expected) - set(observed))
        extra = sorted(set(observed) - set(expected))
        wrong = sorted(pid for pid in set(expected) & set(observed) if expected[pid] != observed[pid])
        raise ValueError(f"native eligible identity mismatch: missing={missing[:10]}, extra={extra[:10]}, wrong={wrong[:10]}")
    baseline = proteins(inputs["backend_baseline_loci"])
    for fam in ("GH5_2", "PL1_8", "PL4_4"):
        old = {pid for pid, family in baseline.items() if family == fam}
        new = {pid for pid, family in observed.items() if family == fam}
        if old != new:
            raise ValueError(f"ordinary-family target identity changed: {fam}")
    print(f"EXISTING_BACKEND_PASS targets={len(observed)} B_edges={len(edges)}")


if __name__ == "__main__":
    main()
