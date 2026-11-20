#!/usr/bin/env python3
"""Resolve the paper's three-family member IDs to pinned GFF CDS segments.

The member table defines the 629-copy universe.  The genome/source crosswalk and
source manifest are explicit, versioned inputs.  Multi-scaffold or mixed-strand
models remain in the accounting output and cannot enter intron calculations.
"""
from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

FAMILIES = {"GH5_2", "PL1_8", "PL4_4"}


def rows(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def norm(value: str) -> str:
    value = value.strip()
    for prefix in ("transcript:", "rna-", "cds-"):
        if value.startswith(prefix):
            value = value[len(prefix):]
    return value[:-2] if value.endswith(".1") else value


def source_alias(source: str, protein: str) -> str:
    if source == "DL_bkragtag":
        import re
        protein = protein.replace("DL_bkragtag_", "D_quercuslanigerum_", 1)
        protein = re.sub(r"(\.1)_[12]$", r"\1", protein)
    return protein


def parse_gff(path: Path, aliases: dict[str, str]):
    features = defaultdict(list)
    with path.open(errors="replace") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "CDS":
                continue
            attrs = dict(part.split("=", 1) for part in fields[8].split(";") if "=" in part)
            matched = None
            for key in ("protein_id", "Parent", "ID"):
                for item in attrs.get(key, "").split(","):
                    matched = aliases.get(norm(item))
                    if matched:
                        break
                if matched:
                    break
            if matched:
                start, end = int(fields[3]), int(fields[4])
                assert 0 < start <= end, (path, matched, start, end)
                features[matched].append((fields[0], fields[6], start, end, fields[7]))
    return features


def build(members: list[dict], genome_sources: list[dict], sources: list[dict], project_root: Path):
    member_rows = [r for r in members if r["family"] in FAMILIES]
    assert len(member_rows) == 629 and len({r["protein_id"] for r in member_rows}) == 629
    crosswalk = {r["genome"]: r["source_file"] for r in genome_sources}
    assert len(crosswalk) == len(genome_sources) == 38
    source_paths = {r["source_file"]: project_root / r["gff_snapshot_path"] for r in sources}
    assert len(source_paths) == len(sources) == 38
    assert {r["genome"] for r in member_rows} == set(crosswalk)
    assert set(crosswalk.values()) == set(source_paths)
    by_source = defaultdict(list)
    for member in member_rows:
        by_source[crosswalk[member["genome"]]].append(member)
    accounting, segments = [], []
    for source in sorted(by_source):
        aliases = {}
        for member in by_source[source]:
            protein = member["protein_id"]
            # The retained DL split-refreeze GFF uses the placement IDs, while
            # older producer records use D_quercuslanigerum IDs for some copies.
            for alias in (protein, source_alias(source, protein)):
                key = norm(alias)
                if key in aliases and aliases[key] != protein:
                    raise ValueError(f"{source}: ambiguous CDS alias {key}")
                aliases[key] = protein
        features = parse_gff(source_paths[source], aliases)
        for member in by_source[source]:
            protein = member["protein_id"]
            cds = features.get(protein, [])
            scaffolds = {x[0] for x in cds}
            strands = {x[1] for x in cds}
            status = "single_scaffold" if cds and len(scaffolds) == len(strands) == 1 else (
                "missing_CDS" if not cds else "cross_scaffold_or_strand")
            accounting.append(dict(protein_id=protein, family=member["family"], genome=member["genome"],
                                   source_file=source, status=status, cds_segments=len(cds),
                                   scaffold=";".join(sorted(scaffolds)), strand=";".join(sorted(strands))))
            if status == "single_scaffold":
                ordered = sorted(cds, key=lambda x: x[2], reverse=next(iter(strands)) == "-")
                for i, (scaffold, strand, start, end, phase) in enumerate(ordered, 1):
                    segments.append(dict(protein_id=protein, family=member["family"], genome=member["genome"],
                                         source_file=source, scaffold=scaffold, strand=strand,
                                         transcript_order=i, cds_start=start, cds_end=end, gff_phase=phase))
    assert len(accounting) == 629 and len({r["protein_id"] for r in accounting}) == 629
    return accounting, segments


def build_blectin(models: list[dict], sources: list[dict], project_root: Path):
    assert len(models) == len({r["protein_id"] for r in models}) == 28
    assert len({r["species"] for r in models}) == 16
    paths = {r["source_file"]: project_root / r["gff_snapshot_path"] for r in sources}
    by_source = defaultdict(list)
    for model in models:
        by_source[model["source_file"]].append(model)
    assert set(by_source) <= set(paths)
    accounting, segments = [], []
    for source, source_models in sorted(by_source.items()):
        aliases = {norm(r["protein_id"]): r["protein_id"] for r in source_models}
        assert len(aliases) == len(source_models), source
        features = parse_gff(paths[source], aliases)
        for model in source_models:
            protein = model["protein_id"]
            cds = features.get(protein, [])
            scaffolds = {x[0] for x in cds}
            strands = {x[1] for x in cds}
            if not cds or len(scaffolds) != 1 or len(strands) != 1:
                raise ValueError(f"{source} {protein}: unresolved coherent B-lectin CDS")
            accounting.append(dict(protein_id=protein, family="B_lectin", genome=model["genome"],
                                   source_file=source, status="single_scaffold", cds_segments=len(cds),
                                   scaffold=next(iter(scaffolds)), strand=next(iter(strands)),
                                   species=model["species"], reference_block=model["reference_block"],
                                   placement_id=model["placement_id"]))
            ordered = sorted(cds, key=lambda x: x[2], reverse=next(iter(strands)) == "-")
            for i, (scaffold, strand, start, end, phase) in enumerate(ordered, 1):
                segments.append(dict(protein_id=protein, family="B_lectin", genome=model["genome"],
                                     source_file=source, scaffold=scaffold, strand=strand,
                                     transcript_order=i, cds_start=start, cds_end=end, gff_phase=phase))
    assert len(accounting) == 28
    return accounting, segments


def write(path: Path, data: list[dict]):
    assert data
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(data[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(data)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--members", type=Path, required=True)
    parser.add_argument("--genome-sources", type=Path, required=True)
    parser.add_argument("--sources", type=Path, required=True)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--blectin-models", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    accounting, segments = build(rows(args.members), rows(args.genome_sources), rows(args.sources), args.project_root)
    write(args.out / "cohort_accounting.tsv", accounting)
    write(args.out / "member_cds_segments.tsv", segments)
    if args.blectin_models:
        b_accounting, b_segments = build_blectin(rows(args.blectin_models), rows(args.sources), args.project_root)
        write(args.out / "blectin_accounting.tsv", b_accounting)
        write(args.out / "blectin_cds_segments.tsv", b_segments)
    from collections import Counter
    print(dict(Counter(r["status"] for r in accounting)))


if __name__ == "__main__":
    main()
