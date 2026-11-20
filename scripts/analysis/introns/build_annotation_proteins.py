#!/usr/bin/env python3
"""Build the accepted intron-mapping proteins from pinned annotation proteomes."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

FAMILIES = ("GH5_2", "PL1_8", "PL4_4")


def read(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def norm(value: str):
    return value[:-2] if value.endswith(".1") else value


def read_fasta(path: Path):
    seqs, key, parts = {}, None, []
    with path.open(errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                if key is not None:
                    seqs[key] = "".join(parts)
                key, parts = line[1:].split()[0], []
            else:
                parts.append(line.strip())
    if key is not None:
        seqs[key] = "".join(parts)
    return seqs


def hash_file(path: Path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--accounting", type=Path, required=True)
    parser.add_argument("--sources", type=Path, required=True)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--gff-proteins", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    cohort = [r for r in read(args.accounting) if r["status"] == "single_scaffold"]
    assert len(cohort) == 603 and Counter(r["family"] for r in cohort) == {"GH5_2": 176, "PL1_8": 212, "PL4_4": 215}
    sources = {r["source_file"]: r for r in read(args.sources)}
    assert len(sources) == 38
    by_source = defaultdict(list)
    for row in cohort:
        by_source[row["source_file"]].append(row)
    assert set(by_source) == set(sources)
    args.out.mkdir(parents=True, exist_ok=False)
    selected = defaultdict(dict)
    qc = []
    for source, members in sorted(by_source.items()):
        record = sources[source]
        path = args.project_root / record["proteome_snapshot_path"]
        if hash_file(path) != record["proteome_sha256"]:
            raise ValueError(f"{source}: annotation proteome hash changed")
        seqs = read_fasta(path)
        normalized = {}
        for key, seq in seqs.items():
            nkey = norm(key)
            if nkey in normalized and normalized[nkey] != seq:
                raise ValueError(f"{source}: ambiguous normalized protein {nkey}")
            normalized[nkey] = seq
        for member in members:
            pid = member["protein_id"]
            sequence = seqs.get(pid) or normalized.get(norm(pid))
            if not sequence:
                raise ValueError(f"{source}: missing annotation protein {pid}")
            selected[member["family"]][pid] = sequence.rstrip("*")
        qc.append(dict(source_file=source, selected=len(members), annotation_proteins=len(seqs)))
    assert sum(len(x) for x in selected.values()) == 603
    differences = []
    for family in FAMILIES:
        gff = read_fasta(args.gff_proteins / f"{family}_current_genome_proteins.faa")
        assert set(gff) == set(selected[family])
        for pid, sequence in selected[family].items():
            differences.append(dict(protein_id=pid, family=family, annotation_aa=len(sequence),
                                    gff_aa=len(gff[pid]), identical=sequence == gff[pid]))
        with (args.out / f"{family}_current_genome_proteins.faa").open("w") as handle:
            for pid, sequence in sorted(selected[family].items()):
                handle.write(f">{pid}\n{sequence}\n")
    with (args.out / "protein_model_comparison.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(differences[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(differences)
    summary = dict(sources=len(qc), copies=len(differences), identical=sum(r["identical"] for r in differences),
                   different=sum(not r["identical"] for r in differences))
    (args.out / "QC.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    main()
