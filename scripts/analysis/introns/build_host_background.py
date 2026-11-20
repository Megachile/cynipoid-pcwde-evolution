#!/usr/bin/env python3
"""Recompute whole-genome CDS composition and intron backgrounds per source.

Uses the same transcript and per-gene rules as the corrected-model candidate's
host metric extractor.  BK is supplied by the locked source roster, where its
NCBI assembly 31/Gnomon GFF is the only BK entry.
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path
from statistics import median

from pyfaidx import Fasta


def read(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def norm(value: str):
    value = value.strip()
    for prefix in ("transcript:", "rna-", "cds-"):
        if value.startswith(prefix):
            value = value[len(prefix):]
    return value[:-2] if value.endswith(".1") else value


def attributes(raw: str):
    return dict(bit.split("=", 1) for bit in raw.split(";") if "=" in bit)


def process(source: dict, focal: set[str], root: Path, index_dir: Path, aliases: dict[str, dict]):
    transcripts = defaultdict(list)
    transcript_gene = {}
    gff = root / source["gff_snapshot_path"]
    with gff.open(errors="replace") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            seqid, _, kind, start, end, _, strand, _, raw = fields
            attr = attributes(raw)
            if kind in ("mRNA", "transcript"):
                tid = attr.get("ID", "")
                if tid:
                    transcript_gene[tid] = attr.get("Parent", tid)
            elif kind == "CDS":
                parent = attr.get("Parent", "")
                if parent:
                    for tid in parent.split(","):
                        transcripts[tid].append((seqid, int(start), int(end), strand,
                                                 attr.get("protein_id", "")))
    assembly = root / source["assembly_snapshot_path"]
    index = index_dir / f"{source['source_file']}.fai"
    assert index.is_file(), index
    genome = Fasta(str(assembly), indexname=str(index), as_raw=True,
                   sequence_always_upper=True, rebuild=False)
    # GFF scaffold names absent from the assembly resolve only through pinned
    # exact-sequence aliases; each alias target must exist with the recorded length.
    for gff_seqid, alias in aliases.items():
        assert gff_seqid not in genome, (source["source_file"], gff_seqid, "alias shadows assembly name")
        target = alias["assembly_seqid"]
        assert target in genome and len(genome[target]) == int(alias["length"]), (source["source_file"], gff_seqid, target)
    complement = str.maketrans("ACGTN", "TGCAN")
    rows = []
    counts = dict(starting_cds_transcripts=len(transcripts), mixed_scaffold_strand=0,
                  seqid_alias_resolved=0, missing_sequence=0, short_cds=0, complete_frame=0,
                  gc3_eligible_transcripts=0, longest_cds_per_gene=0,
                  gc3_eligible_background=0)
    missing_examples = []
    for tid, pieces in transcripts.items():
        if len({(p[0], p[3]) for p in pieces}) != 1:
            counts['mixed_scaffold_strand'] += 1
            continue
        pieces.sort(key=lambda p: p[1])
        seqid, strand = pieces[0][0], pieces[0][3]
        lookup = seqid
        if seqid in aliases:
            lookup = aliases[seqid]["assembly_seqid"]
            counts['seqid_alias_resolved'] += 1
        try:
            fragments = [str(genome[lookup][start - 1:end]) for _, start, end, _, _ in pieces]
        except KeyError:
            counts['missing_sequence'] += 1
            if len(missing_examples) < 10:
                missing_examples.append((tid, seqid))
            continue
        for piece, fragment in zip(pieces, fragments):
            assert len(fragment) == piece[2] - piece[1] + 1, (source["source_file"], tid, piece)
        cds = "".join(fragments)
        if strand == "-":
            cds = cds.translate(complement)[::-1]
        length = len(cds)
        if length < 3:
            counts['short_cds'] += 1
            continue
        valid = [cds[i:i + 3] for i in range(0, length - 2, 3)
                 if set(cds[i:i + 3]) <= set("ACGT")]
        complete = length % 3 == 0
        counts['complete_frame'] += int(complete)
        terminal_stop = complete and cds[-3:] in {"TAA", "TAG", "TGA"}
        counts['gc3_eligible_transcripts'] += int(complete and bool(valid))
        called_bases = sum(base in "ACGT" for base in cds)
        gaps = [b[1] - a[2] - 1 for a, b in zip(pieces, pieces[1:]) if b[1] > a[2] + 1]
        cumulative = 0
        positions = []
        for left, right in zip(pieces, pieces[1:]):
            cumulative += left[2] - left[1] + 1
            if right[1] > left[2] + 1:
                relative = cumulative / length
                positions.append(1 - relative if strand == "-" else relative)
        ids = {norm(tid), norm(pieces[0][4])}
        rows.append(dict(source_file=source["source_file"], transcript_id=tid,
                         gene_id=transcript_gene.get(tid, tid), scaffold=seqid,
                         protein_id=pieces[0][4], cds_nt=length,
                         protein_aa=(length // 3 - int(terminal_stop)) if complete else "",
                         frame_complete=complete, intron_count=len(gaps),
                         intron_density_per_kb=1000 * len(gaps) / length,
                         median_intron_nt=median(gaps) if gaps else "",
                         gc=sum(base in "GC" for base in cds) / called_bases if called_bases else "",
                         gc3=sum(codon[2] in "GC" for codon in valid) / len(valid) if complete and valid else "",
                         valid_codons=len(valid), candidate_ids=";".join(sorted(ids)),
                         _focal_ids=ids & focal, _relative_positions=positions, _gap_lengths=gaps))
    best = {}
    for row in rows:
        gene = row["gene_id"]
        if gene not in best or row["cds_nt"] > best[gene]["cds_nt"]:
            best[gene] = row
    background = list(best.values())
    counts['longest_cds_per_gene'] = len(background)
    counts['gc3_eligible_background'] = sum(row['gc3'] != '' for row in background)
    counts['missing_sequence_examples'] = ';'.join(f'{tid}:{seqid}' for tid, seqid in missing_examples)
    assert (counts['starting_cds_transcripts'] == counts['mixed_scaffold_strand'] +
            counts['missing_sequence'] + counts['short_cds'] + len(rows))
    assert counts['gc3_eligible_background'] <= counts['longest_cds_per_gene'] <= len(rows)
    focal_rows = [row for row in rows if row["_focal_ids"]]
    found = set().union(*(row["_focal_ids"] for row in focal_rows)) if focal_rows else set()
    if found != focal:
        raise ValueError(f"{source['source_file']}: focal ID mismatch: missing={sorted(focal - found)[:10]}; extra={sorted(found - focal)[:10]}")
    return background, focal_rows, counts


def write(path: Path, rows: list[dict]):
    fields = [key for key in rows[0] if not key.startswith("_")]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--accounting", type=Path, required=True)
    parser.add_argument("--blectin-models", type=Path, required=True)
    parser.add_argument("--sources", type=Path, required=True)
    parser.add_argument("--seqid-aliases", type=Path, required=True)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--index-dir", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    aliases = defaultdict(dict)
    for row in read(args.seqid_aliases):
        assert row["gff_seqid"] not in aliases[row["source_file"]], row
        aliases[row["source_file"]][row["gff_seqid"]] = row
    for source_aliases in aliases.values():
        targets = [row["assembly_seqid"] for row in source_aliases.values()]
        assert len(targets) == len(set(targets)), "assembly alias target reused"
    focal = defaultdict(set)
    for row in read(args.accounting):
        if row["status"] == "single_scaffold":
            focal[row["source_file"]].add(norm(row["protein_id"]))
    assert sum(len(ids) for ids in focal.values()) == 603
    blectin = read(args.blectin_models)
    assert len(blectin) == 28 and len({r["protein_id"] for r in blectin}) == 28
    for row in blectin:
        focal[row["source_file"]].add(norm(row["protein_id"]))
    assert sum(len(ids) for ids in focal.values()) == 631
    qc = []
    all_positions = []
    all_densities = []
    all_lengths = []
    with (args.out / "background_intron_details.tsv").open("w", newline="") as handle:
        detail = csv.DictWriter(handle, ["source_file", "transcript_id", "intron_number", "intron_length",
                                         "relative_cds_position"], delimiter="\t", lineterminator="\n")
        detail.writeheader()
        for source in read(args.sources):
            name = source["source_file"]
            background, focal_rows, counts = process(source, focal[name], args.project_root, args.index_dir,
                                                   aliases.get(name, {}))
            for row in background:
                all_positions.extend(row["_relative_positions"])
                all_densities.append(row["intron_density_per_kb"])
                if row["median_intron_nt"] != "":
                    all_lengths.append(row["median_intron_nt"])
                for i, (length, relative) in enumerate(zip(row["_gap_lengths"], row["_relative_positions"]), 1):
                    detail.writerow(dict(source_file=name, transcript_id=row["transcript_id"], intron_number=i,
                                         intron_length=length, relative_cds_position=relative))
            write(args.out / f"{name}.tsv", [{**r, "cohort": "background"} for r in background] +
                  [{**r, "cohort": "focal"} for r in focal_rows])
            qc.append(dict(source_file=name, background=len(background), focal_rows=len(focal_rows),
                           target_ids=len(focal[name]), **counts))
    assert len(qc) == 38 and sum(row["target_ids"] for row in qc) == 631
    write(args.out / "QC.tsv", qc)
    (args.out / "QC.json").write_text(json.dumps(qc, indent=2) + "\n")
    comparator = dict(background_genes=len(all_densities), introns=len(all_positions),
                      first_cds_quarter_n=sum(value <= 0.25 for value in all_positions),
                      first_cds_quarter_fraction=sum(value <= 0.25 for value in all_positions) / len(all_positions),
                      median_relative_cds_position=median(all_positions),
                      median_intron_density_per_kb=median(all_densities),
                      median_of_gene_median_intron_nt=median(all_lengths))
    (args.out / "whole_genome_comparator.json").write_text(json.dumps(comparator, indent=2) + "\n")
    print(json.dumps({"sources": len(qc), "backgrounds": sum(r["background"] for r in qc),
                      "focal_rows": sum(r["focal_rows"] for r in qc)}, sort_keys=True))


if __name__ == "__main__":
    main()
