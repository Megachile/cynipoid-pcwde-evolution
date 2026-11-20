#!/usr/bin/env python3
"""Build three-family CDS proteins, introns, and GC3 from pinned assemblies.

Input CDS segments are derived from pinned GFFs by build_member_cds.py.  Every
single-scaffold member is retained, including zero-intron copies.  Cohort and
sequence checks fail before these tables can be used by a paper build.
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

from Bio.Seq import Seq
from pyfaidx import Fasta

FAMILIES = ("GH5_2", "PL1_8", "PL4_4")
REFERENCES = {"GH5_2": "XP_033220919", "PL1_8": "XP_033225824", "PL4_4": "XP_033222742"}
STOPS = {"TAA", "TAG", "TGA"}


def read(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write(path: Path, rows: list[dict], columns: list[str]):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def reverse_complement(seq: str) -> str:
    return str(Seq(seq).reverse_complement())


def read_fasta(path: Path):
    sequences, key, chunks = {}, None, []
    with path.open() as handle:
        for line in handle:
            if line.startswith('>'):
                if key is not None:
                    assert key not in sequences
                    sequences[key] = ''.join(chunks)
                key, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
    if key is not None:
        assert key not in sequences
        sequences[key] = ''.join(chunks)
    return sequences


def introns(segments: list[dict], genome=None, junctions=None):
    """Use the candidate's transcript-order phase_before convention."""
    ordered = sorted(segments, key=lambda x: int(x["transcript_order"]))
    strand = ordered[0]["strand"]
    total = sum(int(x["cds_end"]) - int(x["cds_start"]) + 1 for x in ordered)
    initial_phase = int(ordered[0]["gff_phase"]) if ordered[0]["gff_phase"] in ("0", "1", "2") else 0
    coding_total = total - initial_phase
    cumulative = 0
    result = []
    for i, (left, right) in enumerate(zip(ordered, ordered[1:]), 1):
        left_start, left_end = int(left["cds_start"]), int(left["cds_end"])
        right_start, right_end = int(right["cds_start"]), int(right["cds_end"])
        cumulative += left_end - left_start + 1
        intron_start, intron_end = (left_end + 1, right_start - 1) if strand == "+" else (right_end + 1, left_start - 1)
        if intron_end < intron_start:
            continue
        if junctions is None:
            seq = str(genome[left["scaffold"]][intron_start - 1:intron_end]).upper()
            if strand == "-":
                seq = reverse_complement(seq)
            donor, acceptor = seq[:2], seq[-2:]
        else:
            donor, acceptor = junctions[i]
            assert len(donor) == len(acceptor) == min(2, intron_end - intron_start + 1)
        phase = (cumulative - initial_phase) % 3
        result.append(dict(intron_number=i, start=intron_start, end=intron_end,
                           length=intron_end - intron_start + 1,
                           phase_before=phase,
                           phase_after=int(right["gff_phase"]) if right["gff_phase"] in ("0", "1", "2") else 0,
                           phase=phase,
                           relative_cds_position=round((cumulative - initial_phase) / coding_total, 4),
                           donor_dinuc=donor, acceptor_dinuc=acceptor,
                           aa_position=round((cumulative - initial_phase) / 3)))
    return result


def build(accounting: list[dict], segments: list[dict], sources: list[dict], project_root: Path,
          index_dir: Path, universe_expected: int, eligible_expected: int,
          frozen_cds=None, frozen_junctions=None):
    source_paths = {r["source_file"]: project_root / r["assembly_snapshot_path"] for r in sources}
    by_pid = defaultdict(list)
    for segment in segments:
        by_pid[segment["protein_id"]].append(segment)
    by_source = defaultdict(list)
    for row in accounting:
        if row["status"] == "single_scaffold":
            by_source[row["source_file"]].append(row)
    assert len(accounting) == universe_expected and sum(len(x) for x in by_source.values()) == eligible_expected
    observations, metrics, proteins = [], [], defaultdict(dict)
    index_dir.mkdir(parents=True, exist_ok=True)
    for source in sorted(by_source):
        if frozen_cds is None:
            assembly = source_paths[source]
            assert assembly.is_file(), assembly
            # Existing shared .fai files can be stale even when the FASTA and GFF
            # are pinned. Build indexes in task-owned storage from the FASTA bytes.
            genome = Fasta(str(assembly), indexname=str(index_dir / f"{source}.fai"),
                           rebuild=True, as_raw=True, sequence_always_upper=True)
        else:
            genome = None
        for member in by_source[source]:
            pid = member["protein_id"]
            cds = sorted(by_pid[pid], key=lambda x: int(x["transcript_order"]))
            assert len(cds) == int(member["cds_segments"]) and cds
            assert {x["scaffold"] for x in cds} == {member["scaffold"]}
            assert {x["strand"] for x in cds} == {member["strand"]}
            if frozen_cds is None:
                parts = [str(genome[x["scaffold"]][int(x["cds_start"]) - 1:int(x["cds_end"])]).upper() for x in cds]
                for segment, part in zip(cds, parts):
                    expected = int(segment["cds_end"]) - int(segment["cds_start"]) + 1
                    assert len(part) == expected, (pid, segment["scaffold"], expected, len(part))
                if member["strand"] == "-":
                    parts = [reverse_complement(part) for part in parts]
                sequence = "".join(parts)
            else:
                sequence = frozen_cds[pid]
                assert len(sequence) == sum(int(x['cds_end']) - int(x['cds_start']) + 1 for x in cds)
            initial_phase = int(cds[0]["gff_phase"]) if cds[0]["gff_phase"] in ("0", "1", "2") else 0
            coding = sequence[initial_phase:]
            coding_codons = [coding[i:i + 3] for i in range(0, len(coding) - 2, 3)]
            protein = str(Seq("".join(coding_codons)).translate())
            internal_stops = protein[:-1].count("*")
            assert internal_stops == 0, (pid, internal_stops)
            proteins[member["family"]][pid] = protein.rstrip("*")
            # The accepted host-composition extractor counts raw GFF CDS
            # triplets; its GC3 reference does not trim the initial GFF phase.
            codons = [sequence[i:i + 3] for i in range(0, len(sequence) - 2, 3)]
            valid = [codon for codon in codons if set(codon) <= set("ACGT")]
            coding_valid = [codon for codon in coding_codons if set(codon) <= set("ACGT")]
            called = sum(base in "ACGT" for base in sequence)
            xs = introns(cds, genome, frozen_junctions.get(pid) if frozen_junctions is not None else None)
            metrics.append(dict(protein_id=pid, family=member["family"], genome=member["genome"],
                                source_file=source, cds_nt=len(sequence), trailing_nt=len(sequence) % 3,
                                internal_stops=internal_stops, gc=sum(base in "GC" for base in sequence) / called if called else "",
                                gc3=sum(codon[2] in "GC" for codon in valid) / len(valid) if valid and len(sequence) % 3 == 0 else "",
                                gc3_phase_corrected=sum(codon[2] in "GC" for codon in coding_valid) / len(coding_valid) if coding_valid else "",
                                intron_count=len(xs), intron_density_per_kb=1000 * len(xs) / len(sequence)))
            base = dict(protein_id=pid, family=member["family"], genome=member["genome"],
                        source_file=source, scaffold=member["scaffold"], strand=member["strand"], intron_count=len(xs))
            if not xs:
                observations.append(base | dict(intron_number=0, start=0, end=0, length=0, phase_before=-1,
                                                phase_after=-1, phase=-1, relative_cds_position=0,
                                                donor_dinuc="", acceptor_dinuc="", aa_position=0))
            for x in xs:
                observations.append(base | {k: x[k] for k in ("intron_number", "start", "end", "length", "phase_before",
                                                          "phase_after", "phase", "relative_cds_position",
                                                          "donor_dinuc", "acceptor_dinuc", "aa_position")})
    assert len(metrics) == eligible_expected and len({r["protein_id"] for r in observations}) == eligible_expected
    return observations, metrics, proteins


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--accounting", type=Path, required=True)
    parser.add_argument("--segments", type=Path, required=True)
    parser.add_argument("--member-table", type=Path, required=True)
    parser.add_argument("--blectin-accounting", type=Path)
    parser.add_argument("--blectin-segments", type=Path)
    parser.add_argument("--sources", type=Path, required=True)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--index-dir", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--focal-cds", type=Path)
    parser.add_argument("--focal-junctions", type=Path)
    args = parser.parse_args()
    assert bool(args.focal_cds) == bool(args.focal_junctions)
    args.out.mkdir(parents=True, exist_ok=False)
    source_rows = read(args.sources)
    frozen_cds = read_fasta(args.focal_cds) if args.focal_cds else None
    frozen_junctions = None
    if args.focal_junctions:
        frozen_junctions = defaultdict(dict)
        for row in read(args.focal_junctions):
            pid, number = row['protein_id'], int(row['intron_number'])
            assert number not in frozen_junctions[pid]
            frozen_junctions[pid][number] = (row['donor_dinuc'], row['acceptor_dinuc'])
        assert len(frozen_cds) == 631
    observations, metrics, proteins = build(read(args.accounting), read(args.segments), source_rows,
                                            args.project_root, args.index_dir, 629, 603,
                                            frozen_cds, frozen_junctions)
    write(args.out / "per_intron_observations.tsv", observations, list(observations[0]))
    write(args.out / "focal_composition.tsv", metrics, list(metrics[0]))
    for family in FAMILIES:
        family_rows = [row for row in observations if row["family"] == family]
        write(args.out / f"{family}_intron_metadata.tsv", family_rows, list(observations[0]))
        with (args.out / f"{family}_current_genome_proteins.faa").open("w") as handle:
            for pid, sequence in sorted(proteins[family].items()):
                handle.write(f">{pid}\n{sequence}\n")
    members = {r["protein_id"]: r for r in read(args.member_table)}
    reference_rows = []
    for family, pid in REFERENCES.items():
        member = members[pid]
        chromosome = member["bk_chromosome"]
        reference_rows.append(dict(family=family, ref_id=pid, ref_species=member["species"],
                                   ref_chromosome=chromosome, ref_locus_class=chromosome,
                                   ref_len_full=len(proteins[family][pid])))
    write(args.out / "ref_proteins.tsv", reference_rows, list(reference_rows[0]))
    if bool(args.blectin_accounting) != bool(args.blectin_segments):
        raise ValueError("B-lectin accounting and segments must be supplied together")
    if args.blectin_accounting:
        b_observations, b_metrics, b_proteins = build(read(args.blectin_accounting), read(args.blectin_segments),
                                                      source_rows, args.project_root, args.index_dir, 28, 28,
                                                      frozen_cds, frozen_junctions)
        write(args.out / "B_lectin_intron_metadata.tsv", b_observations, list(b_observations[0]))
        write(args.out / "blectin_focal_composition.tsv", b_metrics, list(b_metrics[0]))
        with (args.out / "B_lectin_current_genome_proteins.faa").open("w") as handle:
            for pid, sequence in sorted(b_proteins["B_lectin"].items()):
                handle.write(f">{pid}\n{sequence}\n")
    qc = {"eligible": len(metrics), "observed_proteins": len({r["protein_id"] for r in observations}),
          "families": dict(Counter(r["family"] for r in metrics)),
          "zero_intron": sum(r["intron_count"] == 0 for r in metrics)}
    (args.out / "QC.json").write_text(json.dumps(qc, indent=2, sort_keys=True) + "\n")
    print(json.dumps(qc, sort_keys=True))


if __name__ == "__main__":
    main()
