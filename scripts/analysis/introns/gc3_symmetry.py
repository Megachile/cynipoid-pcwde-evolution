#!/usr/bin/env python3
"""Report-only phase-aware host GC3 sensitivity for retained build-8 outputs."""
import csv
import statistics as st
import sys
from collections import defaultdict
from pathlib import Path

from pyfaidx import Fasta


def rows(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def percentile(distribution, item):
    import bisect
    return 100 * (bisect.bisect_left(distribution, item) + bisect.bisect_right(distribution, item)) / (2 * len(distribution))


build = Path(sys.argv[1])
data = Path(sys.argv[2])
out = Path(sys.argv[3])
host_dir = build / "introns/host_background"
comp_dir = build / "introns/gc3"
aliases = defaultdict(dict)
for row in rows(Path("repo/config/intron_seqid_aliases.tsv")):
    aliases[row["source_file"]][row["gff_seqid"]] = row["assembly_seqid"]
base_gc3 = {}
sym_gc3 = {}
missing_total = recoverable_total = focal_fallback = 0
details = []
for host_file in sorted(host_dir.glob("*.tsv")):
    if host_file.name in {"QC.tsv", "background_intron_details.tsv"}:
        continue
    source = host_file.stem
    host = rows(host_file)
    background = [r for r in host if r["cohort"] == "background"]
    missing = {r["transcript_id"]: r for r in background if not r["gc3"]}
    missing_total += len(missing)
    focal_fallback += sum(r["cohort"] == "focal" and not r["gc3"] for r in host)
    base = sorted(float(r["gc3"]) for r in background if r["gc3"])
    recovered = {}
    if missing:
        gff = data / f"repository_inputs/analysis/introns/sources/{source}/gff.gff"
        segments = defaultdict(list)
        with gff.open(errors="replace") as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) != 9 or f[2] != "CDS":
                    continue
                attrs = dict(x.split("=", 1) for x in f[8].split(";") if "=" in x)
                for tid in attrs.get("Parent", "").split(","):
                    if tid in missing:
                        segments[tid].append((f[0], int(f[3]), int(f[4]), f[6], f[7]))
        assembly = data / f"repository_inputs/analysis/introns/sources/{source}/assembly.fasta"
        index = build / f"introns/indexes/{source}.fai"
        genome = Fasta(str(assembly), indexname=str(index), rebuild=False,
                       as_raw=True, sequence_always_upper=True)
        complement = str.maketrans("ACGTN", "TGCAN")
        for tid, pieces in segments.items():
            if len({(p[0], p[3]) for p in pieces}) != 1:
                continue
            pieces.sort(key=lambda p: p[1])
            scaffold = aliases[source].get(pieces[0][0], pieces[0][0])
            fragments = [str(genome[scaffold][left-1:right]) for _, left, right, _, _ in pieces]
            sequence = "".join(fragments)
            if pieces[0][3] == "-":
                sequence = sequence.translate(complement)[::-1]
            first = pieces[-1] if pieces[0][3] == "-" else pieces[0]
            phase = int(first[4]) if first[4] in {"0", "1", "2"} else 0
            coding = sequence[phase:]
            codons = [coding[i:i+3] for i in range(0, len(coding)-2, 3)
                      if set(coding[i:i+3]) <= set("ACGT")]
            if codons:
                recovered[tid] = sum(c[2] in "GC" for c in codons) / len(codons)
        genome.close()
    recoverable_total += len(recovered)
    amended = sorted(base + list(recovered.values()))
    base_gc3[source] = base
    sym_gc3[source] = amended
    details.append((source, len(background), len(missing), len(recovered), len(base), len(amended)))


def summarize(path):
    old, new = defaultdict(list), defaultdict(list)
    for row in rows(path):
        source = row["source_file"]
        score = float(row["gc3"])
        key = (row["family"], row["species"])
        old[key].append(percentile(base_gc3[source], score))
        new[key].append(percentile(sym_gc3[source], score))
    by_family = defaultdict(list)
    for key in old:
        by_family[key[0]].append((key[1], st.median(old[key]), st.median(new[key])))
    return {family: (len(items), st.median(x[1] for x in items), st.median(x[2] for x in items),
                     sum(x[1] > 50 for x in items), sum(x[2] > 50 for x in items))
            for family, items in by_family.items()}


pc = summarize(comp_dir / "pcwde_gc3_per_copy.tsv")
bl = summarize(comp_dir / "blectin_gc3_per_copy.tsv")
with out.open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["family", "species", "primary_median_percentile", "symmetric_median_percentile",
                     "primary_species_above_host", "symmetric_species_above_host",
                     "host_missing_gc3", "host_recovered_phase_aware", "focal_fallback"])
    for family, (species, primary, symmetric, above_old, above_new) in sorted((pc | bl).items()):
        writer.writerow([family, species, primary, symmetric, above_old, above_new,
                         missing_total, recoverable_total, focal_fallback])
with out.with_name("GC3_SYMMETRY_SOURCES.tsv").open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["source", "host_representatives", "missing_gc3", "phase_aware_recovered",
                     "primary_eligible", "symmetric_eligible"])
    writer.writerows(details)
print(f"GC3_SYMMETRY_PASS missing={missing_total} recovered={recoverable_total} focal_fallback={focal_fallback}")
