#!/usr/bin/env python3
"""Summarize corrected PCWDE and 28-copy B-lectin GC3 against pinned comparators."""
from __future__ import annotations

import argparse
import bisect
import csv
import json
import statistics as st
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO
from scipy.stats import wilcoxon

FAMILIES = ("GH5_2", "PL1_8", "PL4_4")
ALIASES = {"Ceroptres masudai": "Ceroptres kovalevi",
           "Eschatocerus acacia": "Eschatocerus acaciae",
           "Isocolus centaureae": "Eubothrus centaureae"}


def read(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write(path: Path, rows: list[dict]):
    assert rows
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def norm(value: str):
    value = value.strip()
    for prefix in ("transcript:", "rna-", "cds-"):
        if value.startswith(prefix):
            value = value[len(prefix):]
    return value[:-2] if value.endswith(".1") else value


def value(raw: str):
    try:
        return float(raw)
    except (TypeError, ValueError):
        return None


def percentile(distribution: list[float], item: float):
    return 100 * (bisect.bisect_left(distribution, item) + bisect.bisect_right(distribution, item)) / (2 * len(distribution))


def paired_wilcoxon(values: list[float]):
    return float(wilcoxon(values).pvalue) if len(values) > 1 and any(x != 0 for x in values) else None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--host-dir", type=Path, required=True)
    parser.add_argument("--member-table", type=Path, required=True)
    parser.add_argument("--blectin-models", type=Path, required=True)
    parser.add_argument("--focal-composition", type=Path, required=True)
    parser.add_argument("--blectin-composition", type=Path, required=True)
    parser.add_argument("--bacterial-grid", type=Path, required=True)
    parser.add_argument("--blectin-bacterial-cds", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    members = [r for r in read(args.member_table) if r["family"] in FAMILIES]
    bmodels = read(args.blectin_models)
    assert len(members) == 629 and len(bmodels) == 28
    cohort = {r["protein_id"]: r for r in read(args.focal_composition)}
    bcohort = {r["protein_id"]: r for r in read(args.blectin_composition)}
    assert len(cohort) == 603 and len(bcohort) == 28
    targets = {}
    for member in members:
        pid = member["protein_id"]
        if pid in cohort:
            targets[(cohort[pid]["source_file"], norm(pid))] = (pid, member["family"], member["species"])
    for model in bmodels:
        pid = model["protein_id"]
        targets[(model["source_file"], norm(pid))] = (pid, "B_lectin", model["species"])
    assert len(targets) == 631
    backgrounds = {}
    joined = []
    for path in sorted(args.host_dir.glob("*.tsv")):
        if path.name in {"QC.tsv", "background_intron_details.tsv"}:
            continue
        source = path.stem
        rows = read(path)
        distributions = {metric: sorted(v for r in rows if r["cohort"] == "background"
                                        if (v := value(r[metric])) is not None)
                         for metric in ("gc3", "intron_density_per_kb", "median_intron_nt")}
        assert all(distributions.values()), source
        backgrounds[source] = distributions
        for row in rows:
            if row["cohort"] != "focal":
                continue
            keys = {key for candidate in row["candidate_ids"].split(";")
                    if (key := (source, norm(candidate))) in targets}
            if not keys:
                continue
            assert len(keys) == 1, (source, row["transcript_id"], keys)
            pid, family, species = targets[next(iter(keys))]
            metric = cohort[pid] if family != "B_lectin" else bcohort[pid]
            out = dict(protein_id=pid, family=family, species=species, source_file=source,
                       background_source=source, transcript_id=row["transcript_id"],
                       background_genes=len(distributions["gc3"]),
                       background_median_gc3=st.median(distributions["gc3"]))
            for name, distribution in distributions.items():
                number = value(row[name])
                if name == "gc3" and number is None:
                    number = value(metric["gc3_phase_corrected"])
                out[name] = number if number is not None else ""
                out[name + "_host_percentile"] = percentile(distribution, number) if number is not None else ""
            joined.append(out)
    by_id = {r["protein_id"]: r for r in joined}
    assert len(joined) == len(by_id) == 631 and set(by_id) == set(cohort) | set(bcohort)
    pcwde = sorted((r for r in joined if r["family"] in FAMILIES), key=lambda r: (r["family"], r["species"], r["protein_id"]))
    blectin = sorted((r for r in joined if r["family"] == "B_lectin"), key=lambda r: (r["species"], r["protein_id"]))
    write(args.out / "pcwde_gc3_per_copy.tsv", pcwde)
    write(args.out / "blectin_gc3_per_copy.tsv", blectin)

    bacteria = defaultdict(dict)
    for row in read(args.bacterial_grid):
        if row["group"] == "bacterial_relative" and row["cohort"] == "closest_relatives":
            bacteria[(row["family"], row["lineage"])][row["bacterial_accession"]] = float(row["gene_gc3"])
    stats, paired_rows = {}, []
    for family in FAMILIES:
        rows = [r for r in pcwde if r["family"] == family]
        by_species = defaultdict(list)
        for row in rows:
            by_species[row["species"]].append(row)
        pairs = []
        for species, genes in sorted(by_species.items()):
            host_pct = st.median(r["gc3_host_percentile"] for r in genes)
            host_dist = st.median(abs(r["gc3"] - r["background_median_gc3"]) * 100 for r in genes)
            bacterial_values = bacteria[(family, ALIASES.get(species, species))]
            assert bacterial_values, (family, species)
            host_mid = st.median(backgrounds[genes[0]["source_file"]]["gc3"])
            bacterial_dist = st.median(abs(v - host_mid) * 100 for v in bacterial_values.values())
            host_signed = st.median((r["gc3"] - r["background_median_gc3"]) * 100 for r in genes)
            bacterial_signed = st.median((v - host_mid) * 100 for v in bacterial_values.values())
            pairs.append((species, host_pct, host_dist, bacterial_dist, host_signed, bacterial_signed))
            paired_rows.append(dict(family=family, species=species, host_percentile=host_pct,
                                    host_distance_pp=host_dist, bacterial_distance_pp=bacterial_dist,
                                    bacterial_accessions=len(bacterial_values),
                                    copy_signed_offset_pp=host_signed, bacterial_signed_offset_pp=bacterial_signed))
        stats[family] = dict(copies=len(rows), species=len(by_species), paired_species=len(pairs),
                             median_species_host_percentile=st.median(x[1] for x in pairs),
                             species_above_host_median=sum(x[1] > 50 for x in pairs),
                             median_density_percentile=st.median(r["intron_density_per_kb_host_percentile"] for r in rows),
                             median_length_percentile=st.median(r["median_intron_nt_host_percentile"] for r in rows if r["median_intron_nt_host_percentile"] != ""),
                             median_copy_signed_offset_pp=st.median(x[4] for x in pairs),
                             median_bacterial_signed_offset_pp=st.median(x[5] for x in pairs),
                             host_distance_pp=st.median(x[2] for x in pairs),
                             bacterial_distance_pp=st.median(x[3] for x in pairs),
                             wilcoxon_bacteria_minus_host_distance_p=paired_wilcoxon([x[3] - x[2] for x in pairs]))
    assert [stats[f]["species"] for f in FAMILIES] == [37, 37, 36]
    write(args.out / "pcwde_species_pairs.tsv", paired_rows)
    (args.out / "pcwde_gc3_stats.json").write_text(json.dumps(stats, indent=2) + "\n")

    seq = str(next(SeqIO.parse(args.blectin_bacterial_cds, "fasta")).seq).upper()
    assert len(seq) == 3609 and set(seq) <= set("ACGT")
    bacterial_gc3 = sum(seq[i] in "GC" for i in range(2, len(seq), 3)) / (len(seq) // 3)
    model_by_id = {r["protein_id"]: r for r in bmodels}
    b_rows = []
    for row in blectin:
        model = model_by_id[row["protein_id"]]
        source = row["source_file"]
        mid = row["background_median_gc3"]
        b_rows.append(dict(**row, placement_id=model["placement_id"], reference_block=model["reference_block"],
                           architecture="chr6 full" if model["reference_block"] == "chr6" else
                           "chr3 reduced" if model["reference_block"] == "chr3" else "other/unplaced",
                           bacterial_accession="WP_067703811.1", bacterial_gc3=bacterial_gc3,
                           copy_signed_offset_pp=100 * (row["gc3"] - mid),
                           bacterial_signed_offset_pp=100 * (bacterial_gc3 - mid),
                           copy_absolute_distance_pp=100 * abs(row["gc3"] - mid),
                           bacterial_host_percentile=percentile(backgrounds[source]["gc3"], bacterial_gc3),
                           bacterial_absolute_distance_pp=100 * abs(bacterial_gc3 - mid)))
    assert len(b_rows) == 28 and len({r["species"] for r in b_rows}) == 16
    write(args.out / "blectin_gc3_comparison.tsv", b_rows)
    bstats = {}
    for label, selected in (("all", b_rows), ("chr6 full", [r for r in b_rows if r["architecture"] == "chr6 full"]),
                            ("chr3 reduced", [r for r in b_rows if r["architecture"] == "chr3 reduced"])):
        by_species = defaultdict(list)
        for row in selected:
            by_species[row["species"]].append(row)
        pairs = [(species, st.median(x["copy_absolute_distance_pp"] for x in genes),
                  st.median(x["bacterial_absolute_distance_pp"] for x in genes),
                  st.median(x["gc3_host_percentile"] for x in genes),
                  st.median(x["copy_signed_offset_pp"] for x in genes),
                  st.median(x["bacterial_signed_offset_pp"] for x in genes))
                 for species, genes in sorted(by_species.items())]
        bstats[label] = dict(copies=len(selected), species=len(pairs),
                             species_median_host_percentile=st.median(x[3] for x in pairs),
                             species_above_host_median=sum(x[3] > 50 for x in pairs),
                             median_copy_signed_offset_pp=st.median(x[4] for x in pairs),
                             median_bacterial_signed_offset_pp=st.median(x[5] for x in pairs),
                             median_copy_distance_pp=st.median(x[1] for x in pairs),
                             median_bacterial_distance_pp=st.median(x[2] for x in pairs),
                             wilcoxon_distance_p_two_sided=paired_wilcoxon([x[2] - x[1] for x in pairs]))
    (args.out / "blectin_gc3_stats.json").write_text(json.dumps(bstats, indent=2) + "\n")
    print(json.dumps({"pcwde": stats, "blectin": bstats}, sort_keys=True))


if __name__ == "__main__":
    main()
