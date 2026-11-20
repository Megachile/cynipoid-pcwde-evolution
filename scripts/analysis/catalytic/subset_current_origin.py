#!/usr/bin/env python3
"""Subset pinned catalytic alignment calls to current full origin-domain IDs.

The source rows contain residues extracted from the pinned reference-augmented
alignments. Filtering rows preserves their alignment columns and residue calls;
only the per-site counts are recomputed here. No alignment or permutation runs.
"""
import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path

FAMILIES = ("GH5_2", "PL1_8", "PL4_4")
KNOWN_REMOVALS = {
    ("PL1_8", "cynipoid|Bassettia_pallida_1"),
    ("PL4_4", "cynipoid|Bassettia_pallida_5"),
    ("PL4_4", "cynipoid|Bassettia_pallida_8"),
}


def read_tsv(path):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader), reader.fieldnames


def write_tsv(path, rows, fields):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def cynipoid_ids(path):
    ids = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            sid = line[1:].split()[0]
            if sid.startswith("cynipoid|") or sid.split("|", 1)[0].endswith("_ENA"):
                ids.append(sid)
    if len(ids) != len(set(ids)):
        raise ValueError(f"duplicate current origin IDs: {path}")
    return set(ids)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    root = args.data_root / "repository_inputs"
    source = root / "tables/authority/catalytic_current"
    detail, detail_fields = read_tsv(source / "catalytic_site_conservation_by_sequence.tsv")
    previous, summary_fields = read_tsv(source / "catalytic_site_conservation_summary.tsv")
    removal_rows, _ = read_tsv(root / "deposit/base/sidecars/removed_cb_origin_bp_tips.tsv")
    recorded = {(r["tree"].split("origin_context_")[1].split("_full")[0], r["tip"])
                for r in removal_rows if "origin_context_" in r["tree"] and "_full.treefile" in r["tree"]}
    if not KNOWN_REMOVALS <= recorded:
        raise ValueError(f"removed IDs absent from tree removal record: {KNOWN_REMOVALS - recorded}")
    current = {family: cynipoid_ids(root / f"deposit/base/sequences/origin_context_{family}_full.faa")
               for family in FAMILIES}
    keys = defaultdict(set)
    for row in detail:
        if row["family"] not in current:
            raise ValueError(f"unexpected family: {row['family']}")
        keys[(row["alignment_set"], row["family"])].add(row["sequence_id"])
    removed = {(family, sid) for family in FAMILIES for sid in keys[("origin", family)] - current[family]}
    reverse = {(family, sid) for family in FAMILIES for sid in current[family] - keys[("origin", family)]}
    if removed != KNOWN_REMOVALS or reverse:
        raise ValueError(f"unexpected origin membership: removed={sorted(removed)} missing={sorted(reverse)}")
    for family in FAMILIES:
        if keys[("de_novo", family)] != keys[("origin", family)]:
            raise ValueError(f"de novo/origin ID mismatch in {family}")
        if not keys[("locus", family)] <= keys[("origin", family)]:
            raise ValueError(f"locus has IDs absent from origin in {family}")
    kept = [row for row in detail if row["sequence_id"] in current[row["family"]]]
    by_site = defaultdict(list)
    for row in kept:
        by_site[(row["alignment_set"], row["family"], row["site_label"])].append(row)
    rebuilt = []
    for old in previous:
        rows = by_site[(old["alignment_set"], old["family"], old["site_label"])]
        if len(rows) != len({r["sequence_id"] for r in rows}):
            raise ValueError(f"duplicate sequence at site: {old['family']} {old['site_label']}")
        counts = Counter(r["observed_aa"] for r in rows if r["site_status"] == "observed")
        expected = old["reference_residue"][0]
        total, observed, conserved = len(rows), sum(counts.values()), counts[expected]
        row = dict(old)
        row.update(n_total_sequences=str(total), n_site_observed=str(observed),
                   n_missing_or_gap=str(total-observed), n_expected_residue=str(conserved),
                   pct_expected_among_observed=f"{100*conserved/observed:.3f}" if observed else "NA",
                   pct_expected_among_all=f"{100*conserved/total:.3f}" if total else "NA",
                   observed_residue_counts=";".join(f"{aa}:{n}" for aa, n in sorted(counts.items())))
        rebuilt.append(row)
    args.out.mkdir(parents=True, exist_ok=False)
    write_tsv(args.out / "catalytic_site_conservation_by_sequence.tsv", kept, detail_fields)
    write_tsv(args.out / "catalytic_site_conservation_summary.tsv", rebuilt, summary_fields)
    print("CATALYTIC_SUBSET_PASS", {f: len(current[f]) for f in FAMILIES}, sorted(removed))


if __name__ == "__main__":
    main()
