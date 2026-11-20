#!/usr/bin/env python3
"""Drop exact CB-origin BP records from tree-feeding FASTA and TSV inputs.

Lines of retained records are copied byte-for-byte; no realignment, column
trimming or reformatting.  Each file must lose exactly the expected number of
records.  Originals are saved under outputs/prune/before/.
Usage: drop_rows.py DATA_ROOT OUT_EVIDENCE.tsv
"""
from __future__ import annotations

import csv
import shutil
import sys
from pathlib import Path

PL18 = {("PL1_8", "cynipoid|Bassettia_pallida_1")}
PL44 = {("PL4_4", "cynipoid|Bassettia_pallida_5"), ("PL4_4", "cynipoid|Bassettia_pallida_8")}
FAMILY_TIPS = PL18 | PL44
BLECTIN = "focal|transcript|Bassettia_pallida|BP_TRINITY_DN2383_c0_g1__Blectin_d1"

B = "repository_inputs/deposit/base/"
F = "repository_inputs/figures/"
FASTA = {
    B + "sequences/complete_locus_PL1_8.faa": {"cynipoid|Bassettia_pallida_1"},
    B + "alignments/complete_locus_PL1_8.aln.faa": {"cynipoid|Bassettia_pallida_1"},
    B + "sequences/complete_locus_PL4_4.faa": {"cynipoid|Bassettia_pallida_5", "cynipoid|Bassettia_pallida_8"},
    B + "alignments/complete_locus_PL4_4.aln.faa": {"cynipoid|Bassettia_pallida_5", "cynipoid|Bassettia_pallida_8"},
    B + "sequences/origin_context_PL1_8_full.faa": {"cynipoid|Bassettia_pallida_1"},
    B + "alignments/origin_context_PL1_8_full.aln.faa": {"cynipoid|Bassettia_pallida_1"},
    B + "sequences/origin_context_PL4_4_full.faa": {"cynipoid|Bassettia_pallida_5", "cynipoid|Bassettia_pallida_8"},
    B + "alignments/origin_context_PL4_4_full.aln.faa": {"cynipoid|Bassettia_pallida_5", "cynipoid|Bassettia_pallida_8"},
    B + "sequences/B_lectin_domains.faa": {BLECTIN},
    B + "alignments/B_lectin_domains.aln.faa": {BLECTIN},
    F + "figure_s4/B_lectin_domains.faa": {BLECTIN},
    F + "figure_s4/B_lectin_domains.aln.faa": {BLECTIN},
}
# path -> (key columns, set of key tuples)
TSV = {
    B + "sidecars/complete_locus_PL1_8_tip_labels.tsv": (("raw_tip_id",), {(t,) for _, t in PL18}),
    B + "sidecars/complete_locus_PL4_4_tip_labels.tsv": (("raw_tip_id",), {(t,) for _, t in PL44}),
    B + "sidecars/origin_context_full_tip_mapping.tsv": (("family", "tree_tip_id"), FAMILY_TIPS),
    B + "sidecars/B_lectin_architecture.tsv": (("tip_id",), {(BLECTIN,)}),
    B + "sidecars/B_lectin_domain_manifest.tsv": (("tip_id",), {(BLECTIN,)}),
    F + "figure_s4/architecture.tsv": (("tip_id",), {(BLECTIN,)}),
    F + "figure_s4/domain_manifest.tsv": (("tip_id",), {(BLECTIN,)}),
    F + "figure6/all_family_rna_evidence_overlay.tsv": (("family", "tip"), FAMILY_TIPS),
    F + "figure6/all_family_rna_evidence_audit.tsv": (("family", "tip"), FAMILY_TIPS),
    F + "figure6/approved_roster_projection.tsv": (("family", "current_full_tip"), FAMILY_TIPS),
    F + "figure6/label_overrides.tsv": (("family", "raw_tree_tip"), FAMILY_TIPS),
    F + "figure6/raw_tip_to_normalized_display_label.tsv": (("family", "raw_tree_tip"), FAMILY_TIPS),
    F + "figure6/tree_tip_identifier_crosswalk.tsv": (("family", "tree_tip"), FAMILY_TIPS),
    F + "figure6/accepted_old_vs_current_row_audit.tsv": (("family", "current_tree_tip"), FAMILY_TIPS),
    F + "figure_s3/main_figure6_rna_evidence_overlay.tsv": (("family", "tip"), FAMILY_TIPS),
    F + "figure7/interface.tsv": (("one_best_tip_id",), {("L037",)}),
    F + "figure7/domain_lane_alignment.tsv": (("locus_id",), {("L037",)}),
}
# Rows present only for some keys (e.g. reduced roster) — expected drop counts.
EXPECTED_TSV_DROPS = {
    F + "figure6/approved_roster_projection.tsv": 1,
    F + "figure6/accepted_old_vs_current_row_audit.tsv": 1,
}


def save(root: Path, rel: str, before: Path) -> None:
    target = before / rel
    if not target.exists():
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(root / rel, target)


def drop_fasta(path: Path, ids: set[str]) -> tuple[int, int, list[str]]:
    lines = path.read_text(encoding="utf-8").splitlines(keepends=True)
    out, dropped, total, skipping = [], [], 0, False
    for line in lines:
        if line.startswith(">"):
            total += 1
            rid = line[1:].strip().split()[0]
            skipping = rid in ids
            if skipping:
                dropped.append(rid)
        if not skipping:
            out.append(line)
    if sorted(dropped) != sorted(ids):
        raise SystemExit(f"FASTA drop mismatch {path}: {dropped}")
    path.write_text("".join(out), encoding="utf-8")
    return total, total - len(dropped), dropped


def drop_tsv(path: Path, cols: tuple[str, ...], keys: set[tuple]) -> tuple[int, int, list[str]]:
    lines = path.read_text(encoding="utf-8").splitlines(keepends=True)
    header = lines[0].rstrip("\r\n").split("\t")
    idx = [header.index(c) for c in cols]
    out, dropped = [lines[0]], []
    for line in lines[1:]:
        cells = line.rstrip("\r\n").split("\t")
        key = tuple(cells[i] for i in idx)
        if key in keys:
            dropped.append("|".join(key))
        else:
            out.append(line)
    path.write_text("".join(out), encoding="utf-8")
    return len(lines) - 1, len(out) - 1, dropped


def main() -> None:
    root, evidence = Path(sys.argv[1]), Path(sys.argv[2])
    before = evidence.parent / "before"
    rows = []
    for rel, ids in FASTA.items():
        save(root, rel, before)
        n0, n1, dropped = drop_fasta(root / rel, ids)
        rows.append((rel, "fasta", n0, n1, ";".join(sorted(dropped))))
    for rel, (cols, keys) in TSV.items():
        save(root, rel, before)
        n0, n1, dropped = drop_tsv(root / rel, cols, keys)
        expected = EXPECTED_TSV_DROPS.get(rel, len(keys))
        if len(dropped) != expected or len(set(dropped)) != len(dropped):
            raise SystemExit(f"TSV drop mismatch {rel}: {dropped} expected {expected}")
        rows.append((rel, "tsv:" + ",".join(cols), n0, n1, ";".join(sorted(dropped))))
    with evidence.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["file", "kind", "records_before", "records_after", "dropped_keys"])
        w.writerows(rows)
    for r in rows:
        print(*r, sep="\t")


if __name__ == "__main__":
    main()
