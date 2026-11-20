#!/usr/bin/env python3
"""Build the reader synteny workbook from the corrected block and flank products."""
import argparse
import csv
from collections import defaultdict
from pathlib import Path

from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill


def rows(path):
    with Path(path).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main():
    parser = argparse.ArgumentParser()
    for name in ("blocks", "states", "flanks", "roster", "counts", "out"):
        parser.add_argument("--" + name, required=True)
    args = parser.parse_args()
    blocks, states, flanks, roster = (rows(getattr(args, name)) for name in ("blocks", "states", "flanks", "roster"))
    count_rows = rows(args.counts)
    counts_by_block = {(r["family"], r["block_id"]): r for r in count_rows}
    if len(count_rows) != len(blocks) or len(counts_by_block) != len(blocks):
        raise SystemExit("Block count intermediate does not close against blocks")
    if not blocks or len(states) != len(blocks) * 49 or len(roster) != 49:
        raise SystemExit("Workbook inputs do not close at blocks x 49 genomes")
    by_state = {(r["family"], r["block_id"], r["genome"]): r for r in states}
    if len(by_state) != len(states):
        raise SystemExit("Duplicate block/genome state")
    by_flank = defaultdict(list)
    for row in flanks:
        by_flank[row["family"], row["block_id"], row["species"]].append(row)
    expected = {pid for block in blocks for pid in block["accepted_locus_ids"].split(";") if pid}
    if len(flanks) != len(expected) or {r["accepted_locus_id"] for r in flanks} != expected:
        raise SystemExit("Workbook flanks do not cover every eligible assigned copy")
    cross = [r for r in count_rows if r["cross_genome"] == "yes"]
    wb = Workbook()
    index = wb.active
    index.title = "Index"
    index.append(["Block", "Family", "Reference locus", "Detected genomes", "Occupied", "Empty", "Not detected", "Detailed sheet"])
    for number, block in enumerate(blocks, 1):
        key = block["family"], block["block_id"]
        counts = counts_by_block[key]
        detail_name = f"Block {number}" if counts["cross_genome"] == "yes" else ""
        index.append([f"Block {number}", key[0], block["accepted_locus_ids"].split(";")[0],
                      int(counts["detected_genomes"]), int(counts["occupied_genomes"]), int(counts["empty_genomes"]), int(counts["not_detected_genomes"]), detail_name])
        if not detail_name:
            continue
        sheet = wb.create_sheet(detail_name)
        sheet.append([f"{detail_name} — {key[0]} conserved locus"])
        sheet.append([f"Detected in {counts['detected_genomes']}/49 genomes: {counts['occupied_genomes']} occupied, {counts['empty_genomes']} empty, {counts['not_detected_genomes']} not detected."])
        sheet.append([])
        flanks_header = [f"U{i}" for i in range(12, 0, -1)] + [f"D{i}" for i in range(1, 13)]
        sheet.append(["Species", "Genome/accession", "Chromosome/scaffold", *flanks_header[:12], "TARGET LOCUS", *flanks_header[12:]])
        for genome in roster:
            state = by_state[key + (genome["builder_genome_key"],)]
            copies = sorted(by_flank[key + (genome["species"].replace(" ", "_"),)], key=lambda r: r["accepted_locus_id"])
            # Only block members have a focal flank row. Empty and not-detected
            # cells have none; neither do cells occupied by the post-detection
            # region scan, which list their target hits without flanks.
            member = state["search_route"] == "accepted_block_member"
            if state["status"] == "PRESENT" and member and not copies:
                raise SystemExit(f"Occupied block lacks a flank row: {key} {genome['builder_genome_key']}")
            if state["status"] == "PRESENT" and not member and not state["target_unit_ids"]:
                raise SystemExit(f"Region-scan occupancy lacks a target: {key} {genome['builder_genome_key']}")
            if (state["status"] != "PRESENT" or not member) and copies:
                raise SystemExit(f"Nonmember block cell has a flank row: {key} {genome['builder_genome_key']}")
            for copy in copies or [None]:
                if copy:
                    target = copy["accepted_locus_id"]
                elif state["status"] == "PRESENT":
                    target = state["target_unit_ids"].replace(";", ", ") + " (in search region)"
                else:
                    target = "empty" if state["status"] == "EMPTY" else "not detected"
                sheet.append([genome["species"], genome["accession_or_project"], copy["scaffold"] if copy else state["scaffold"],
                              *([copy[f"U{i}"] for i in range(12, 0, -1)] if copy else [""] * 12),
                              target,
                              *([copy[f"D{i}"] for i in range(1, 13)] if copy else [""] * 12)])
        sheet.freeze_panes = "D5"
        sheet.auto_filter.ref = f"A4:AB{sheet.max_row}"
    for sheet in wb:
        for cell in sheet[1]:
            cell.font = Font(bold=True, color="FFFFFF")
            cell.fill = PatternFill("solid", fgColor="234F70")
    index.freeze_panes = "A2"
    index.auto_filter.ref = f"A1:H{index.max_row}"
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    wb.save(out)
    print(f"PASS workbook: {len(blocks)} indexed blocks, {len(cross)} detailed blocks, 49 genomes")


if __name__ == "__main__":
    main()
