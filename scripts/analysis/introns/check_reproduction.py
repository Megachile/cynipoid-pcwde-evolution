#!/usr/bin/env python3
"""Check fixed cohort and site controls before rendering intron results."""
from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from pathlib import Path


def rows(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--accounting", type=Path, required=True)
    parser.add_argument("--observations", type=Path, required=True)
    parser.add_argument("--table-s5", type=Path, required=True)
    parser.add_argument("--host-qc", type=Path, required=True)
    parser.add_argument("--gc3-stats", type=Path, required=True)
    parser.add_argument("--focal-composition", type=Path, required=True)
    args = parser.parse_args()
    accounting = rows(args.accounting)
    assert len(accounting) == 629 and len({r["protein_id"] for r in accounting}) == 629
    eligible = [r for r in accounting if r["status"] == "single_scaffold"]
    assert len(eligible) == 603 and len(accounting) - len(eligible) == 26
    observations = [r for r in rows(args.observations) if r["family"] in {"GH5_2", "PL1_8", "PL4_4"} and int(r["intron_number"]) > 0]
    assert len(observations) == 939
    assert sum(float(r["relative_cds_position"]) <= .25 for r in observations) == 799
    assert math.isclose(statistics.median(float(r["relative_cds_position"]) for r in observations), .0547, abs_tol=1e-6)
    focal = [r for r in rows(args.focal_composition) if r["family"] in {"GH5_2", "PL1_8", "PL4_4"}]
    assert len(focal) == 603
    assert math.isclose(statistics.median(float(r["intron_density_per_kb"]) for r in focal), 1.0416666666666667, abs_tol=1e-6)
    s5 = rows(args.table_s5)
    part_a = {r["family_id"]: r for r in s5 if r["part"] == "A"}
    for family, expected in {"GH5_2": (176, 164, 229), "PL1_8": (212, 201, 303),
                             "PL4_4": (215, 196, 407), "B_lectin": (28, 6, 35)}.items():
        row = part_a[family]
        assert tuple(int(row[k]) for k in ("eligible_proteins", "proteins_with_introns", "total_introns")) == expected, family
    sites = {(r["family_id"], int(r["reference_position_aa"])): r for r in s5 if r["part"] == "B"}
    for key, (count, phase) in {("GH5_2", 18): (154, "1 (151/154)"),
                                 ("PL1_8", 9): (88, "0 (86/88)"),
                                 ("PL4_4", 37): (32, "1 (27/32)"),
                                 ("PL4_4", 83): (85, "1 (42/85)")}.items():
        row = sites[key]
        assert int(row["observation_count"]) == count and row["phase"] == phase, key
    host = rows(args.host_qc)
    assert len(host) == 38 and sum(int(r["target_ids"]) for r in host) == 631
    assert sum(int(r["focal_rows"]) for r in host) == 631
    stats = json.loads(args.gc3_stats.read_text())
    assert [stats[f]["species"] for f in ("GH5_2", "PL1_8", "PL4_4")] == [37, 37, 36]
    assert [stats[f]["copies"] for f in ("GH5_2", "PL1_8", "PL4_4")] == [176, 212, 215]
    # The BK GNOMON 31 host background replaces the candidate's BK Helixer
    # background. Its GH5_2 rank moves 65.36 -> 47.63 and the family count
    # above the host median moves 32 -> 31. The 649,754-gene background also
    # includes the diagnosed sequence-ID alias repair.
    for family, percentile, above in (
        ("GH5_2", 61.138118387189024, 31),
        ("PL1_8", 25.6538585728124, 1),
        ("PL4_4", 35.47036132895277, 2),
    ):
        assert math.isclose(stats[family]["median_species_host_percentile"], percentile, abs_tol=1e-6), family
        assert stats[family]["species_above_host_median"] == above, family
    blectin_path = args.gc3_stats.with_name("blectin_gc3_stats.json")
    blectin = json.loads(blectin_path.read_text())["all"]
    assert (blectin["copies"], blectin["species"], blectin["species_above_host_median"]) == (28, 16, 16)
    # The 90.6287 -> 88.7116 shift is solely the decided BK background swap.
    assert math.isclose(blectin["species_median_host_percentile"], 88.7115792284225, abs_tol=1e-6)
    print("INTRON_REPRODUCTION_PASS 629/603/26, 939/799, Table S5 sites, 38/631 hosts, 37/37/36 GC3 species")


if __name__ == "__main__":
    main()
