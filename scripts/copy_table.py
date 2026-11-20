"""Shared PCWDE copy table: the single source of per-genome copy counts.

The copy table (built by `make` from `repository_inputs/members/`, see scripts/members/) has one row per locus;
`protein_ids` lists the gene copies (protein models after fragment merging)
at that locus. Figure 1, Figure 2 and Table S2 all count copies from here.

Environment: project Python 3 (standard library only).
"""

import csv
from collections import Counter, defaultdict
from pathlib import Path

FAMILIES = ("GH5_2", "PL1_8", "PL4_4", "B-lectin")
CHROMOSOMES = tuple(f"chr{i}" for i in range(1, 11))
REQUIRED = {"species", "family", "locus_protein_id", "protein_ids", "bk_chromosome"}
# Database species names in the copy table -> names used in the paper.
PAPER_NAMES = {
    "Ceroptres masudai": "Ceroptres kovalevi",
    "Eschatocerus acacia": "Eschatocerus acaciae",
    "Isocolus centaureae": "Eubothrus centaureae",
}


def read_copies(path):
    """Return one record per gene copy (species as named in the paper)."""
    with Path(path).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not REQUIRED.issubset(reader.fieldnames or []):
            raise SystemExit(f"copy table lacks required columns: {path}")
        loci = list(reader)
    copies, seen = [], set()
    for locus in loci:
        chromosome = locus["bk_chromosome"]
        if locus["family"] not in FAMILIES:
            raise SystemExit(f"copy table has an unknown family: {locus['family']}")
        if chromosome != "unplaced" and chromosome not in CHROMOSOMES:
            raise SystemExit(f"copy table has an unknown chromosome: {chromosome}")
        for protein in (p.strip() for p in locus["protein_ids"].split(",")):
            if not protein:
                continue
            if protein in seen:
                raise SystemExit(f"copy table lists a protein twice: {protein}")
            seen.add(protein)
            species = PAPER_NAMES.get(locus["species"], locus["species"])
            copies.append({"species": species, "family": locus["family"],
                           "protein_id": protein, "locus_id": locus["locus_protein_id"],
                           "chromosome": chromosome, "scaffold": locus.get("scaffold", ""),
                           "placement_route": locus.get("placement_route", "")})
    return copies


def read_tblastn_presence(path):
    """Return the (species, family) pairs present only through a validated tBLASTn locus (presence rule b)."""
    with Path(path).open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    pairs = {(PAPER_NAMES.get(r["species"], r["species"]), r["family"])
             for r in rows if r["final_state"] == "admitted_validated_rescue"}
    unknown = {family for _, family in pairs} - set(FAMILIES)
    if unknown:
        raise SystemExit(f"tBLASTn presence table has unknown families: {sorted(unknown)}")
    return pairs


def copy_cells(copies):
    """Return {(species, family): Counter(total, unplaced, chr1..chr10)}."""
    cells = defaultdict(Counter)
    for copy in copies:
        cell = cells[copy["species"], copy["family"]]
        cell["total"] += 1
        cell[copy["chromosome"]] += 1
    return cells
