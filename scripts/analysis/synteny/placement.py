#!/usr/bin/env python3
"""Single placement rule for accepted genome-derived targets (all four families).

Methods: "Using the Ragtag or JCVI B. kinseyi chromosome assignments, each target
hit in any family was assigned to one of the 10 B. kinseyi chromosomes, or marked
as unplaced if found on an unassigned scaffold."  Figure 2 caption: RagTag
assignment for scaffold-scale assemblies, JCVI synteny region match for
chromosome-scale assemblies; B. kinseyi itself uses its NCBI chromosomes.

Figure 2, Table S3 and the Results placed/total counts all read the table built
with `place()`; no consumer re-derives a chromosome label.
"""
from __future__ import annotations

UNPLACED = "unplaced"
ROUTES = ("RagTag", "JCVI", "NCBI-BK")

# B. kinseyi GCF_010883055.1 chromosomes: GenBank CM and RefSeq NC names.
BK_CHROMOSOMES = {
    **{f"CM0213{38 + i}.1": f"chr{i + 1}" for i in range(10)},
    **{f"NC_0466{57 + i}.1": f"chr{i + 1}" for i in range(10)},
}


def bk_chromosome(seqid: str) -> str | None:
    """BK chromosome for a BK sequence name or a RagTag object named after one."""
    if not seqid:
        return None
    name = seqid[: -len("_RagTag")] if seqid.endswith("_RagTag") else seqid
    return BK_CHROMOSOMES.get(name)


def place_ragtag(ragtag_object: str | None, unplaced_reason: str = "unassigned_ragtag_scaffold") -> tuple[str, str]:
    """Scaffold-scale route: the target's RagTag object is a BK chromosome object or not."""
    chrom = bk_chromosome(ragtag_object or "")
    if chrom:
        return chrom, "ragtag_bk_chromosome_object"
    return UNPLACED, unplaced_reason


def place_ncbi_bk(seqid: str) -> tuple[str, str]:
    """B. kinseyi NCBI models: the annotated chromosome itself."""
    chrom = bk_chromosome(seqid)
    if chrom:
        return chrom, "ncbi_bk_chromosome"
    return UNPLACED, "ncbi_bk_unplaced_scaffold"


def place_jcvi(contained_segment_chromosomes: set[str], partner_chromosomes: set[str]) -> tuple[str, str]:
    """Chromosome-scale route (coordinator rule of 2026-10-19).

    A target contained in simple-anchor segments of exactly one BK chromosome
    takes it; segments of different chromosomes -> unplaced.  Outside every
    segment, a native chromosome whose anchors pair with one BK chromosome
    takes it (whole-chromosome region match); otherwise unplaced.
    """
    if len(contained_segment_chromosomes) == 1:
        return next(iter(contained_segment_chromosomes)), "contained_jcvi_region"
    if len(contained_segment_chromosomes) > 1:
        return UNPLACED, "conflicting_jcvi_regions"
    if len(partner_chromosomes) == 1:
        return next(iter(partner_chromosomes)), "single_partner_chromosome"
    return UNPLACED, "outside_jcvi_region"


def place(route: str, **evidence) -> tuple[str, str]:
    if route == "RagTag":
        return place_ragtag(evidence.get("ragtag_object"), evidence.get("unplaced_reason", "unassigned_ragtag_scaffold"))
    if route == "NCBI-BK":
        return place_ncbi_bk(evidence["seqid"])
    if route == "JCVI":
        return place_jcvi(set(evidence["contained_segment_chromosomes"]), set(evidence["partner_chromosomes"]))
    raise ValueError(f"unknown placement route {route!r}")


def is_placed(label: str) -> bool:
    return label.startswith("chr")
