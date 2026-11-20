"""The PCWDE locus rule: group gene copies into loci from their ordered-gene positions.

Same-family copies on one scaffold join when at most MAX_LOCUS_GAP genes (any
gene, member or not) lie between them in the ordered gene list; a copy present
on several scaffold segments links its loci across segments (union-find). The
locus is anchored to the scaffold holding most of its copies and named
`{genome}_{scaffold}_pos{min}[_to_{max}]` from positions on that scaffold.

This is the joining logic of `find_loci_segment_aware` and `_make_locus` in
`scripts/analysis/synteny/01_build_synteny.py`, as run by the 2026-10-20
placement-conformance job 1344052, without that module's live chromosome
lookup. Equivalence on the 369 three-family loci is recorded in
scripts/deposit_inputs/README.md.

Environment: project Python 3 (standard library only).
"""
from collections import Counter

MAX_LOCUS_GAP = 5


def group_loci(occurrences, family):
    """occurrences: iterable of (genome, scaffold, rank, position, protein_id);
    family: {(genome, protein_id): family}. Returns a list of locus dicts."""
    parent = {key: key for key in family}

    def root(key):
        while parent[key] != key:
            parent[key] = parent[parent[key]]
            key = parent[key]
        return key

    def join(a, b):
        a, b = root(a), root(b)
        if a != b:
            parent[max(a, b)] = min(a, b)

    by_scaffold, seen = {}, {}
    for genome, scaffold, rank, position, protein in occurrences:
        key = (genome, protein)
        if key not in family:
            raise SystemExit(f"order row for a protein with no member row: {key}")
        by_scaffold.setdefault((genome, scaffold), []).append((rank, position, key))
        seen.setdefault(key, []).append((genome, scaffold, rank, position))
    missing = set(family) - set(seen)
    if missing:
        raise SystemExit(f"members without an ordered position: {sorted(missing)[:5]}")
    for entries in by_scaffold.values():
        entries.sort()
        by_family = {}
        for rank, position, key in entries:
            by_family.setdefault(family[key], []).append((rank, key))
        for ranked in by_family.values():
            for (left, a), (right, b) in zip(ranked, ranked[1:]):
                if right - left - 1 <= MAX_LOCUS_GAP:
                    join(a, b)
    components = {}
    for key in family:
        components.setdefault(root(key), []).append(key)
    loci = []
    for members in components.values():
        families = {family[key] for key in members}
        if len(families) != 1:
            raise SystemExit("cross-family locus component")
        counts = Counter((g, s) for key in members for g, s, _, _ in seen[key])
        genome, scaffold = sorted(counts, key=lambda k: (-counts[k], k))[0]
        positions = [p for key in members for g, s, _, p in seen[key] if (g, s) == (genome, scaffold)]
        locus_id = f"{genome}_{scaffold}_pos{min(positions)}" + (f"_to_{max(positions)}" if len(positions) > 1 else "")
        loci.append({"locus_id": locus_id, "family": families.pop(), "genome": genome, "scaffold": scaffold,
                     "protein_ids": sorted(key[1] for key in members)})
    return sorted(loci, key=lambda r: r["locus_id"])
