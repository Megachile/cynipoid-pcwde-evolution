"""Shared S4 counting units and comprehensive, post-detection copy coverage.

The legacy unplaced PL1_8 rescue is a printed S4 unit, with no detected
region on its scaffold. It enters the denominator once, with no block hit.
Model-less copies enter membership only through retained in-region hits;
shared chromosome or an occupied companion alone never implies containment.
"""
from collections import Counter, defaultdict


def comprehensive_coverage(copies, members, states, species_by_genome):
    units = {r['locus_id']: r for r in copies}
    if len(units) != len(copies):
        raise ValueError('duplicate S4 counting unit')
    expected = {r['accepted_locus_id'] for r in members}
    rescue = 'NODE_38372_length_2218_cov_4.820177:1323-346'
    if set(units) - expected != {rescue} or expected - set(units):
        raise ValueError('S4/synteny roster mismatch beyond legacy rescue')
    edges = defaultdict(set)
    reasons = defaultdict(set)
    def add(unit, family, block, species, reason):
        if unit not in units:
            raise ValueError(f'unknown accepted copy in block: {unit}')
        row = units[unit]
        if row['family_id'] != family or row['species_id'] != species:
            raise ValueError(f'copy family/species mismatch: {unit}')
        edges[(family, block)].add((unit, species))
        reasons[(unit, block)].add(reason)
    public_family = lambda f: 'B-lectin' if f == 'beta_lectin/GH30_lineage' else f
    for row in members:
        if row['assignment_status'] == 'assigned':
            unit = row['accepted_locus_id']
            add(unit, public_family(row['family']), row['block_id'], units[unit]['species_id'], 'block_defining_member')
    for row in states:
        targets = [x for x in row['target_unit_ids'].split(';') if x]
        if targets and row['status'] != 'PRESENT':
            raise ValueError('target membership in non-PRESENT cell')
        for unit in targets:
            add(unit, public_family(row['family']), row['block_id'], species_by_genome[row['genome']], 'retained_target_hit')
    multi = {key for key, hits in edges.items() if len({species for _, species in hits}) >= 2}
    by_unit = defaultdict(set)
    for key, hits in edges.items():
        for unit, _ in hits:
            by_unit[unit].add(key)
    evidence = []
    for unit, row in sorted(units.items()):
        keys = by_unit[unit]
        evidence.append(dict(locus_id=unit, family_id=row['family_id'], species_id=row['species_id'],
            counting_unit='legacy_raw_PL1_8_rescue' if unit == rescue else 'accepted_biological_copy',
            blocks=';'.join(sorted(block for _, block in keys)),
            multi_genome_blocks=';'.join(sorted(block for family, block in keys if (family, block) in multi)),
            included=int(bool(keys & multi)),
            evidence=';'.join(sorted({reason for _, block in keys for reason in reasons[(unit, block)]})) or 'no_retained_in_region_membership'))
    totals = Counter(r['family_id'] for r in evidence)
    covered = Counter(r['family_id'] for r in evidence if r['included'])
    if any(covered[f] > n for f, n in totals.items()):
        raise ValueError('coverage exceeds shared S4 total')
    return evidence, totals, covered, multi
