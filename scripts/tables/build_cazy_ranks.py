#!/usr/bin/env python3
"""Rank accepted PCWDE copies and scan-derived copies against fixed CAZy units."""

import argparse
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path

FOCAL = ('GH5_2', 'PL1_8', 'PL4_4')


def read(path):
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def write(path, fields, rows):
    with path.open('w', newline='') as handle:
        w = csv.writer(handle, delimiter='\t', lineterminator='\n')
        w.writerow(fields)
        w.writerows(rows)


def rank_full(count, native, reference_n):
    if count == 0:
        return '', '', ''
    above = sum(n > count for n in native.values())
    tied = sum(n == count for n in native.values())
    first = above + 1
    mid = first + tied / 2
    # The focal unit has been removed from native. A fixed N+1 denominator
    # makes the paired routes directly comparable, including an absent route.
    percentile = 100 * (reference_n + 1 - mid) / reference_n if reference_n else 100
    return first, mid, percentile


def rank(count, native, reference_n):
    first, mid, percentile = rank_full(count, native, reference_n)
    return (first, f'{mid:.2f}', f'{percentile:.2f}') if first else ('', '', '')


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--scan-counts', type=Path, required=True)
    p.add_argument('--member-table', type=Path, required=True,
                   help='frozen accepted member table, one row per gene copy')
    p.add_argument('--manifest', type=Path, required=True)
    p.add_argument('--member-crosswalk', type=Path, required=True,
                   help='sequence-backed accepted member to processed scan gene map')
    p.add_argument('--out', type=Path, required=True)
    args = p.parse_args()
    full_rows = []
    scan = args.scan_counts
    species = {r['species'] for r in read(args.manifest)}
    if len(species) != 64:
        raise ValueError('expected 64 Figure 1 species')
    by_species = defaultdict(dict)
    for r in read(scan / 'unit_counts.tsv'):
        by_species[r['species']][r['rank_unit']] = int(r['gene_count'])
    crosswalk_rows = read(args.member_crosswalk)
    crosswalk = {(r['species'], r['family'], r['accepted_protein_id']): r for r in crosswalk_rows}
    if len(crosswalk) != len(crosswalk_rows):
        raise ValueError('duplicate accepted member in crosswalk')
    gene_units = defaultdict(set)
    for r in read(scan / 'gene_units.tsv'):
        gene_units[(r['species'], r['gene_id'])].add(r['rank_unit'])
    accepted = Counter()
    accepted_genes = defaultdict(set)
    missing_hits, different_hits = [], []
    seen_member = set()
    for r in read(args.member_table):
        if r['family'] not in FOCAL:
            continue
        sp = r['species'].replace(' ', '_')
        if sp not in species:
            raise ValueError(f'accepted copy outside the Figure 1 roster: {sp}')
        family = r['family']
        member_key = (sp, family, r['protein_id'])
        if member_key in seen_member:
            raise ValueError(f'duplicate accepted member {member_key}')
        seen_member.add(member_key)
        accepted[(sp, family)] += 1
        item = crosswalk.get(member_key)
        if item is None:
            raise ValueError(f'accepted member missing from crosswalk: {member_key}')
        hit_genes = set(filter(None, item['scan_gene_ids'].split(',')))
        accepted_genes[sp].update(hit_genes)
        if item['dbcan_sub_hit'] == 'no':
            missing_hits.append((sp, family, r['protein_id']))
        elif not hit_genes:
            raise ValueError(f'hit member lacks mapped scan gene: {member_key}')
        elif not any(family in gene_units[(sp, gene)] for gene in hit_genes):
            different_hits.append((sp, family, r['protein_id'],
                                   ','.join(sorted({unit for gene in hit_genes for unit in gene_units[(sp, gene)]}))))
    native = defaultdict(lambda: defaultdict(set))
    for (sp, gene), units in gene_units.items():
        if gene in accepted_genes[sp]:
            continue
        for unit in units:
            if unit not in FOCAL:
                native[sp][unit].add(gene)
    rank_rows, native_rows, component_rows, component_native_rows, share_rows, primary_rows = [], [], [], [], [], []
    components = defaultdict(dict)
    for r in read(scan / 'component_counts.tsv'):
        components[r['species']][r['ecami_group']] = int(r['gene_count'])
    mapping = {r['ecami_group']: r['rank_unit'] for r in read(scan / 'group_mapping.tsv')}
    native_components = defaultdict(lambda: defaultdict(set))
    for r in read(scan / 'component_gene_units.tsv'):
        sp, gene, group = r['species'], r['gene_id'], r['ecami_group']
        if gene not in accepted_genes[sp] and mapping[group] not in FOCAL:
            native_components[sp][group].add(gene)
    for sp in sorted(species):
        background = {unit: len(genes) for unit, genes in native[sp].items()}
        for unit, count in sorted(background.items()):
            native_rows.append((sp, unit, count))
        n = len(background)
        family_only = sum(re.search(r'_\d+$', unit) is None for unit in background)
        share_rows.append((sp, n, family_only, f'{family_only/n:.4f}' if n else ''))
        for family in FOCAL:
            scan_count = by_species[sp].get(family, 0)
            accepted_count = accepted[(sp, family)]
            # Only the focal count changes. The native label roster and counts
            # are identical for both routes, even when one focal route is zero.
            for route, count in (('accepted', accepted_count), ('scan', scan_count)):
                first, mid, pct = rank(count, background, n)
                rank_rows.append((sp, family, route, count, first, mid, n + 1, pct,
                                  family_only, f'{family_only/n:.4f}' if n else ''))
        component_background = {group: len(genes) for group, genes in native_components[sp].items()}
        component_native_rows.extend((sp, group, count) for group, count in sorted(component_background.items()))
        for family in FOCAL:
            groups = sorted(g for g, unit in mapping.items() if unit == family)
            largest = max(groups, key=lambda g: (components[sp].get(g, 0), g)) if groups else ''
            for route, count in (('largest_group', components[sp].get(largest, 0)),
                                 ('pooled_scan_upper', by_species[sp].get(family, 0)),
                                 ('accepted_member', accepted[(sp, family)])):
                first, mid_full, pct_full = rank_full(count, component_background, len(component_background))
                mid, pct = (f'{mid_full:.2f}', f'{pct_full:.2f}') if first else ('', '')
                full_rows.append((sp, family, route, count, first, mid_full, len(component_background) + 1, pct_full))
                primary_rows.append((sp, family, route, largest if route == 'largest_group' else '',
                                     count, first, mid, len(component_background) + 1, pct))
        for group, count in components[sp].items():
            if mapping[group] in FOCAL or group == 'GH9_e27':
                reference = {g: c for g, c in component_background.items() if g != group}
                first, mid, pct = rank(count, reference, len(reference))
                component_rows.append((sp, group, mapping[group], count, first, mid,
                                       len(reference) + 1, pct))
        gh9_unit = mapping.get('GH9_e27', 'GH9')
        gh9_unit_count = by_species[sp].get(gh9_unit, 0)
        gh9_background = {unit: count for unit, count in background.items() if unit != gh9_unit}
        first, mid, pct = rank(gh9_unit_count, gh9_background, len(gh9_background))
        rank_rows.append((sp, gh9_unit, 'GH9_e27_native_unit', gh9_unit_count,
                          first, mid, len(gh9_background) + 1, pct, family_only,
                          f'{family_only/n:.4f}' if n else ''))
    args.out.mkdir(parents=True, exist_ok=False)
    write(args.out / 'primary_ecami_precision.tsv',
          ('species', 'unit', 'route', 'count', 'integer_rank', 'midrank',
           'rank_denominator', 'midrank_percentile_full_precision'), full_rows)
    if len(seen_member) != len(crosswalk):
        raise ValueError('crosswalk has extra or missing accepted members')
    write(args.out / 'paired_ranks.tsv', ('species','unit','route','count','integer_rank','midrank',
          'rank_denominator','midrank_percentile','family_level_native_units','family_level_unit_share'), rank_rows)
    write(args.out / 'native_counts.tsv', ('species','unit','gene_count'), native_rows)
    write(args.out / 'component_ranks.tsv', ('species','ecami_group','rank_unit','count',
          'integer_rank','midrank','rank_denominator','midrank_percentile'), component_rows)
    write(args.out / 'native_component_counts.tsv', ('species','ecami_group','gene_count'), component_native_rows)
    write(args.out / 'primary_ecami_ranks.tsv', ('species','unit','route','focal_ecami_group',
          'count','integer_rank','midrank','rank_denominator','midrank_percentile'), primary_rows)
    write(args.out / 'accepted_without_scan_hit.tsv', ('species','family','protein_id'), missing_hits)
    write(args.out / 'accepted_different_scan_unit.tsv', ('species','family','protein_id','scan_units'), different_hits)
    write(args.out / 'family_level_share.tsv', ('species','native_units','family_level_units','share'), share_rows)
    print(f'CAZY_RANKS_PASS genomes={len(species)} paired_rows={len(rank_rows)} '
          f'accepted_no_hit={len(missing_hits)} accepted_different_unit={len(different_hits)}')


if __name__ == '__main__':
    main()
