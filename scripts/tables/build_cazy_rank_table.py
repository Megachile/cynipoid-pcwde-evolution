#!/usr/bin/env python3
"""Make the reader table from the paired, build-derived CAZy ranks."""
import argparse
import csv
import hashlib
import json
import statistics
from pathlib import Path


def read(path):
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--paired-ranks', type=Path, required=True)
    p.add_argument('--primary-ranks', type=Path, required=True)
    p.add_argument('--numbers', type=Path, required=True)
    p.add_argument('--lineage-groups', type=Path, required=True)
    p.add_argument('--lineage-precision', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    with a.numbers.open(newline='') as handle:
        labels = {r['table_id']: r['label'] for r in csv.DictReader(handle, delimiter='\t')}
    label = labels['cazy_rank']
    rows = [r for r in read(a.paired_ranks) if r['unit'] in {'GH5_2', 'PL1_8', 'PL4_4'}]
    by_key = {(r['species'], r['unit'], r['route']): r for r in rows}
    primary_rows = read(a.primary_ranks)
    primary = {(r['species'], r['unit'], r['route']): r for r in primary_rows}
    if len(primary_rows) != 576 or len(primary) != len(primary_rows):
        raise SystemExit('expected 576 unique eCAMI route rows')
    keys = sorted({(r['species'], r['unit']) for r in rows})
    if len(rows) != 384 or len(keys) != 192 or len(by_key) != len(rows):
        raise SystemExit('expected 192 unique paired genome-family comparisons')
    fields = ['species', 'family', 'largest_ecami_group', 'largest_group_count',
              'largest_group_rank', 'largest_group_percentile', 'pooled_scan_count',
              'pooled_scan_rank', 'pooled_scan_percentile', 'accepted_count',
              'accepted_rank', 'accepted_percentile', 'mixed_scan_rank',
              'mixed_accepted_rank', 'native_family_only_share']
    output = []
    count_changes = rank_changes = 0
    for species, unit in keys:
        accepted = by_key[species, unit, 'accepted']
        scan = by_key[species, unit, 'scan']
        if accepted['rank_denominator'] != scan['rank_denominator']:
            raise SystemExit(f'route denominator differs: {species} {unit}')
        count_changes += accepted['count'] != scan['count']
        rank_changes += accepted['integer_rank'] != scan['integer_rank']
        largest = primary[species, unit, 'largest_group']
        pooled = primary[species, unit, 'pooled_scan_upper']
        member = primary[species, unit, 'accepted_member']
        if len({largest['rank_denominator'], pooled['rank_denominator'], member['rank_denominator']}) != 1:
            raise SystemExit(f'eCAMI route denominator differs: {species} {unit}')
        if int(pooled['count']) != int(scan['count']) or int(member['count']) != int(accepted['count']):
            raise SystemExit(f'eCAMI and mixed route counts differ: {species} {unit}')
        output.append([species, unit, largest['focal_ecami_group'], largest['count'],
                       largest['integer_rank'], largest['midrank_percentile'],
                       pooled['count'], pooled['integer_rank'], pooled['midrank_percentile'],
                       member['count'], member['integer_rank'], member['midrank_percentile'],
                       scan['integer_rank'], accepted['integer_rank'], accepted['family_level_unit_share']])
    if (count_changes, rank_changes) != (32, 29):
        raise SystemExit(f'unexpected route sensitivity: {count_changes} counts, {rank_changes} ranks')
    lineage_rows = read(a.lineage_groups)
    if len(lineage_rows) != 64 or len({r['species'] for r in lineage_rows}) != 64:
        raise SystemExit('lineage authority is not the 64-genome panel')
    def group(row):
        if row['lineage_group'] != 'outside lineage summaries':
            return row['lineage_group']
        return 'other Figitidae' if row['source_family'] == 'Figitidae' else 'outside reader cohort'
    groups = {r['species']: group(r) for r in lineage_rows}
    if set(groups) != {species for species, _ in keys}:
        raise SystemExit('lineage authority and eCAMI rank species differ')
    group_order = ('Cynipini', 'non-Cynipini Cynipidae', 'Diplolepididae', 'E. primus',
                   'other Figitidae')
    if [sum(g == lineage for g in groups.values()) for lineage in group_order] != [12, 22, 3, 1, 10]:
        raise SystemExit('Table S6 cynipoid reader cohort changed')
    if sum(g == 'outside reader cohort' for g in groups.values()) != 16:
        raise SystemExit('Table S6 excluded reference cohort changed')
    reader_fields = ['lineage_group', 'family', 'genomes_n', 'genomes_with_family_present',
                     'median_rank', 'median_midrank_percentile', 'rank_range']
    reader = []
    for lineage in group_order:
        species = [s for s, g in groups.items() if g == lineage]
        for family in ('GH5_2', 'PL1_8', 'PL4_4'):
            found = [primary[s, family, 'largest_group'] for s in species]
            present = [r for r in found if r['integer_rank']]
            ranks = [int(r['integer_rank']) for r in present]
            percentiles = [float(r['midrank_percentile']) for r in present]
            if any(int(r['count']) <= 0 for r in present) or any(int(r['count']) > 0 for r in found if not r['integer_rank']):
                raise SystemExit(f'primary rank presence mismatch: {lineage} {family}')
            reader.append(dict(lineage_group=lineage, family=family, genomes_n=len(species),
                               genomes_with_family_present=len(present),
                               median_rank=f'{statistics.median(ranks):g}' if ranks else '',
                               median_midrank_percentile=f'{statistics.median(percentiles):.2f}' if percentiles else '',
                               rank_range=f'{min(ranks)}–{max(ranks)}' if ranks else ''))
    if len(reader) != 15 or sum(int(r['genomes_n']) for r in reader) != 144:
        raise SystemExit('Table S6 reader row coverage changed')
    precision_rows = read(a.lineage_precision)
    precision = {(r['lineage_group'], r['family'], r['route']): r for r in precision_rows}
    if len(precision_rows) != 45 or len(precision) != 45:
        raise SystemExit('expected 45 unique lineage precision rows')
    display_rows = []
    for r in reader:
        full = precision.get((r['lineage_group'], r['family'], 'largest_group'))
        if full is None or any(full[k] != str(r[k]) for k in ('genomes_n', 'genomes_with_family_present')):
            raise SystemExit('missing or inconsistent lineage precision join')
        display_rows.append({k: full[k] for k in ('lineage_group', 'family', 'median_midrank_percentile_full_precision')})
    a.out.mkdir(parents=True, exist_ok=False)
    with (a.out / 'Table_S14_precision.tsv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(display_rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader(); writer.writerows(display_rows)
    with (a.out / 'per_genome_cazy_ranks.tsv').open('w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(fields)
        writer.writerows(output)
    with (a.out / 'Table_S14.tsv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=reader_fields, delimiter='\t', lineterminator='\n')
        writer.writeheader(); writer.writerows(reader)
    (a.out / 'Table_S14_caption.txt').write_text(
        f'Table {label}. Lineage-level present gene-count ranks of three PCWDE families among native '
        'dbCAN-sub eCAMI groups in 48 cynipoid genomes. The primary focal count is its largest '
        'mapped eCAMI group. Medians, percentiles and rank ranges include only genomes with the family present; '
        'the genome and present-genome counts show the denominator. Percentiles use midranks. '
        'eCAMI groups are sequence clusters, not curated CAZy subfamilies. Full per-genome '
        'and secondary-route ranks, including the broader 64-genome scan, are in the data deposit. These are present-copy descriptions, not rates.\n')
    def digest(path):
        return hashlib.sha256(path.read_bytes()).hexdigest()
    (a.out / 'INPUT_OUTPUT_SHA256.json').write_text(json.dumps({
        'inputs': {name: digest(path) for name, path in
                   (('paired_ranks', a.paired_ranks), ('primary_ranks', a.primary_ranks),
                    ('lineage_groups', a.lineage_groups), ('lineage_precision', a.lineage_precision), ('table_numbers', a.numbers))},
        'outputs': {name: digest(a.out / name) for name in
                    ('Table_S14.tsv', 'per_genome_cazy_ranks.tsv', 'Table_S14_caption.txt', 'Table_S14_precision.tsv')},
    }, indent=2, sort_keys=True) + '\n')
    print(f'CAZY_TABLE_PASS comparisons={len(keys)} count_changes={count_changes} rank_changes={rank_changes}')


if __name__ == '__main__':
    main()
