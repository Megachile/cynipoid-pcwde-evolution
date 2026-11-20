#!/usr/bin/env python3
"""Descriptive CAZy ranks and concrete native comparators by Figure 1 lineage."""

import argparse
import csv
import statistics
from collections import defaultdict
from pathlib import Path

FOCAL = ('GH5_2', 'PL1_8', 'PL4_4')
GROUPS = ('Diplolepididae', 'non-Cynipini Cynipidae', 'Cynipini', 'E. primus')


def read(path):
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def write(path, fields, rows):
    with path.open('w', newline='') as handle:
        w = csv.writer(handle, delimiter='\t', lineterminator='\n')
        w.writerow(fields)
        w.writerows(rows)


def fmt(value):
    return f'{value:.2f}'


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--lineages', type=Path, required=True)
    p.add_argument('--paired-ranks', type=Path, required=True)
    p.add_argument('--native-counts', type=Path, required=True)
    p.add_argument('--component-ranks', type=Path, required=True)
    p.add_argument('--primary-ranks', type=Path, required=True)
    p.add_argument('--native-components', type=Path, required=True)
    p.add_argument('--group-mapping', type=Path, required=True)
    p.add_argument('--primary-precision', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    lineage_rows = read(a.lineages)
    group_of = {r['species']: r['lineage_group'] for r in lineage_rows}
    if len(lineage_rows) != 64 or len(group_of) != 64:
        raise ValueError('expected 64 unique Figure 1 lineages')
    grouped = {g: sorted(s for s, v in group_of.items() if v == g) for g in GROUPS}
    if any(not v for v in grouped.values()):
        raise ValueError('empty required lineage group')
    paired = read(a.primary_ranks)
    precision = read(a.primary_precision)
    key = lambda r: (r['species'], r['unit'], r['route'])
    if len(precision) != 576 or len({key(r) for r in precision}) != len(precision):
        raise ValueError('expected 576 unique full-precision primary rows')
    if {key(r) for r in precision} != {key(r) for r in paired}:
        raise ValueError('missing precision join')
    legacy = {key(r): r for r in paired}
    for r in precision:
        old = legacy[key(r)]
        for field in ('count', 'integer_rank', 'rank_denominator'):
            if r[field] != old[field]:
                raise ValueError('precision/legacy rank mismatch')
        if r['integer_rank'] and f"{float(r['midrank_percentile_full_precision']):.2f}" != old['midrank_percentile']:
            raise ValueError('precision serialization differs')
    reader_group = {r['species']: (r['lineage_group'] if r['lineage_group'] != 'outside lineage summaries'
                    else 'other Figitidae' if r['source_family'] == 'Figitidae' else 'outside reader cohort')
                    for r in lineage_rows}
    full_summary = []
    for group in (*GROUPS, 'other Figitidae'):
        cohort = {s for s, g in reader_group.items() if g == group}
        for route in ('largest_group', 'pooled_scan_upper', 'accepted_member'):
            for family in FOCAL:
                present = [r for r in precision if r['species'] in cohort and
                           r['unit'] == family and r['route'] == route and r['integer_rank']]
                values = [float(r['midrank_percentile_full_precision']) for r in present]
                full_summary.append((group, family, route, len(cohort), len(present),
                                     statistics.median(values) if values else ''))
    native = {(r['species'], r['ecami_group']): int(r['gene_count']) for r in read(a.native_components)}
    labels = {r['ecami_group']: r['rank_unit'] for r in read(a.group_mapping)}
    native_units = {unit for _, unit in native}
    summary, closest, top, gh9 = [], [], [], []
    for group in GROUPS:
        genome_set = set(grouped[group])
        for route in ('largest_group', 'pooled_scan_upper', 'accepted_member'):
            for family in FOCAL:
                present = [r for r in paired if r['species'] in genome_set and
                           r['unit'] == family and r['route'] == route and int(r['count']) > 0]
                counts = [int(r['count']) for r in present]
                median_count = statistics.median(counts) if counts else None
                summary.append((group, family, route, len(genome_set), len(present),
                                fmt(median_count) if counts else '',
                                fmt(statistics.median(int(r['integer_rank']) for r in present)) if counts else '',
                                fmt(statistics.median(float(r['midrank_percentile']) for r in present)) if counts else ''))
                if counts:
                    cohort = [r['species'] for r in present]
                    medians = {unit: statistics.median(native.get((sp, unit), 0) for sp in cohort)
                               for unit in native_units
                               if any(native.get((sp, unit), 0) > 0 for sp in cohort)}
                    if medians:
                        distance = min(abs(m - median_count) for m in medians.values())
                        closest.extend((group, family, route, f'{unit} ({labels[unit]})', fmt(median_count),
                                        fmt(value), fmt(distance), len(cohort))
                                       for unit, value in sorted(medians.items())
                                       if abs(value - median_count) == distance)
        scores = []
        for unit in native_units:
            values = [native.get((sp, unit), 0) for sp in grouped[group]]
            if any(values):
                scores.append((unit, statistics.median(values), statistics.mean(values),
                               sum(v > 0 for v in values)))
        scores.sort(key=lambda row: (-row[1], -row[2], row[0]))
        top.extend((group, i, f'{unit} ({labels[unit]})', fmt(median), fmt(mean), n_present, len(genome_set))
                   for i, (unit, median, mean, n_present) in enumerate(scores[:10], 1))
        for route, source in [('GH9_e27_native_unit', read(a.paired_ranks)),
                              ('GH9_e27_component', read(a.component_ranks))]:
            selected = [r for r in source if r['species'] in genome_set and
                        (r['route'] == route if route == 'GH9_e27_native_unit'
                         else r['ecami_group'] == 'GH9_e27') and int(r['count']) > 0]
            gh9.append((group, route, len(genome_set), len(selected),
                        fmt(statistics.median(int(r['count']) for r in selected)) if selected else '',
                        fmt(statistics.median(int(r['integer_rank']) for r in selected)) if selected else '',
                        fmt(statistics.median(float(r['midrank_percentile']) for r in selected)) if selected else ''))
    a.out.mkdir(parents=True, exist_ok=False)
    write(a.out / 'lineage_precision.tsv', ('lineage_group', 'family', 'route', 'genomes_n',
          'genomes_with_family_present', 'median_midrank_percentile_full_precision'), full_summary)
    write(a.out / 'lineage_summary.tsv', ('lineage_group', 'focal_family', 'route',
          'group_genomes', 'positive_genomes', 'median_focal_count', 'median_integer_rank',
          'median_midrank_percentile'), summary)
    write(a.out / 'closest_native.tsv', ('lineage_group', 'focal_family', 'route',
          'native_unit', 'median_focal_count', 'median_native_count', 'absolute_difference',
          'cohort_genomes'), closest)
    write(a.out / 'top_native.tsv', ('lineage_group', 'order', 'native_unit',
          'median_count', 'mean_count', 'positive_genomes', 'group_genomes'), top)
    write(a.out / 'gh9_summary.tsv', ('lineage_group', 'route', 'group_genomes',
          'positive_genomes', 'median_count', 'median_integer_rank',
          'median_midrank_percentile'), gh9)
    print(f'CAZY_LINEAGE_PASS groups={len(GROUPS)} focal_rows={len(summary)} '
          f'closest={len(closest)} top={len(top)}')


if __name__ == '__main__':
    main()
