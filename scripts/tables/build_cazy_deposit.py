#!/usr/bin/env python3
"""Write a dense, traceable genome-by-rank-unit count deposit."""

import argparse
import csv
import re
from pathlib import Path


def read(path):
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--manifest', type=Path, required=True)
    p.add_argument('--unit-counts', type=Path, required=True)
    p.add_argument('--paired-ranks', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    species = {r['species'] for r in read(a.manifest)}
    if len(species) != 64:
        raise ValueError('expected 64 Figure 1 genomes')
    scan_rows = read(a.unit_counts)
    units = {r['rank_unit'] for r in scan_rows} | {'GH5_2', 'PL1_8', 'PL4_4'}
    counts = {(r['species'], r['rank_unit']): int(r['gene_count']) for r in scan_rows}
    if len(counts) != len(scan_rows):
        raise ValueError('duplicate genome by rank-unit scan count')
    focal = {(r['species'], r['unit']): int(r['count'])
             for r in read(a.paired_ranks)
             if r['route'] == 'accepted' and r['unit'] in {'GH5_2', 'PL1_8', 'PL4_4'}}
    if len(focal) != 64 * 3:
        raise ValueError('accepted focal count roster is incomplete')
    a.out.parent.mkdir(parents=True, exist_ok=True)
    with a.out.open('w', newline='') as handle:
        w = csv.writer(handle, delimiter='\t', lineterminator='\n')
        w.writerow(('species', 'rank_unit', 'unit_type', 'scan_gene_count',
                    'accepted_focal_gene_count'))
        for sp in sorted(species):
            for unit in sorted(units):
                w.writerow((sp, unit, 'named_subfamily' if re.search(r'_\d+$', unit)
                            else 'family_level', counts.get((sp, unit), 0),
                            focal.get((sp, unit), '')))
    print(f'CAZY_DEPOSIT_PASS genomes={len(species)} units={len(units)} '
          f'rows={len(species)*len(units)}')


if __name__ == '__main__':
    main()
