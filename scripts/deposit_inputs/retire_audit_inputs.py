#!/usr/bin/env python3
"""Retire five superseded public inputs from a staged deposit."""
import csv
import sys
from pathlib import Path

RETIRED = {
    'figures/expression/metadata.tsv',
    *(f'tables/authority/table_s2/{family}_combined.faa' for family in ('GH5_2','PL1_8','PL4_4')),
    'deposit/base/data/selection/selection_membership_and_public_source.tsv',
    *(f'deposit/base/sequences/origin_context_{family}.faa' for family in ('GH5_2','PL1_8','PL4_4')),
}

def filter_manifest(path, removed):
    with path.open(newline='') as handle:
        rows = list(csv.reader(handle, delimiter='\t'))
    header, body = rows[0], rows[1:]
    if header[1] not in ('size_bytes','relative_path'):
        raise SystemExit(f'unexpected inventory schema: {path}')
    key = header.index('relative_path')
    found = {row[key] for row in body if row[key] in removed}
    if found != removed:
        raise SystemExit(f'inventory does not contain expected retired paths: {removed-found}')
    with path.open('w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerows([header, *(row for row in body if row[key] not in removed)])

def retire(inputs):
    for relative in sorted(RETIRED):
        path = inputs / relative
        if not path.is_file() or path.is_symlink():
            raise SystemExit(f'missing or unsafe retired input: {path}')
        path.unlink()
    filter_manifest(inputs/'figures/INPUT_MANIFEST.tsv', {'expression/metadata.tsv'})
    filter_manifest(inputs/'INPUT_MANIFEST.tsv', RETIRED)
    parent = inputs/'INPUT_MANIFEST.tsv'
    with parent.open(newline='') as handle:
        rows = list(csv.reader(handle, delimiter='\t'))
    nested = inputs/'figures/INPUT_MANIFEST.tsv'
    matches = [row for row in rows[1:] if row[0] == 'figures/INPUT_MANIFEST.tsv']
    if len(matches) != 1:
        raise SystemExit('expected one nested figure inventory pin')
    matches[0][1] = str(nested.stat().st_size)
    with parent.open('w', newline='') as handle:
        csv.writer(handle, delimiter='\t', lineterminator='\n').writerows(rows)
    print(f'RETIRED_INPUTS_PASS files={len(RETIRED)}')

if __name__ == '__main__':
    retire(Path(sys.argv[1]))
