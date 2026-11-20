#!/usr/bin/env python3
"""Validate the frozen intron/GC3 producer-consumer boundary before compute."""
import argparse
import csv
import hashlib
import json
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'release'))
from data_root import validate_reader_inputs


def rows(path):
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def require_columns(path, names):
    with path.open() as handle:
        found = set(handle.readline().rstrip('\n').split('\t'))
    assert set(names) <= found, (path, sorted(set(names) - found))


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--data-root', type=Path, required=True)
    p.add_argument('--sources', type=Path, required=True)
    p.add_argument('--members', type=Path, required=True)
    p.add_argument('--blectin-models', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--validate-inputs-only', action='store_true')
    p.add_argument('--frozen', type=Path)
    a = p.parse_args()
    assert not a.out.exists(), f'output destination already exists: {a.out}'
    repo = Path(__file__).resolve().parents[3]
    validate_reader_inputs(repo, a.data_root)
    sources = rows(a.sources)
    assert len(sources) == len({r['source_file'] for r in sources}) == 38
    if a.frozen:
        producer = json.loads((a.frozen / 'PRODUCER.json').read_text())
        from manifest_bridge import validate_bridge
        validate_bridge(repo, a.sources, producer['source_manifest_sha256'])
        assert (producer['sources'], producer['focal_copies'], producer['cohort'],
                producer['excluded'], producer['blectin']) == (38, 631, 603, 26, 28)
        pinned = rows(a.frozen / 'FROZEN_MANIFEST.tsv')
        assert len({r['path'] for r in pinned}) == len(pinned)
        for r in pinned:
            path = (a.frozen / r['path']).resolve()
            assert path.is_relative_to(a.frozen.resolve()) and path.is_file()
            assert path.stat().st_size == int(r['bytes']) and digest(path) == r['digest'], path
    else:
        projected = (a.data_root / 'PROJECTED_INPUTS.tsv').exists()
        for r in sources:
            for kind in ('gff', 'assembly', 'proteome'):
                path = (a.data_root / r[kind + '_snapshot_path']).resolve()
                assert path.is_relative_to(a.data_root.resolve()), path
                assert path.is_file(), path
                # Projected GFF headers differ from their raw source pins.
                if not projected or kind != 'gff':
                    assert digest(path) == r[kind + '_sha256'], path
    if a.validate_inputs_only:
        print('INTRON_INPUT_FORM_PASS sources=38 frozen=' + str(bool(a.frozen)))
        return
    base = a.data_root / 'repository_inputs/analysis/introns/frozen' if a.frozen else a.data_root / 'repository_inputs/analysis/introns/sources'
    grid = base / 'bacterial_grid.tsv'
    cds = base / 'blectin_bacterial_cds.fna'
    assert digest(grid) == '983a790cd4c5a1cfa3e9cbd9fd928d80fee9d660114ab4c8453eb82c5f56ce56'
    assert digest(cds) == '09bc2bed98e4dc6af53aee753507dbcaea6b18dea16aa765820277fa36381d37'
    require_columns(grid, ('group', 'cohort', 'family', 'lineage', 'bacterial_accession', 'gene_gc3'))
    require_columns(a.members, ('protein_id', 'family', 'species', 'genome'))
    require_columns(a.blectin_models, ('protein_id', 'species', 'source_file'))
    members = [r for r in rows(a.members) if r['family'] in {'GH5_2', 'PL1_8', 'PL4_4'}]
    assert len(members) == len({r['protein_id'] for r in members}) == 629
    bmodels = rows(a.blectin_models)
    assert len(bmodels) == len({r['protein_id'] for r in bmodels}) == 28
    assert len({r['species'] for r in bmodels}) == 16
    from summarize_gc3 import ALIASES
    expected = {(r['family'], ALIASES.get(r['species'], r['species'])) for r in members}
    assert len(expected) == 110, len(expected)
    observed = {(r['family'], r['lineage']) for r in rows(grid)
                if r['group'] == 'bacterial_relative' and r['cohort'] == 'closest_relatives'}
    assert expected <= observed, sorted(expected - observed)
    sequence = ''.join(line.strip() for line in cds.read_text().splitlines() if not line.startswith('>'))
    assert len(sequence) == 3609 and set(sequence) <= set('ACGT')
    print(f'PREFLIGHT_PASS sources=38 pinned_files=116 comparator_keys={len(expected)} members=629 blectin=28')


if __name__ == '__main__':
    main()
