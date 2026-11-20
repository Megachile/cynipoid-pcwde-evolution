#!/usr/bin/env python3
"""Add accepted annotation proteins and bacterial comparators to an A2 freeze."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
from pathlib import Path


PINNED = {
    'bacterial_grid.tsv': '983a790cd4c5a1cfa3e9cbd9fd928d80fee9d660114ab4c8453eb82c5f56ce56',
    'blectin_bacterial_cds.fna': '09bc2bed98e4dc6af53aee753507dbcaea6b18dea16aa765820277fa36381d37',
}


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--frozen', type=Path, required=True)
    p.add_argument('--data-root', type=Path, required=True)
    p.add_argument('--build', type=Path, required=True)
    a = p.parse_args()
    producer = json.loads((a.frozen / 'PRODUCER.json').read_text())
    assert producer['focal_copies'] == 631 and producer['cohort'] == 603
    with (a.frozen / 'FROZEN_MANIFEST.tsv').open(newline='') as handle:
        old = list(csv.DictReader(handle, delimiter='\t'))
    assert len(old) == 85 and all(digest(a.frozen / r['path']) == r['sha256'] for r in old)
    destination = a.frozen / 'annotation_proteins'
    destination.mkdir(exist_ok=False)
    source = a.build / 'introns/annotation_proteins'
    qc = json.loads((source / 'QC.json').read_text())
    assert (qc['sources'], qc['copies'], qc['identical'], qc['different']) == (38, 603, 594, 9)
    for name in ('GH5_2_current_genome_proteins.faa', 'PL1_8_current_genome_proteins.faa',
                 'PL4_4_current_genome_proteins.faa', 'protein_model_comparison.tsv', 'QC.json'):
        shutil.copyfile(source / name, destination / name)
    old_sources = a.data_root / 'repository_inputs/analysis/introns/sources'
    for name, expected in PINNED.items():
        item = old_sources / name
        assert digest(item) == expected
        shutil.copyfile(item, a.frozen / name)
    paths = sorted(path for path in a.frozen.rglob('*') if path.is_file()
                   and path.name not in ('FROZEN_MANIFEST.tsv', 'PRODUCER.json'))
    rows = [dict(path=str(path.relative_to(a.frozen)), bytes=path.stat().st_size, digest=digest(path))
            for path in paths]
    assert len(rows) == 92 and all(r['bytes'] <= 200_000_000 for r in rows)
    with (a.frozen / 'FROZEN_MANIFEST.tsv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, ('path', 'bytes', 'digest'), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    producer['augmentation'] = 'upstream/introns/augment_layer1.py'
    producer['annotation_protein_qc'] = qc
    (a.frozen / 'PRODUCER.json').write_text(json.dumps(producer, indent=2, sort_keys=True) + '\n')
    print(json.dumps(dict(files=len(rows), bytes=sum(r['bytes'] for r in rows)), sort_keys=True))


if __name__ == '__main__':
    main()
