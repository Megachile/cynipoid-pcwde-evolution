#!/usr/bin/env python3
"""Validate and stage the frozen intron layer-one inputs for a paper build."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
from pathlib import Path


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--frozen', type=Path, required=True)
    p.add_argument('--kind', choices=('member_cds', 'host_background', 'annotation_proteins'), required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    with (a.frozen / 'FROZEN_MANIFEST.tsv').open(newline='') as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    selected = [r for r in rows if r['path'].startswith(a.kind + '/')]
    assert selected and len({r['path'] for r in rows}) == len(rows)
    a.out.mkdir(parents=True, exist_ok=False)
    for row in selected:
        rel = Path(row['path'])
        assert not rel.is_absolute() and '..' not in rel.parts
        source = a.frozen / rel
        assert source.is_file() and source.stat().st_size == int(row['bytes']), source
        assert digest(source) == row['digest'], source
        if a.kind == 'host_background' and rel.parts[1] == 'intron_details':
            continue
        shutil.copyfile(source, a.out / rel.relative_to(a.kind))
    if a.kind == 'host_background':
        parts = sorted((a.frozen / a.kind / 'intron_details').glob('*.tsv'))
        assert len(parts) == 38
        with (a.out / 'background_intron_details.tsv').open('wb') as output:
            header = None
            for part in parts:
                with part.open('rb') as input_handle:
                    first = input_handle.readline()
                    if header is None:
                        header = first
                        output.write(first)
                    else:
                        assert first == header
                    shutil.copyfileobj(input_handle, output)
    print(json.dumps(dict(kind=a.kind, files=len(selected)), sort_keys=True))


if __name__ == '__main__':
    main()
