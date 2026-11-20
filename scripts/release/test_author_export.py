#!/usr/bin/env python3
"""Fail-closed author-export selection using the real manifest loader."""
import csv
import shutil
import tempfile
from pathlib import Path
from author_export import supplement_export, KIND

repo = Path(__file__).resolve().parents[2]
source = supplement_export(repo)
assert supplement_export(repo, source) == source
with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'config').mkdir()
    pin = root / 'config/upstream_locked_inputs.tsv'
    with (repo / 'config/upstream_locked_inputs.tsv').open() as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        fields, rows = reader.fieldnames, list(reader)
    row = next(r for r in rows if r['kind'] == KIND)
    dest = root / row['project_path']
    dest.parent.mkdir(parents=True)
    shutil.copyfile(source, dest)
    def write(records):
        with pin.open('w') as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t')
            writer.writeheader(); writer.writerows(records)
    def reject(records, expected, explicit=None):
        write(records)
        try:
            supplement_export(root, explicit)
        except SystemExit as exc:
            assert expected in str(exc), str(exc)
        else:
            raise AssertionError('invalid export accepted')
    reject([], 'missing or duplicate')
    reject([row, row], 'missing or duplicate')
    reject([dict(row, project_path='authority/missing.md')], 'missing supplement')
    reject([dict(row, project_path='../escape.md')], 'unsafe')
    reject([dict(row, project_path=str(source))], 'unsafe')
    reject([dict(row, sha256='0'*64)], 'sha256 or path mismatch')
    reject([row], 'sha256 or path mismatch', source)
    link = root / 'outside.md'; link.symlink_to(source)
    reject([dict(row, project_path='outside.md')], 'unsafe')
print('AUTHOR_EXPORT_TEST_PASS cases=10')
