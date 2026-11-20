"""Verify a location-only manifest conversion against immutable original provenance."""
import csv
import hashlib
import json
from pathlib import Path
from source_aliases import source_path, digest, inventory


def validate_bridge(repo, portable, expected_original):
    contract = json.loads((repo / 'config/intron_manifest_bridge.json').read_text())
    if expected_original != contract['original_sha256']:
        raise ValueError('frozen producer manifest identity mismatch')
    original = source_path(Path('/nonexistent'), contract['original_alias'])
    if digest(original) != expected_original or digest(portable) != contract['portable_sha256']:
        raise ValueError('manifest bridge digest mismatch')
    def read(path):
        with path.open(newline='') as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            return reader.fieldnames, list(reader)
    af, a = read(original); bf, b = read(portable)
    if af != bf or len(a) != len(b) or len(a) != contract['rows']:
        raise ValueError('manifest bridge schema or membership mismatch')
    expected = {}
    for change in contract['changes']:
        key = (change['row'], change['field'])
        if key in expected or not change['field'].endswith('_project_path'):
            raise ValueError('duplicate or non-location bridge change')
        expected[key] = change
    observed = set(); aliases = inventory()
    for i, (before, after) in enumerate(zip(a, b)):
        for field in af:
            if before[field] == after[field]:
                continue
            key = (i, field); change = expected.get(key)
            if (change is None or hashlib.sha256(before[field].encode()).hexdigest() != change['original_sha256']
                    or after[field] != change['portable']):
                raise ValueError('manifest bridge altered row, parameter or location')
            name = after[field]
            if name not in aliases:
                raise ValueError('manifest bridge missing alias')
            payload = source_path(Path('/nonexistent'), name)
            pin = before[field.replace('_project_path', '_sha256')]
            if digest(payload) != pin:
                raise ValueError('manifest bridge payload mismatch')
            observed.add(key)
    if observed != set(expected):
        raise ValueError('manifest bridge incomplete mapping')
    return {'rows': len(a), 'location_changes': len(observed), 'original_sha256': expected_original}
