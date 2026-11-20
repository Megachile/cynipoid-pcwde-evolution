"""Resolve portable source names through the committed immutable input inventory."""
import csv
import hashlib
import os
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
FIELDS = ['logical_id', 'relative_path', 'sha256', 'bytes']


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def inventory():
    with (REPO / 'config/source_aliases.tsv').open(newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if reader.fieldnames != FIELDS:
            raise ValueError('source alias inventory schema mismatch')
        rows = list(reader)
    if not rows or len({r['logical_id'] for r in rows}) != len(rows):
        raise ValueError('empty or duplicate source aliases')
    return {r['logical_id']: r for r in rows}


def source_path(base, relative):
    """Aliases never fall back to a live workspace or an unregistered file."""
    name = str(relative)
    if not name.startswith('source_inputs/'):
        return Path(base) / relative
    row = inventory().get(name)
    if row is None:
        raise ValueError(f'missing source alias: {name}')
    configured = os.environ.get('PCWDE_SOURCE_INPUTS')
    if not configured:
        raise ValueError('PCWDE_SOURCE_INPUTS is required for portable source aliases')
    root = Path(configured)
    if root.is_symlink() or not root.is_dir():
        raise ValueError('unsafe source alias root')
    root = root.resolve(strict=True)
    rel = Path(row['relative_path'])
    if rel.is_absolute() or '..' in rel.parts:
        raise ValueError('unsafe source alias inventory path')
    path = root / rel
    if path.is_symlink() or not path.is_file() or not path.resolve().is_relative_to(root):
        raise ValueError(f'missing or escaping pinned source: {name}')
    if path.stat().st_size != int(row['bytes']) or digest(path) != row['sha256']:
        raise ValueError(f'pinned source hash mismatch: {name}')
    return path
