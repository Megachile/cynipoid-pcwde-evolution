"""Resolve the unique pinned supplement author export, confined to the repository."""
from pathlib import Path
from verify_repository import read_tsv
from data_root import sha256

KIND = 'class U supplement draft Markdown export'

def supplement_export(repo, explicit=None):
    repo = Path(repo).resolve()
    rows = [r for r in read_tsv(repo / 'config/upstream_locked_inputs.tsv') if r['kind'] == KIND]
    if len(rows) != 1:
        raise SystemExit('missing or duplicate class U supplement draft export pin')
    row = rows[0]
    relative = Path(row['project_path'])
    if relative.is_absolute() or '..' in relative.parts:
        raise SystemExit('unsafe supplement draft export path')
    source = (repo / relative).resolve()
    if not source.is_relative_to(repo):
        raise SystemExit('unsafe supplement draft export path')
    if not source.is_file():
        raise SystemExit('missing supplement draft export')
    if sha256(source) != row['sha256']:
        raise SystemExit('supplement draft export sha256 or path mismatch')
    if explicit is not None and Path(explicit).resolve() != source:
        raise SystemExit('supplement draft export sha256 or path mismatch')
    return source
