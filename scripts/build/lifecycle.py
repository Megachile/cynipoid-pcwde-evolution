"""Confined cleanup of explicitly declared, verified transient staging trees."""
import hashlib
import json
from pathlib import Path
import shutil


def inventory(root):
    rows = []
    for p in sorted(root.rglob('*')):
        if p.is_symlink():
            raise ValueError('SOURCE_SYMLINK_FAIL')
        if p.is_file():
            h = hashlib.sha256()
            with p.open('rb') as stream:
                for block in iter(lambda: stream.read(1024 * 1024), b''):
                    h.update(block)
            rows.append((str(p.relative_to(root)), p.stat().st_size, h.hexdigest()))
    return rows


def consume_source(cfg, repo, source, staged, original, receipt):
    """After independent copy verification, consume only the declared task source.

    Default callers never delete their source. Opt-in callers supply a sealed
    full inventory and an explicit source path. The source must be an immediate
    task directory next to the checkout; protected/shared data cannot qualify.
    """
    spec = cfg.get('bounded_staging')
    if not spec:
        return
    parent = repo.parent.resolve()
    declared = Path(spec['source'])
    seal = Path(spec['inventory'])
    forbidden = {repo.resolve(), Path(cfg['destination']).resolve(),
                 Path(cfg['receipt']).resolve(), parent / 'genomes', parent / 'good_data'}
    if source != declared or source.parent != parent or source in forbidden or \
            declared.is_symlink() or seal.is_symlink() or seal.parent != receipt.parent or \
            not source.name.startswith('staged_input_'):
        raise ValueError('SOURCE_CLEANUP_CONFINEMENT_FAIL')
    if hashlib.sha256(seal.read_bytes()).hexdigest() != spec['inventory_sha256']:
        raise ValueError('SOURCE_CLEANUP_SEAL_FAIL')
    expected = [tuple(row) for row in json.loads(seal.read_text())]
    if expected != original or inventory(source) != expected or inventory(staged) != expected:
        raise ValueError('SOURCE_CLEANUP_IDENTITY_FAIL')
    event = {'path': str(source), 'bytes': sum(row[1] for row in expected),
             'files': len(expected), 'inventory_sha256': spec['inventory_sha256'],
             'gate': 'source, copied stage and sealed inventory all identical'}
    (receipt / 'source_cleanup_pending.json').write_text(json.dumps(event, indent=2) + '\n')
    shutil.rmtree(source)
    if source.exists():
        raise ValueError('SOURCE_CLEANUP_INCOMPLETE')
    event['disposition'] = 'validated deletion'
    (receipt / 'source_cleanup.json').write_text(json.dumps(event, indent=2) + '\n')


def consume_comparison(path, repo, expected, log, receipt):
    """Release a successful independent restage only after full equality."""
    if path.parent != repo.parent or not path.name.endswith('.stage') or path.is_symlink():
        raise ValueError('COMPARE_CLEANUP_CONFINEMENT_FAIL')
    text = log.read_text()
    if 'STAGE_COMPARE_PASS' not in text or 'differing=0' not in text:
        raise ValueError('COMPARE_CLEANUP_GATE_FAIL')
    if inventory(path) != expected:
        raise ValueError('COMPARE_CLEANUP_CHANGED')
    event = {'path': str(path), 'bytes': sum(row[1] for row in expected),
             'files': len(expected), 'gate': 'full STAGE_COMPARE_PASS and final inventory'}
    (receipt / 'comparison_cleanup_pending.json').write_text(json.dumps(event, indent=2) + '\n')
    shutil.rmtree(path)
    event['disposition'] = 'validated deletion'
    (receipt / 'comparison_cleanup.json').write_text(json.dumps(event, indent=2) + '\n')
