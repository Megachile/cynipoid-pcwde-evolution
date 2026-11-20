"""Reject historical or private source pointers across every tracked file."""
import json
import re
import subprocess
from pathlib import Path

FORBIDDEN = re.compile(r'workflow/runs/|\b20\d{6}__|\b(?:RESULT|STATUS|PLAN)\.md|/(?:carc|easley|users|root)/|\bCARC_LOG(?:\.md)?\b|\bSCHEDULER_ATTEMPT(?:\.json)?\b|\bneeds_review\b|\bpipeline_complete\b')
VALIDATOR = 'scripts/release/portable_sources.py'
# Only the exact existing forbidden-pattern declarations may name process tokens.
PATTERN_LINES = {
    'scripts/release/validate_data_deposit.py': '(r"(?i)claude|codex|chatgpt|large language model|needs_review|pipeline_complete", "no agent/process marks"),',
    'scripts/tables/build_table_s6.py': 'forbidden = (b"DCAN", b"Dryocosmus cannoni", b"/" + b"carc/", b"/" + b"easley/", b"needs_review", b"pipeline_complete")',
}


def allowed_pattern(name, line):
    if name == VALIDATOR:
        return line.startswith('FORBIDDEN = re.compile(') or line.strip() in {f'{key!r}: {value!r},' for key, value in PATTERN_LINES.items()}
    return line.strip() == PATTERN_LINES.get(name)



def scan(repo):
    names = subprocess.check_output(['git', '-C', str(repo), 'ls-files'], text=True).splitlines()
    census = []
    for name in names:
        path = repo / name
        if path.is_symlink() or not path.is_file():
            raise ValueError(f'unsafe tracked file: {name}')
        text = path.read_bytes().decode('utf-8', errors='replace')
        hits = [{'line': n, 'text': line} for n, line in enumerate(text.splitlines(), 1) if FORBIDDEN.search(line) and not allowed_pattern(name, line)]
        census.append({'path': name, 'bytes': path.stat().st_size, 'validator_patterns_only': False, 'hits': hits})
    bad = [row for row in census if row['hits']]
    return {'pass': not bad, 'tracked_files': len(names), 'unallowlisted_files': len(bad), 'census': census}


def validate(repo):
    report = scan(repo)
    if not report['pass']:
        raise ValueError('forbidden source pointers: ' + ', '.join(r['path'] for r in report['census'] if r['hits'] and not r['validator_patterns_only']))
    print(f"PORTABLE_SOURCES_PASS files={report['tracked_files']}")
    return report


if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(); p.add_argument('--repo', type=Path, required=True); p.add_argument('--report', type=Path)
    a = p.parse_args(); report = validate(a.repo)
    if a.report:
        a.report.write_text(json.dumps(report, indent=2) + '\n')
