#!/usr/bin/env python3
"""Bound the reader payload after assembly, before validation and archive release."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from check_closure import sha


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--payload', type=Path, required=True)
    p.add_argument('--baseline', type=Path, required=True)
    p.add_argument('--receipt', type=Path, required=True)
    a = p.parse_args()
    with a.baseline.open(newline='') as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    assert len(rows) == 1
    baseline = rows[0]
    base_bytes = int(baseline['payload_bytes'])
    files = [path for path in a.payload.rglob('*') if path.is_file()]
    assert files and not any(path.is_symlink() for path in a.payload.rglob('*'))
    candidate_bytes = sum(path.stat().st_size for path in files)
    oversized = [(str(path.relative_to(a.payload)), path.stat().st_size)
                 for path in files if path.stat().st_size > 200_000_000]
    receipt = dict(baseline_main_sha=baseline['main_sha'], baseline_payload_bytes=base_bytes,
                   candidate_payload_bytes=candidate_bytes, growth_bytes=candidate_bytes-base_bytes,
                   max_file_bytes=200_000_000, max_growth_bytes=1_000_000_000,
                   files=len(files), oversized=oversized)
    a.receipt.parent.mkdir(parents=True, exist_ok=True)
    a.receipt.write_text(json.dumps(receipt, indent=2, sort_keys=True) + '\n')
    closure = a.receipt.parent.parent / 'closure'
    token = closure / 'build_token'
    if token.is_file():
        (closure / 'receipt_check_deposit_size.json').write_text(json.dumps({
            'token': token.read_text().strip(),
            'producer': 'scripts/release/check_deposit_size.py',
            'outputs': {a.receipt.relative_to(a.receipt.parent).as_posix(): sha(a.receipt)},
        }, indent=1, sort_keys=True) + '\n')
    if oversized or candidate_bytes > base_bytes + 1_000_000_000:
        raise SystemExit('DEPOSIT_SIZE_FAIL ' + json.dumps(receipt, sort_keys=True))
    print('DEPOSIT_SIZE_PASS ' + json.dumps(receipt, sort_keys=True))


if __name__ == '__main__':
    main()
