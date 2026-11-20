#!/usr/bin/env python3
"""Post-computation packaging fixture on an isolated retained-build copy.

The caller supplies a writable copy of the retained A4 build and its pinned data
root. No genome assemblies or upstream computations are read or repeated here.
"""
import argparse
import csv
import os
from pathlib import Path
import subprocess
import sys

from data_root import validate_reader_inputs
from package_archive import validate_archive

REPO = Path(__file__).resolve().parents[2]


def run(*args, env=None, phrase=None):
    result = subprocess.run([str(x) for x in args], env=env, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if result.returncode or (phrase and phrase not in result.stdout):
        raise SystemExit(f"PACKAGING_FIXTURE_FAIL {' '.join(map(str, args))}:\n{result.stdout}")
    print(result.stdout, end="")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--data-root', type=Path, required=True)
    parser.add_argument('--build', type=Path, required=True,
                        help='writable isolated copy of retained A4 build, with data_deposit absent')
    args = parser.parse_args()
    data = args.data_root.resolve()
    build = args.build.resolve()
    if not data.is_dir() or not build.is_dir() or (build / 'data_deposit').exists():
        raise SystemExit('PACKAGING_FIXTURE_FAIL requires data root and isolated build without deposit')
    with (REPO / 'config/logical_inputs.tsv').open(newline='') as handle:
        logical = list(csv.DictReader(handle, delimiter='\t'))
    expected = {r['figshare_relative_path'] for r in logical}
    if len(logical) != len(expected) or not logical:
        raise SystemExit('PACKAGING_FIXTURE_FAIL logical input inventory invalid')
    if len(validate_reader_inputs(REPO, data)) != len(logical):
        raise SystemExit('PACKAGING_FIXTURE_FAIL pinned input count differs')
    closure = build / 'closure'
    if not (closure / 'build_token').is_file() or not (closure / 'inputs.tsv').is_file():
        raise SystemExit('PACKAGING_FIXTURE_FAIL missing retained closure inputs/token')
    env = {**os.environ, 'PCWDE_REPO_ROOT': str(REPO), 'PCWDE_DATA_ROOT': str(data),
           'PCWDE_DEPOSIT_OUT': str(build / 'data_deposit'), 'PYTHONDONTWRITEBYTECODE': '1'}
    python = sys.executable
    run(python, '-B', REPO / 'scripts/release/build_data_deposit.py', env=env,
        phrase='package_files=')
    run(python, '-B', REPO / 'scripts/release/check_deposit_size.py', '--payload',
        build / 'data_deposit/pcwde_supporting_data', '--baseline',
        REPO / 'config/deposit_size_baseline.tsv', '--receipt',
        build / 'data_deposit/DEPOSIT_SIZE.json', env=env, phrase='DEPOSIT_SIZE_PASS')
    run(python, '-B', REPO / 'scripts/release/check_closure.py', '--deposit',
        '--data-root', data, '--build', build, env=env, phrase='CLOSURE_DEPOSIT_PASS')
    run(python, '-B', REPO / 'scripts/release/validate_data_deposit.py', env=env,
        phrase='reader-data, archive, and scientific-value gates passed')
    if validate_archive(build / 'data_deposit') <= len(logical):
        raise SystemExit('PACKAGING_FIXTURE_FAIL archive member coverage')
    # A new generated file must still fail closure, even when its ctime is current.
    rogue = build / 'data_deposit/UNDECLARED_GENERATED.txt'
    rogue.write_text('fixture rejection\n')
    result = subprocess.run([python, '-B', str(REPO / 'scripts/release/check_closure.py'),
                             '--deposit', '--data-root', str(data), '--build', str(build)],
                            env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    rogue.unlink()
    if result.returncode == 0 or 'no deposit-layout rule explains this file' not in result.stdout:
        raise SystemExit('PACKAGING_FIXTURE_FAIL undeclared generated file was not rejected:\n' + result.stdout)
    print(f'PACKAGING_COMMAND_FIXTURE_PASS logical_inputs={len(logical)} generated_receipts=yes '
          'closure=yes validator=yes archive=yes negative_undeclared=yes')


if __name__ == '__main__':
    main()
