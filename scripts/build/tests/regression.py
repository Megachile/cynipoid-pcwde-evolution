#!/usr/bin/env python3
"""Small real-argv regression for the checked paper build entry points.

Run under the project group. All runner staging uses the committed tiny deposit.
"""
import argparse
import hashlib
import json
import os
import shutil
from pathlib import Path
import subprocess
import sys
import tempfile


def call(*argv):
    return subprocess.run(argv, text=True, capture_output=True,
                          env=dict(os.environ, PYTHONDONTWRITEBYTECODE='1'))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--repo', type=Path, required=True)
    ap.add_argument('--candidate-sha', required=True)
    ap.add_argument('--project-group', type=int, required=True)
    ap.add_argument('--scratch', type=Path, required=True)
    ap.add_argument('--production-inputs', type=Path, required=True)
    ap.add_argument('--unit-only', action='store_true')
    args = ap.parse_args()
    args.scratch.mkdir(parents=True, exist_ok=True)
    args.scratch = args.scratch.resolve()
    with tempfile.TemporaryDirectory(dir=args.scratch) as tmp:
        root = Path(tmp)
        repo = root / 'checkout'
        subprocess.run(['git', 'clone', '-q', '--shared', str(args.repo), str(repo)], check=True)
        subprocess.run(['git', '-C', str(repo), 'checkout', '-q', '--detach', args.candidate_sha], check=True)
        base = subprocess.check_output(['git', '-C', str(repo), 'rev-parse', args.candidate_sha + '^'], text=True).strip()
        main = root / 'main.git'
        subprocess.run(['git', 'clone', '-q', '--bare', '--shared', str(args.repo), str(main)], check=True)
        subprocess.run(['git', '--git-dir', str(main), 'update-ref', 'refs/heads/main', base], check=True)
        deposit = root / 'deposit'
        fixture = Path(__file__).with_name('fixtures')
        subprocess.run(['cp', '-a', str(fixture), str(deposit)], check=True)
        seed = deposit / 'pcwde_supporting_data/build_inputs/companion.txt'
        config = dict(repo=str(repo), base_sha=base, sha=args.candidate_sha,
                      main_git_dir=str(main), deposit=str(deposit), destination=str(root / 'output'),
                      receipt=str(root / 'receipt'), prior_data=str(root), workspace=str(root),
                      python=sys.executable, envs=str(root), provenance_sbatch='sbatch',
                      project_gid=args.project_group)
        cfg = root / 'config.json'
        launcher = repo / 'scripts/build/launch.py'
        runner = repo / 'scripts/build/run_job.py'
        results = []

        def check(name, expected, mode='fixture', extra=(), program=launcher):
            cfg.write_text(json.dumps(config))
            argv = [sys.executable, '-B', str(program), str(cfg), mode, *extra]
            if program == launcher:
                argv.append('--preflight-only')
            result = call(*argv)
            output = result.stdout + result.stderr
            ok = (result.returncode == 0 if expected == 'FIXTURE_GATE_PASS' else
                  (result.returncode == 0 and expected in output) if expected == 'LAUNCH_PREFLIGHT_PASS' else
                  (result.returncode != 0 and expected in output))
            results.append(dict(case=name, pass_=ok, argv=argv, exit=result.returncode,
                                stdout=result.stdout, stderr=result.stderr))
            print(('PASS ' if ok else 'FAIL ') + name + ' exit=' + str(result.returncode), flush=True)
            (args.scratch / 'regression.json').write_text(json.dumps(results, indent=2) + '\n')
            if Path(config['receipt']).exists():
                shutil.copytree(config['receipt'], args.scratch / ('receipt_' + name), dirs_exist_ok=True)
            if not ok:
                raise RuntimeError(name + ': ' + output)

        check('positive', 'LAUNCH_PREFLIGHT_PASS')
        config['project_gid'] += 1
        check('wrong_group', 'PROJECT_GID_FAIL')
        config['project_gid'] = args.project_group
        subprocess.run(['git', '--git-dir', str(main), 'update-ref', 'refs/heads/main', args.candidate_sha], check=True)
        check('main_drift', 'MAIN_PIN_FAIL')
        subprocess.run(['git', '--git-dir', str(main), 'update-ref', 'refs/heads/main', base], check=True)
        config['sha'] = base
        check('sha_drift', 'CHECKOUT_SHA_FAIL')
        config['sha'] = args.candidate_sha
        (repo / 'dirty.fixture').write_text('dirty\n')
        check('dirty_checkout', 'CHECKOUT_DIRTY')
        (repo / 'dirty.fixture').unlink()
        ignored = repo / 'scripts/build/__pycache__'
        ignored.mkdir()
        (ignored / 'launch.cpython-311.pyc').write_bytes(b'ignored bytecode')
        check('ignored_bytecode', 'CHECKOUT_DIRTY')
        (ignored / 'launch.cpython-311.pyc').unlink()
        ignored.rmdir()
        check('production_override', 'PRODUCTION_OVERRIDE_FAIL', 'production', ('--deposit-override', str(deposit)))
        (root / 'output').mkdir()
        check('existing_destination', 'DEST_EXISTS_OR_RELATIVE')
        (root / 'output').rmdir()
        config['destination'] = str(root.parent / 'outside')
        check('outside_destination', 'DEST_OUTSIDE_CHECKOUT_PARENT')
        config['destination'] = str(root / 'output')
        (root / 'output').symlink_to(deposit, target_is_directory=True)
        check('symlink_destination', 'DEST_EXISTS_OR_RELATIVE')
        (root / 'output').unlink()
        (deposit / 'linked').symlink_to(seed)
        check('source_symlink', 'SOURCE_SYMLINK_FAIL', program=runner)
        (deposit / 'linked').unlink()
        if args.unit_only:
            print(f'REGRESSION_UNIT_PASS {len(results)} cases')
            return
        (root / 'receipt').rename(root / 'receipt_symlink')
        source_objects = [deposit, *deposit.rglob('*')]
        stage_bytes = sum(p.stat().st_size for p in source_objects if p.is_file())
        print(f'STAGING_PEAK source_bytes={stage_bytes} staged_bytes={stage_bytes} '
              f'source_objects={len(source_objects)} staged_objects={len(source_objects) + 1}', flush=True)
        # Preserve the original underspecified integration failure and its verifier log.
        check('tiny_production_verifier_rejected', 'COMMAND_FAIL', program=runner)
        (root / 'receipt').rename(root / 'receipt_tiny')
        # Unit guards above use the tiny deposit; integration requires real sealed inputs.
        config['deposit'] = str(args.production_inputs.resolve(strict=True))
        check('staging_and_companions', 'FIXTURE_GATE_PASS', program=runner)
        receipt = root / 'receipt'
        if (receipt / 'result.txt').read_text().strip() != 'FIXTURE_GATE_PASS':
            raise RuntimeError('missing fixture result')
        inventory = json.loads((receipt / 'source_inventory.json').read_text())
        for relative, size, digest in inventory:
            payload = (args.production_inputs / relative).read_bytes()
            if len(payload) != size or hashlib.sha256(payload).hexdigest() != digest:
                raise RuntimeError('source inventory mismatch: ' + relative)
        if 'STAGE_PASS' not in (receipt / 'preflight.txt').read_text():
            raise RuntimeError('missing stage pass')
        output = root / 'output'
        copied = output / 'copied0755'
        if copied.stat().st_mode & 0o777 != 0o755 or (copied / 'seed').read_bytes() != b'seed\n' or \
                (copied / 'companion').read_bytes() != b'companion\n' or \
                (copied / 'companion').stat().st_gid != args.project_group:
            raise RuntimeError('companion mode/content/group mismatch')
        print('PASS staged_hash_agreement and companion_mode_gid', flush=True)
        report = args.scratch / 'regression.json'
        report.write_text(json.dumps(results, indent=2) + '\n')
        print(f'REGRESSION_PASS {len(results)} cases report={report}')


if __name__ == '__main__':
    main()
