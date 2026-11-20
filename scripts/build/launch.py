#!/usr/bin/env python3
"""Checked entry for the canonical paper build Slurm job."""
import argparse
import json
import hashlib
import os
from pathlib import Path
import subprocess
import sys

REQUIRED = {'repo', 'base_sha', 'sha', 'main_git_dir', 'deposit',
            'destination', 'receipt', 'prior_data', 'workspace', 'python',
            'envs', 'provenance_sbatch', 'project_gid'}

def git(*args):
    return subprocess.check_output(['git', *map(str, args)], text=True).strip()

def check(config, mode, override=None, group=True):
    if mode not in ('fixture', 'production'):
        raise ValueError('MODE_FAIL')
    if mode == 'production' and (override or os.environ.get('DEPOSIT_OVERRIDE')):
        raise ValueError('PRODUCTION_OVERRIDE_FAIL')
    if not REQUIRED <= config.keys():
        raise ValueError('CONFIG_FIELDS_FAIL')
    repo = Path(config['repo']).resolve(strict=True)
    dest = Path(config['destination'])
    if not dest.is_absolute() or dest.is_symlink() or dest.exists():
        raise ValueError('DEST_EXISTS_OR_RELATIVE')
    if dest.name in ('', '.', '..') or dest.resolve() != repo.parent / dest.name:
        raise ValueError('DEST_OUTSIDE_CHECKOUT_PARENT')
    if group and os.getgid() != int(config['project_gid']):
        raise ValueError('PROJECT_GID_FAIL')
    if git('--git-dir', config['main_git_dir'], 'rev-parse', 'main') != config['base_sha']:
        raise ValueError('MAIN_PIN_FAIL')
    if git('-C', repo, 'rev-parse', 'HEAD') != config['sha']:
        raise ValueError('CHECKOUT_SHA_FAIL')
    if os.environ.get('PYTHONDONTWRITEBYTECODE') != '1':
        raise ValueError('BYTECODE_ENV_FAIL')
    if git('-C', repo, 'status', '--porcelain', '--ignored'):
        raise ValueError('CHECKOUT_DIRTY')
    source = Path(override or config['deposit']).resolve(strict=True)
    if not (source / 'pcwde_supporting_data/repository_inputs/INPUT_MANIFEST.tsv').is_file():
        raise ValueError('MANIFEST_MISSING')
    if config.get('bounded_staging'):
        spec = config['bounded_staging']
        seal = Path(spec['inventory'])
        if (Path(spec['source']) != source or source.parent != repo.parent or
                not source.name.startswith('staged_input_') or seal.is_symlink() or
                seal.parent != repo.parent):
            raise ValueError('SOURCE_CLEANUP_CONFINEMENT_FAIL')
        if hashlib.sha256(seal.read_bytes()).hexdigest() != spec['inventory_sha256']:
            raise ValueError('SOURCE_CLEANUP_SEAL_FAIL')
    return repo, dest, source

def main():
    p = argparse.ArgumentParser()
    p.add_argument('config')
    p.add_argument('mode', choices=('fixture', 'production'))
    p.add_argument('--deposit-override')
    p.add_argument('--preflight-only', action='store_true')
    a = p.parse_args()
    cfg = json.loads(Path(a.config).read_text())
    check(cfg, a.mode, a.deposit_override)
    if a.preflight_only:
        print('LAUNCH_PREFLIGHT_PASS')
        return 0
    job = Path(__file__).with_name('job.sbatch')
    options = []
    if a.mode == 'fixture':
        options = (['--time=00:30:00', '--cpus-per-task=2', '--mem=8G']
                   if cfg.get('bounded_staging') else
                   ['--time=00:10:00', '--cpus-per-task=1', '--mem=1G'])
    argv = [cfg['provenance_sbatch'], *options, '--', str(job), str(Path(a.config).resolve()),
            a.mode]
    if a.deposit_override:
        argv.append(str(Path(a.deposit_override).resolve()))
    print('SUBMIT_ARGV=' + json.dumps(argv), flush=True)
    return subprocess.call(argv)

if __name__ == '__main__':
    sys.exit(main())
