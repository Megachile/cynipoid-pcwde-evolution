#!/usr/bin/env python3
"""Run a checked paper build on a compute node with durable receipts."""
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

from launch import check

from lifecycle import inventory, consume_source, consume_comparison

def run(argv, log, env=None):
    with log.open('w') as out:
        out.write('ARGV=' + json.dumps(list(map(str, argv))) + '\n')
        out.flush()
        p = subprocess.run(list(map(str, argv)), stdout=out, stderr=subprocess.STDOUT, env=env)
        out.write('EXIT=' + str(p.returncode) + '\n')
    if p.returncode:
        raise RuntimeError(f'COMMAND_FAIL {log} exit={p.returncode}')

def main():
    config_path, mode, *rest = sys.argv[1:]
    if len(rest) > 1:
        raise ValueError('ARGUMENTS_FAIL')
    cfg = json.loads(Path(config_path).read_text())
    repo, dest, source = check(cfg, mode, rest[0] if rest else None)
    receipt = Path(cfg['receipt'])
    receipt.mkdir(parents=True, exist_ok=False)
    local = None
    try:
        original = inventory(source)
        (receipt / 'source_inventory.json').write_text(json.dumps(original) + '\n')
        py = cfg['python']
        env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1')
        env.pop('DEPOSIT_OVERRIDE', None)
        # Sequence the independently reconstructed stage before build products
        # exist, so it never overlaps the complete output deposit.
        if cfg.get('bounded_staging'):
            comparison = Path(str(dest) + '.stage')
            if comparison.exists():
                raise ValueError('COMPARISON_EXISTS')
            run(['make', '-C', repo, 'stage-check', f'PYTHON={py}',
                 f'BUILD={comparison}', f'DATA_ROOT={source / "pcwde_supporting_data"}',
                 f'BASE={cfg["prior_data"]}', f'WORKSPACE={cfg["workspace"]}'],
                receipt / 'stage_check.log', env)
            if inventory(source) != original:
                raise ValueError('SOURCE_CHANGED_DURING_COMPARE')
            consume_comparison(comparison, repo, inventory(comparison),
                               receipt / 'stage_check.log', receipt)
        local = Path(tempfile.mkdtemp(prefix='pcwde-copytable.', dir='/tmp'))
        staged = local / 'deposit'
        shutil.copytree(source, staged, copy_function=shutil.copy2)
        if inventory(staged) != original or inventory(source) != original:
            raise ValueError('STAGE_OR_SOURCE_CHANGED')
        if mode == "production":
            consume_source(cfg, repo, source, staged, original, receipt)
        data = staged / 'pcwde_supporting_data'
        (receipt / 'preflight.txt').write_text(f'STAGE_PASS source={source} staged={staged} sha={cfg["sha"]}\n')
        py = cfg['python']
        env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', PCWDE_SOURCE_INPUTS=str(data))
        env.pop('DEPOSIT_OVERRIDE', None)
        if mode == 'fixture':
            run(['make', '-C', repo, 'verify', f'PYTHON={py}', f'DATA_ROOT={data}'],
                receipt / 'fixture_verify.log', env)
            run([py, '-B', repo / 'scripts/release/test_check_closure.py'],
                receipt / 'fixture_closure.log', env)
            run([py, '-B', repo / 'scripts/release/figure_build_gate.py',
                 '--repo', repo, '--data-root', data, '--build', dest],
                receipt / 'fixture.log', env)
            # Exercise the staged S9 input and its S10 agreement gate on the
            # canonical launcher's copied deposit, without producing paper data.
            import csv
            import json as _json
            import hashlib
            analysis = data / 'repository_inputs/tables/authority/supplement_new/catalytic_clustering/ANALYSIS.json'
            locked = repo / 'scripts/analysis/catalytic/clustering/LOCKED_ANALYSIS.json'
            if hashlib.sha256(analysis.read_bytes()).digest() != hashlib.sha256(locked.read_bytes()).digest():
                raise ValueError('S9_STAGED_BYTES_FAIL')
            summary = dest / 'fixture_catalytic_summary.tsv'
            rows = []
            for family, record in _json.loads(analysis.read_text())['families'].items():
                if family not in ('GH5_2', 'PL1_8', 'PL4_4'):
                    continue
                for site, call in record['sites'].items():
                    rows.append(dict(alignment_set='de_novo', family=family,
                                     site_label=site, site_role='catalytic',
                                     n_expected_residue=call['mapped_retained'],
                                     n_site_observed=call['mapped_total']))
            with summary.open('w') as handle:
                writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
                writer.writeheader(); writer.writerows(rows)
            cmd = [py, '-B', repo / 'scripts/tables/build_tables_s9_s13.py',
                   '--data-root', data, '--out', dest / 'fixture_tables',
                   '--analysis', analysis, '--catalytic-summary', summary]
            run(cmd, receipt / 'fixture_s9_positive.log', env)
            rows[0]['n_expected_residue'] += 1
            with summary.open('w') as handle:
                writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
                writer.writeheader(); writer.writerows(rows)
            negative = subprocess.run([str(x) for x in cmd], cwd=repo, env=env,
                                      text=True, capture_output=True)
            (receipt / 'fixture_s9_negative.log').write_text(negative.stdout + negative.stderr)
            if negative.returncode == 0 or 'S9/S10 count disagreement' not in negative.stderr + negative.stdout:
                raise ValueError('S9_NEGATIVE_GATE_FAIL')
            # Prove that the closure reader records this staged input as read.
            # Keep the fixture deposit to one file so unrelated production
            # inputs do not need artificial reads or declarations.
            import importlib.util
            closure_path = repo / 'scripts/release/check_closure.py'
            spec = importlib.util.spec_from_file_location('check_closure_fixture', closure_path)
            closure = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(closure)
            mini = local / 'closure_input_fixture'
            mini_input = mini / 'repository_inputs/tables/authority/supplement_new/catalytic_clustering/ANALYSIS.json'
            mini_input.parent.mkdir(parents=True)
            shutil.copy2(analysis, mini_input)
            declarations = local / 'closure_empty_declarations.tsv'
            declarations.write_text('path_pattern\tclass\treason\n')
            closure.reset(mini, local / 'closure_input_result')
            mini_input.read_bytes()
            closure.check_inputs(mini, local / 'closure_input_result', declared=declarations)
            if 'ANALYSIS.json\tread\t' not in (local / 'closure_input_result/inputs.tsv').read_text():
                raise ValueError('S9_CLOSURE_READ_FAIL')
            shutil.copy2(local / 'closure_input_result/inputs.tsv', receipt / 'fixture_closure_inputs.tsv')
            intent = json.loads((dest / 'figure_render_intent.json').read_text())
            if intent['commit'] != cfg['sha']:
                raise ValueError('FIGURE_INTENT_FAIL')
            src = dest / 'source0755'
            src.mkdir(mode=0o755)
            (src / 'seed').write_bytes(b'seed\n')
            copied = dest / 'copied0755'
            shutil.copytree(src, copied, copy_function=shutil.copy2)
            copied.chmod(copied.stat().st_mode & ~0o2000)
            (copied / 'companion').write_bytes(b'companion\n')
            if (copied / 'seed').read_bytes() != b'seed\n' or \
                    copied.stat().st_mode & 0o777 != 0o755 or \
                    (copied / 'companion').stat().st_gid != int(cfg['project_gid']):
                raise ValueError('COPYTREE_GROUP_MODE_FAIL')
            (receipt / 'result.txt').write_text('FIXTURE_GATE_PASS\n')
        else:
            run(['make', '-C', repo, 'verify', f'PYTHON={py}', f'DATA_ROOT={data}'],
                receipt / 'verify.log', env)
            if 'VERIFY_PASS' not in (receipt / 'verify.log').read_text():
                raise ValueError('VERIFY_MARKER_MISSING')
            fc = local / 'fc'
            run(['bash', repo / 'scripts/release/local_env.bash', fc,
                 Path(cfg['envs']) / 'busco_env', Path(cfg['envs']) / 'ggtree_env'],
                receipt / 'environment.log', env)
            env.update(TEXMFVAR=str(local / 'texmfvar'),
                       TEXMFCACHE=str(local / 'texmfcache'),
                       XDG_CACHE_HOME=str(local / 'xdg_cache'))
            for name in ('TEXMFVAR', 'TEXMFCACHE', 'XDG_CACHE_HOME'):
                Path(env[name]).mkdir()
            target = cfg.get('make_target', 'manuscript')
            if target not in ('manuscript', 'paper'):
                raise ValueError('MAKE_TARGET_FAIL')
            make = ['make', '-C', repo, target, f'PYTHON={fc / "bin/python"}',
                    f'RSCRIPT={fc / "bin/Rscript"}', f'BUILD={dest}', f'DATA_ROOT={data}',
                    f'GS={Path(cfg["envs"]) / "meme_env/bin/gs"}']
            run(make, receipt / 'make.log', env)
            if not cfg.get('bounded_staging'):
                run(['make', '-C', repo, 'stage-check', f'PYTHON={fc / "bin/python"}',
                     f'BUILD={dest}.stage', f'DATA_ROOT={data}',
                     f'BASE={cfg["prior_data"]}', f'WORKSPACE={cfg["workspace"]}'],
                    receipt / 'stage_check.log', env)
            if subprocess.check_output(['git', '-C', repo, 'status', '--porcelain', '--ignored']).strip():
                raise ValueError('CHECKOUT_DIRTY_AFTER_BUILD')
            (dest / 'CHECKED_OUT_SHA').write_text(cfg['sha'] + '\n')
            if (dest / 'CHECKED_OUT_SHA').stat().st_gid != int(cfg['project_gid']):
                raise ValueError('OUTPUT_GID_FAIL')
            (receipt / 'result.txt').write_text('BUILD_PASS sha=' + cfg['sha'] + '\n')
        shutil.rmtree(local)
        (receipt / 'scratch_cleanup.txt').write_text(f'validated deletion {local}\n')
    except Exception as exc:
        (receipt / 'failure.txt').write_text(repr(exc) + '\n')
        raise

if __name__ == '__main__':
    main()
