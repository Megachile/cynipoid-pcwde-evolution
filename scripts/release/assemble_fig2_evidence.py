#!/usr/bin/env python3
"""Build-bound Figure 2 review previews and evidence packet."""
import argparse
import csv
import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path

BASE = 'd72d6d02d0360d3eb4c8aa27ca24f59903a3d7d4'
DRAFT_SHA = '40b1a497ef6308d64e069af547f2a156e87b67ce8ed0ccbaf39ad104c8cf5e63'
FAMILIES = ('GH5_2', 'PL1_8', 'PL4_4')

def digest(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def need(p):
    p = Path(p)
    if not p.is_file() or p.stat().st_size == 0:
        raise SystemExit(f'missing or empty evidence: {p}')
    return p

def git(repo, *args):
    return subprocess.check_output(['git', '-C', str(repo), *args], text=True).strip()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('action', choices=('preflight', 'preview', 'packet'))
    for arg in ('repo', 'build', 'figure', 'out', 'baseline', 'closure', 'build-log', 'frozen', 'expected-sha', 'comparison', 'invocation'):
        ap.add_argument('--' + arg)
    a = ap.parse_args()
    repo = Path(a.repo).resolve(); build = Path(a.build).resolve()
    draft = need(repo / 'config/fig2_preview_draft.md')
    if digest(draft) != DRAFT_SHA:
        raise SystemExit('preview draft hash mismatch')
    if 'PCWDE outgroup placement trees' not in draft.read_text():
        raise SystemExit('preview caption missing')
    if a.expected_sha and git(repo, 'rev-parse', 'HEAD') != a.expected_sha:
        raise SystemExit('checkout SHA mismatch')
    if a.action == 'preflight':
        for p in ('scripts/trees/prune_origin_display_trees.R', 'scripts/release/check_origin_key_nodes.py', 'scripts/figures/mbe_page_preview.py', 'config/origin_display_tips.tsv', 'docs/ORIGIN_DISPLAY_TREES.md'):
            need(repo / p)
        rows = list(csv.DictReader(need(repo / 'config/logical_inputs.tsv').open(), delimiter='\t'))
        for f in FAMILIES:
            matches = [r for r in rows if r['figshare_relative_path'].endswith(f'origin_context_{f}_full.treefile')]
            if len(matches) != 1 or len(matches[0]['sha256']) != 64:
                raise SystemExit(f'full-tree pin missing: {f}')
        print('FIG2_PREFLIGHT_PASS sha=' + git(repo, 'rev-parse', 'HEAD'))
        return
    if a.action == 'preview':
        fig = need(a.figure)
        if a.baseline:
            receipt = json.loads(need(Path(a.baseline) / 'figure_build_receipt.json').read_text())
            if receipt.get('commit') != BASE:
                raise SystemExit('baseline build commit mismatch')
            if digest(fig) != receipt['figures']['origin_trees/origin_trees.pdf']:
                raise SystemExit('baseline figure hash mismatch')
        out = Path(a.out); out.parent.mkdir(parents=True, exist_ok=True)
        subprocess.run([sys.executable, str(repo / 'scripts/figures/mbe_page_preview.py'), '--fig', str(fig), '--figure-number', '2', '--caption-match', 'PCWDE outgroup placement trees', '--manuscript', str(draft), '--out', str(out)], check=True)
        for suffix in ('.pdf', '.png', '.json'): need(str(out) + suffix)
        print('FIG2_PREVIEW_PASS ' + str(out))
        return
    if git(repo, 'rev-parse', BASE) != BASE or git(repo, 'status', '--porcelain'):
        raise SystemExit('checkout/base not clean')
    sha = git(repo, 'rev-parse', 'HEAD')
    if a.expected_sha != sha: raise SystemExit('packet checkout SHA mismatch')
    closure = need(a.closure); log = need(a.build_log)
    if 'BUILD_PASS sha=' + sha not in log.read_text(): raise SystemExit('build pass absent')
    if 'CLOSURE' not in closure.read_text() and '== 1.' not in closure.read_text(): raise SystemExit('closure comparison absent')
    base = Path(a.baseline); frozen = Path(a.frozen)
    receipt = json.loads(need(base / 'figure_build_receipt.json').read_text())
    if receipt.get('commit') != BASE: raise SystemExit('baseline identity mismatch')
    required = [build / 'trees/origin_key_nodes.tsv', build / 'trees/origin_display_support_check.tsv', build / 'numbers.json', build / 'figures/origin_trees/origin_trees.pdf', build / 'figures/origin_trees/origin_trees_full.pdf']
    for f in FAMILIES:
        for v in ('fig2', 's1'):
            required.append(build / f'trees/origin_display_{f}_{v}.treefile')
            required.append(frozen / f'{f}_{"pruned" if v == "fig2" else "s1"}.treefile')
    for p in required: need(p)
    for stem in ('fig2_candidate_page', 'fig2_before_page'):
        for ext in ('.png', '.pdf', '.json'): need(build / 'previews' / (stem + ext))
    rows = list(csv.DictReader((build / 'trees/origin_key_nodes.tsv').open(), delimiter='\t'))
    if len(rows) != 12 or {(r['family'], r['view']) for r in rows} != {(f, v) for f in FAMILIES for v in ('fig2', 's1')}:
        raise SystemExit('key-node coverage mismatch')
    if 'ORIGIN_KEY_CONSISTENCY_PASS' not in log.read_text() or 'CLOSURE_DEPOSIT_PASS' not in log.read_text():
        raise SystemExit('required checks absent from build log')
    packet = build / 'fig2_evidence'; packet.mkdir(exist_ok=True)
    sources = [repo / p for p in ('scripts/trees/prune_origin_display_trees.R', 'scripts/trees/compare_frozen_origin_views.R', 'scripts/figures/origin_trees.R', 'scripts/release/check_origin_key_nodes.py', 'scripts/release/assemble_fig2_evidence.py', 'config/origin_display_tips.tsv', 'config/logical_inputs.tsv', 'docs/ORIGIN_DISPLAY_TREES.md')]
    comparison = need(a.comparison)
    invocation = need(a.invocation)
    if sha not in invocation.read_text(): raise SystemExit('invocation SHA mismatch')
    sources += required + [closure, log, comparison, invocation, base / 'figure_build_receipt.json']
    sources += [build / 'previews' / (stem + ext) for stem in ('fig2_candidate_page', 'fig2_before_page') for ext in ('.png', '.pdf', '.json')]
    sources += [repo / 'Makefile']
    manifest = []
    for i, src in enumerate(sources):
        dst = packet / f'{i:03d}_{src.name}'
        shutil.copy2(src, dst)
        manifest.append((dst.name, dst.stat().st_size, digest(dst), str(src)))
    diff = git(repo, 'diff', '--stat', BASE + '..' + sha)
    (packet / 'diff_stat.txt').write_text(diff + '\n')
    names = git(repo, 'diff', '--name-only', BASE + '..' + sha).splitlines()
    (packet / 'changed_files.tsv').write_text('path\tlayer\n' + ''.join(p + '\t' + ('layer_2_producer' if p.startswith(('scripts/trees/', 'config/origin_display', 'docs/ORIGIN')) else 'layer_3_renderer_or_build' if p.startswith(('scripts/figures/', 'scripts/release/', 'Makefile')) else 'input_manifest_or_presentation') + '\n' for p in names))
    (packet / 'identity.json').write_text(json.dumps({'candidate': sha, 'baseline': BASE, 'preview_draft_sha256': DRAFT_SHA}, indent=2) + '\n')
    for p in ('diff_stat.txt', 'changed_files.tsv', 'identity.json'):
        f = packet / p; manifest.append((p, f.stat().st_size, digest(f), str(f)))
    with (packet / 'INVENTORY.tsv').open('w') as out:
        out.write('file\tbytes\tsha256\tsource\n')
        for row in manifest: out.write('\t'.join(map(str, row)) + '\n')
    print('FIG2_PACKET_PASS sha=' + sha + ' files=' + str(len(manifest)))
if __name__ == '__main__': main()
