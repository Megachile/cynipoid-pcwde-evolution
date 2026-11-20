#!/usr/bin/env python3
"""Export the supplement assembler's section text as an editable authoring seed."""
import argparse
import csv
import html
import json
import re
from pathlib import Path

p = argparse.ArgumentParser()
p.add_argument('--build', type=Path, required=True)
p.add_argument('--repo', type=Path, required=True)
a = p.parse_args()
sections = json.loads((a.build / 'supplement/section_text.json').read_text())
numbers = json.loads((a.build / 'numbers.json').read_text())
with (a.build / 'table_s2/candidate/supplementary_tables/Table_S2.tsv').open() as f:
    s2 = list(csv.DictReader(f, delimiter='\t'))
with (a.build / 'table_s6/supplementary_tables/Table_S6.tsv').open() as f:
    s6 = list(csv.DictReader(f, delimiter='\t'))
with (a.build / 'figures/b_lectin_tree/evidence/drawn_placements.tsv').open() as f:
    locus_tips = list(csv.DictReader(f, delimiter='\t'))
with (a.build / 'figures/b_lectin_domain_tree/primary/sidecars/Figure_S4_domain_manifest.tsv').open() as f:
    domain_records = list(csv.DictReader(f, delimiter='\t'))
with (a.repo / 'config/supplement_seed_checks.tsv').open() as f:
    checks = list(csv.DictReader(f, delimiter='\t'))

def value(source):
    if source == 'pending':
        return 'pending until synteny step B merges'
    if source.startswith('numbers:'):
        keys = source.split(':', 1)[1].split(',')
        return ', '.join(f'{key}={numbers[key.split(".")[0]][key.split(".")[1]]}' for key in keys)
    if source == 'b_lectin:locus_tips':
        focal = sum(r['tip_role'] != 'outgroup' for r in locus_tips)
        outgroup = sum(r['tip_role'] == 'outgroup' for r in locus_tips)
        return f'Figure 5 locus-tree tips: {focal} cynipoid loci + {outgroup} outgroup; Table S2 counts all identified proteins, a different unit'
    if source == 'b_lectin:domain_records':
        focal = sum(r['tip_class'] == 'focal' for r in domain_records)
        outgroup = sum(r['tip_class'] == 'bacterial' for r in domain_records)
        if focal + outgroup != len(domain_records):
            raise ValueError('unexpected Figure S4 domain tip class')
        return f'Figure S4 domain-tree records: {len(domain_records)} total ({focal} cynipoid, {outgroup} outgroup)'
    if source == 'table_s2:B_lectin':
        (row,) = [r for r in s2 if r['family_id'] == 'B_lectin']
        return f"{row['focal_count']} focal proteins + {row['outgroup_count']} outgroup; {row['phylogeny_unit_count']}"
    if source == 'table_s6:tissues':
        tissues = sorted({r['tissue'] for r in s6 if r['life_stage'] == 'larva'})
        return 'larval tissues: ' + ', '.join(tissues)
    if source == 'table_s6:adult_species':
        species = {r['species'] for r in s6 if r['life_stage'] == 'adult' and r['tissue'] in ('venom_gland', 'accessory_gland')}
        return f'{len(species)} species with adult venom/accessory rows'
    if source == 'table_s6:Biorhiza_pallida':
        tissues = sorted({r['tissue'] for r in s6 if r['species'] == 'Biorhiza pallida'})
        return 'Biorhiza pallida Table S6 rows (archival): ' + (', '.join(tissues) if tissues else 'none')
    raise ValueError(source)

def render(ss):
    lines = []
    for s in ss:
        lines += ['## ' + s['heading'], '']
        if s['author_slot']:
            lines += [s['author_slot'], '']
        for paragraph in s['paragraphs']:
            lines += [paragraph, '']
        for item in s['items']:
            kind, ident = item.split(':', 1)
            if kind == 'figure':
                with (a.repo / 'config/figure_numbers.tsv').open() as f:
                    (row,) = [r for r in csv.DictReader(f, delimiter='\t') if r['figure_id'] == ident]
                label = 'Figure ' + row['label']
            else:
                with (a.repo / 'config/table_numbers.tsv').open() as f:
                    (row,) = [r for r in csv.DictReader(f, delimiter='\t') if r['table_id'] == ident]
                label = 'Table ' + row['label']
            lines += [f'**{label}.** [Item label; source artwork/table remains in the assembled supplement]', '']
    return '\n'.join(lines).rstrip() + '\n'

baseline = render(sections)
for check in checks:
    matches = [(s, i) for s in sections if s['section_id'] == check['section_id'] for i, para in enumerate(s['paragraphs']) if check['anchor'] in para]
    if len(matches) != 1:
        raise SystemExit(f"check {check['item']} anchor occurs {len(matches)} times")
    s, i = matches[0]
    fixed = ('; ' + check['fixed_note']) if check.get('fixed_note') else ''
    note = f" [CHECK #{check['item']}: {check['label']}: {value(check['value_source'])}{fixed}]"
    s['paragraphs'][i] = s['paragraphs'][i].replace(check['anchor'], check['anchor'] + note, 1)
seed = render(sections)
stripped = re.sub(r' \[CHECK #\d+: [^\]]*\]', '', seed)
if stripped != baseline:
    raise SystemExit('PREFIX_TEST_FAIL: notes alter assembled section text')
for s in sections:
    if ('## ' + s['heading']) not in seed:
        raise SystemExit('HEADING_TEST_FAIL')
out = a.build / 'supplement_draft_seed'
out.mkdir(exist_ok=True)
(out / 'supplement_draft_seed.md').write_text(seed)
body = '<!doctype html><html lang="en"><meta charset="utf-8"><title>Supplement draft seed</title><body>'
for block in seed.split('\n\n'):
    if block.startswith('## '):
        body += '<h2>' + html.escape(block[3:]) + '</h2>'
    else:
        body += '<p>' + html.escape(block) + '</p>'
(out / 'supplement_draft_seed.html').write_text(body + '</body></html>\n')
(out / 'prefix_test.txt').write_text('PASS: stripping CHECK notes reproduces assembled section text exactly; headings and labels derive from assembler intermediate\n')
print(f'SEED_PASS notes={len(checks)} sections={len(sections)}')
