#!/usr/bin/env python3
"""Render the supplement prose intermediate from the repo-owned template and numbers.json.

authority/supplement/supplement_source.md is Adam's pinned Google Doc export with
each result number replaced by a numbers.json token (see supplement_template.py).
The build fails on an unknown token, any [CHECK note, a numeral that is neither a
token nor listed in config/supplement_literals.tsv, a stale literal row, or a
parameter that disagrees with its producer value. The pinned export itself is
still sha256-checked as the class U input the template was derived from."""
import csv
import hashlib
import json
import re
from pathlib import Path
import argparse
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
import supplement_template as st  # noqa: E402
from author_export import supplement_export

p = argparse.ArgumentParser()
p.add_argument('--repo', type=Path, required=True)
p.add_argument('--export', type=Path)
p.add_argument('--template', type=Path)
p.add_argument('--numbers', type=Path, required=True)
p.add_argument('--literals', type=Path)
p.add_argument('--report', type=Path)
p.add_argument('--out', type=Path, required=True)
a = p.parse_args()
repo = a.repo.resolve()
source = supplement_export(repo, a.export)
digest = hashlib.sha256(source.read_bytes()).hexdigest()
template = (a.template or repo / 'authority/supplement/supplement_source.md').resolve()
literals = st.read_literals(a.literals or repo / 'config/supplement_literals.tsv')
numbers = json.loads(a.numbers.read_text())
with (repo / 'config/supplement_sections.tsv').open() as f:
    expected = list(csv.DictReader(f, delimiter='\t'))
with (repo / 'config/figure_numbers.tsv').open() as f:
    figures = {r['figure_id']: r['label'] for r in csv.DictReader(f, delimiter='\t')}
with (repo / 'config/table_numbers.tsv').open() as f:
    tables = {r['table_id']: r['label'] for r in csv.DictReader(f, delimiter='\t')}
pending = repo / "authority/source_refresh_pending.json"
if pending.exists() and json.loads(pending.read_text())["exceptions"]:
    raise SystemExit("author source refresh has pending decisions; full build prohibited")
# Item placement is derived independently from the immutable author export.
source_sections = re.split(r"^## S[0-9]+\. .*?$", source.read_text().replace('**Table S2.**', '**Table S2a.** **Table S2b.**'), flags=re.M)[1:]
if len(source_sections) != len(expected):
    raise SystemExit("author export section coverage mismatch")
seen = set()
for section, body in zip(expected, source_sections):
    labels = re.findall(r"\*\*((?:Figure|Table) S[0-9]+[ab]?)\.\*\*", body)
    for label in labels:
        if label in seen:
            raise SystemExit(f"duplicate author item label: {label}")
        seen.add(label)
    wanted = []
    for item in section['items'].split(','):
        if not item:
            continue
        kind, stable_id = item.split(':')
        label = (figures if kind == 'figure' else tables)[stable_id]
        wanted.extend(['Table S2a', 'Table S2b'] if stable_id == 'expression' else
                      [('Figure ' if kind == 'figure' else 'Table ') + label])
    if labels != wanted:
        raise SystemExit(f"author item placement differs for {section['section_id']}: {labels} != {wanted}")
raw = template.read_text()
for context in (
    'assigned to CynOG (S5) by hmmsearch',
    'membership in the CynOG (S5) containing',
    'similar median copy number (Table S7).',
):
    if raw.count(context) != 1:
        raise SystemExit(f'supplement subject pointer missing or ambiguous: {context}')
if 'Table S7 and Figure S7' in raw:
    raise SystemExit('retired Figure S7 pointer in current supplement')
lines = raw.splitlines()
starts = [(i, line[3:]) for i, line in enumerate(lines) if line.startswith('## ')]
titles = [s['number'] + '. ' + s['heading'] for s in expected]
got = [title for _, title in starts]
source_titles = titles
for title in source_titles:
    if got.count(title) != 1:
        raise SystemExit(f'missing or duplicate supplement heading: {title}')
if got != source_titles:
    raise SystemExit('supplement headings reordered or unexpected: ' + repr(got))
source_labels = {}
for section in expected:
    labels = []
    for item in section['items'].split(','):
        if not item:
            continue
        kind, stable_id = item.split(':')
        labels.extend(['Table S2a', 'Table S2b'] if stable_id == 'expression' else
                      [('Figure ' if kind == 'figure' else 'Table ') + (figures if kind == 'figure' else tables)[stable_id]])
    source_labels[section['section_id']] = labels

sections = []
failures = []
token_rows = []
unknown = {r['section'] for r in literals} - {s['section_id'] for s in expected}
if unknown:
    failures.append(f'literal rows for unknown sections: {sorted(unknown)}')
for j, section in enumerate(expected):
    body = '\n'.join(lines[starts[j][0]+1: starts[j+1][0] if j+1 < len(starts) else len(lines)]).strip()
    blocks = [b.strip().replace('  \n', '\n') for b in re.split(r'\n\s*\n', body) if b.strip()]
    if section['author_slot']:
        if not blocks or blocks.pop(0).replace('\\[', '[').replace('\\]', ']') != section['author_slot']:
            raise SystemExit(f'missing author slot: {titles[j]}')
    items = [x for x in section['items'].split(',') if x]
    labels = [f'**{label}.** [Item label; source artwork/table remains in the assembled supplement]'
              for label in source_labels[section['section_id']]]
    def decode_doc_notes(block):
        return re.sub(r'\\\[CHECK .*?\\\]',
                      lambda m: re.sub(r'\\([\[\]#_+])', r'\1', m.group()),
                      block)
    normalized = [decode_doc_notes(b) for b in blocks]
    normalized = [b.replace('\\[', '[').replace('\\]', ']') if b.startswith('**Figure ') or b.startswith('**Table ') else b for b in normalized]
    found = [b for b in normalized if b.startswith('**Figure ') or b.startswith('**Table ')]
    if found != labels:
        raise SystemExit(f'supplement item labels mismatch: {titles[j]}')
    paragraphs = [b for b in normalized if b not in labels]
    failures += [f'{titles[j]}: [CHECK note in prose' for b in paragraphs if '[CHECK' in b.replace('\\[', '[')]
    gate_errors, _ = st.check_section(section['section_id'], paragraphs, literals, numbers, repo)
    failures += gate_errors
    rendered = []
    for b in paragraphs:
        text, used, errors = st.render(b, numbers)
        failures += [f"{titles[j]}: {e}" for e in errors]
        token_rows += [(section['section_id'], k, f, v) for k, f, v in used]
        rendered.append(text)
    sections.append(dict(section_id=section['section_id'], heading=titles[j], author_slot=section['author_slot'], paragraphs=rendered, items=items))
# Exact paragraph adoption is checked separately from numeric validity.
from author_outgroup_counts import project
from author_result_numbers import project as project_results
projected_author, _ = project_results(source.read_text(), numbers, 'supplement')
projected_author, _ = project(projected_author, numbers, 'supplement')
author_paragraphs = [p.strip().replace('  \n', '\n') for p in re.split(r'\n\s*\n', projected_author)
                     if p.strip() and not p.strip().startswith(('##', '**Figure ', '**Table '))]
rendered_paragraphs = [p for section in sections for p in section['paragraphs']]
if rendered_paragraphs != author_paragraphs:
    failures.append('rendered supplement prose differs from the pinned author export')
if a.report:
    a.report.parent.mkdir(parents=True, exist_ok=True)
    a.report.write_text(''.join(f + '\n' for f in failures))
if failures:
    sys.stderr.write(''.join(f'SUPPLEMENT_NUMBERS_FAIL {f}\n' for f in failures))
    raise SystemExit(f'supplement prose gate failed: {len(failures)} problem(s)')
a.out.parent.mkdir(parents=True, exist_ok=True)
with (a.out.parent / 'supplement_tokens.tsv').open('w', newline='') as handle:
    w = csv.writer(handle, delimiter='\t', lineterminator='\n')
    w.writerow(['section', 'key', 'format', 'rendered'])
    w.writerows(token_rows)
a.out.write_text(json.dumps(sections, ensure_ascii=False, indent=2) + '\n')
print(f'SUPPLEMENT_SECTIONS_PASS sections={len(sections)} tokens={len(token_rows)} literals={len(literals)} input_sha256={digest} template_sha256={hashlib.sha256(template.read_bytes()).hexdigest()} method_sha256={hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}')
