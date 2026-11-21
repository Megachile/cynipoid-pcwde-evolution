#!/usr/bin/env python3
"""Check current-draft pointer subjects and the figure-page boundary without rendering."""
import json
from pathlib import Path
import tempfile
from unittest.mock import patch
import fitz
import build_manuscript as bm

raw = bm.SOURCE.read_text()
pending = json.loads((bm.PAPER / 'authority/source_decisions.json').read_text())['exceptions']
if pending:
    try:
        bm.validate_frozen_authority()
    except RuntimeError as exc:
        assert 'pending decisions' in str(exc)
    else:
        raise AssertionError('pending source adoption allowed')
else:
    bm.validate_frozen_authority()
with tempfile.TemporaryDirectory() as tmp, patch('author_result_numbers.project', side_effect=lambda text, numbers, document: (text, [])):
    stats_root = Path(tmp) / 'introns/gc3'
    stats_root.mkdir(parents=True)
    (Path(tmp) / 'numbers.json').write_text(json.dumps({'supplement': {'outgroup_parent_proteins_total': 123}}))
    prose = bm.apply_authorized_supplement_references(raw, stats_root)
    try:
        bm.apply_authorized_supplement_references(raw.replace('(Table S6)', '(Table S5)'), stats_root)
    except RuntimeError:
        pass
    else:
        raise AssertionError('wrong-topic, in-range table pointer accepted')
assert 'Table S2)' not in prose
assert 'Table S2b)' in prose
assert 'supplementary text S7' in prose
assert '(Table S6)' in raw
# The author wording is preserved independently of mechanical pointer changes.
assert 'apparent chromosome translocation is best understood' not in prose
assert 'rate of interchromosomal duplication' in prose
assert bm.tex_escape('**2**') == r'\textbf{2}'
row = dict(content_id='fixture', number='1')
with tempfile.TemporaryDirectory() as tmp:
    aux = Path(tmp) / 'manuscript.aux'
    aux.write_text(r'\newlabel{fig:fixture}{{1}{3}}')
    pdf = Path(tmp) / 'fixture.pdf'
    with fitz.open() as document:
        for text in ('Text mentions Figure 1.', 'Second page', 'Figure 1. Caption', 'Fourth page'):
            page = document.new_page()
            page.insert_text((72, 72), text)
        document.save(pdf)
    assert bm.check_figure_pages(pdf, aux, [row])[0]['figure_page'] == 3
    aux.write_text(r'\newlabel{fig:fixture}{{1}{4}}')
    try:
        bm.check_figure_pages(pdf, aux, [row])
    except RuntimeError as exc:
        assert 'exceeds 2 pages' in str(exc)
    else:
        raise AssertionError('late figure accepted')
print('MAIN_TEXT_POINTERS_AND_FIGURE_PAGE_BOUNDARY_PASS')
