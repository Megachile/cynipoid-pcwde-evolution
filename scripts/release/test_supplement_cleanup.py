#!/usr/bin/env python3
"""Bounded label/removal and lineage-split fixture for the supplement cleanup."""
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'scripts/analysis/chromosome_gain'))
from ecami_score import gain_lineage  # noqa: E402


def rows(path):
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


tables = {r['table_id']: r['label'] for r in rows(ROOT / 'config/table_numbers.tsv')}
figures = {r['figure_id']: r['label'] for r in rows(ROOT / 'config/figure_numbers.tsv')}
sections = rows(ROOT / 'config/supplement_sections.tsv')
assert len(tables) == 10 and set(tables.values()) == {f'S{i}' for i in range(1, 11)}
assert 'di_sort_limits' not in tables and 'locus_conflict' not in tables
assert 'chromosome_gain_trees' not in figures
assert {r['label'] for r in rows(ROOT / 'config/figure_numbers.tsv') if r['label'].startswith('S')} == {f'S{i}' for i in range(1, 7)}
items = [item for r in sections for item in r['items'].split(',') if item]
assert len(items) == len(set(items)) == len(tables) + 6
assert {f'table:{key}' for key in tables} | {f'figure:{key}' for key, value in figures.items() if value.startswith('S')} == set(items)
assert next(r for r in sections if r['section_id'] == 'S-B')['items'].endswith('table:expression')
assert next(r for r in sections if r['section_id'] == 'S-F')['items'] == 'table:gc3,table:intron'
assert 'chromosome_gain_distribution' not in figures and 'figure:chromosome_gain_distribution' not in items
# Tables appear in numeric label order across the assembled sections.
shown = [int(tables[i.split(':', 1)[1]][1:]) for i in items if i.startswith('table:')]
assert shown == sorted(shown), shown
assert all(r['section_id'] != 'additional' for r in sections)
clade = {'A', 'B'}
events = [{'descendant_tips': 'A'}, {'descendant_tips': 'A,B'},
          {'descendant_tips': 'A,C'}, {'descendant_tips': 'C'}]
classes = [gain_lineage(event, clade) for event in events]
assert classes == ['cynipini_gains', 'cynipini_gains', 'elsewhere_gains', 'elsewhere_gains']
assert classes.count('cynipini_gains') + classes.count('elsewhere_gains') == len(events)
try:
    gain_lineage({'descendant_tips': ''}, clade)
except ValueError:
    pass
else:
    raise AssertionError('empty gain branch accepted')
print('SUPPLEMENT_CLEANUP_FIXTURE_PASS tables=10 figures=6 events=4')
