#!/usr/bin/env python3
"""Write the existing intron/GC3 statistics as the technical table for the mapped supplementary table."""
import argparse
import csv
import json
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument('--pcwde', type=Path, required=True)
parser.add_argument('--blectin', type=Path, required=True)
parser.add_argument('--out', type=Path, required=True)
args = parser.parse_args()
pcwde = json.loads(args.pcwde.read_text())
blectin = json.loads(args.blectin.read_text())
if set(pcwde) != {'GH5_2', 'PL1_8', 'PL4_4'} or set(blectin) != {'all', 'chr6 full', 'chr3 reduced'}:
    raise SystemExit('GC3_TABLE_FAIL unexpected summary strata')
# Technical table: original schema and unrounded source values. The reader
# projection (headings and display rounding) is made by render_tables.py and
# checked against this table by check_intron_products.py.
fields = ['cohort', 'stratum', 'copies', 'species', 'median_host_gc3_percentile',
          'species_above_host_median', 'median_copy_host_distance_pp',
          'median_bacterial_host_distance_pp', 'paired_wilcoxon_p',
          'median_copy_signed_offset_pp', 'median_bacterial_signed_offset_pp']
rows = []
def nonnegative_distances(values):
    if float(values[6]) < 0 or float(values[7]) < 0:
        raise SystemExit('GC3_TABLE_FAIL absolute distance is negative')
    return values
for family in ('GH5_2', 'PL1_8', 'PL4_4'):
    s = pcwde[family]
    rows.append(nonnegative_distances([family, 'all', s['copies'], s['species'], s['median_species_host_percentile'],
                 s['species_above_host_median'], s['host_distance_pp'], s['bacterial_distance_pp'],
                 s['wilcoxon_bacteria_minus_host_distance_p'], s['median_copy_signed_offset_pp'],
                 s['median_bacterial_signed_offset_pp']]))
for stratum in ('all', 'chr6 full', 'chr3 reduced'):
    s = blectin[stratum]
    rows.append(nonnegative_distances(['B-lectin', stratum, s['copies'], s['species'],
                 s['species_median_host_percentile'], s['species_above_host_median'],
                 s['median_copy_distance_pp'], s['median_bacterial_distance_pp'],
                 s['wilcoxon_distance_p_two_sided'], s['median_copy_signed_offset_pp'],
                 s['median_bacterial_signed_offset_pp']]))
args.out.mkdir(parents=True, exist_ok=True)
with (Path(__file__).resolve().parents[2] / 'config/table_numbers.tsv').open(newline='') as handle:
    final_labels = {row['table_id']: row['label'] for row in csv.DictReader(handle, delimiter='\t')}
gc3_label = final_labels.get('gc3')
if not gc3_label or not gc3_label.startswith('S') or not gc3_label[1:].isdigit():
    raise SystemExit('GC3_TABLE_FAIL missing final GC3 label')
with (args.out / 'Table_S16.tsv').open('w', newline='') as handle:
    writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
    writer.writerow(fields)
    writer.writerows(rows)
(args.out / 'Table_S16_caption.txt').write_text('')
print(f'GC3_TABLE_PASS rows=6 mapped_label={gc3_label}')
