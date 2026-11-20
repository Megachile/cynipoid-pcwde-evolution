#!/usr/bin/env python3
"""Reconstruct the 99-row S1 resource inventory from the pinned source panel.

The panel records inclusion, accession, preparation and most metrics. Eleven
BUSCO cells come from a frozen catalog export. Display wording follows the
accepted paper presentation; no inclusion or measurement is recalculated.
"""
import argparse
import csv
from pathlib import Path

FIELDS = ('resource_class', 'species_display', 'source_accession',
          'resource_preparation', 'panel_status_exclusion_reason',
          'assembly_n50_bp', 'assembly_busco_pct', 'proteome_busco_pct')
CLASS = {'genome': 'Genome', 'transcriptome': 'Transcriptome-derived proteome',
         'reference_proteome': 'Other reference proteome'}
PREP = {'NCBI / genome assembly': 'Genome assembly',
        'ENA (assembled) / genome assembly': 'Genome assembly',
        'NCBI GenBank / genome assembly': 'Genome assembly',
        'USDA / genome assembly': 'USDA PBARC genome assembly',
        'RNA-seq transcriptome assembly / Trinity / EviGene / TransDecoder': 'This study: Trinity / EviGene / TransDecoder',
        'NCBI reference proteome / reference proteome/genome assembly': 'NCBI reference proteome',
        'i5K OGS1.0 proteome / reference proteome/genome assembly': 'i5K OGS1.0 proteome'}
NAMES = {'Ceroptres sp. (D. ignotum inquiline)': 'Ceroptres sp. ex Druon ignotum',
         'Synergus sp. (B. kinseyi inquiline)': 'Synergus sp. ex Belonocnema kinseyi',
         'Synergus sp. (D. ignotum inquiline)': 'Synergus sp. ex Druon ignotum',
         'Synergus sp. (D. quercuslanigerum inquiline)': 'Synergus sp. ex Druon quercuslanigerum'}

def read(path):
    with Path(path).open(newline='', encoding='utf-8') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))

def status(row, retained_accessions):
    state, reason = row['status'], row['exclusion_reason']
    if reason.startswith('Below 40% BUSCO completeness:'):
        suffix = '; higher-quality conspecific assembly retained' if 'duplicate species' in reason else ''
        return 'Excluded — <40% complete BUSCOs' + suffix
    if 'replaced by Euceroptres_primus_USDA' in reason:
        return 'Excluded — superseded by retained USDA PBARC assembly'
    if 'retained Leptopilina sampling' in reason:
        return 'Excluded — redundant; retained Leptopilina sampling was sufficient'
    if 'redundant (' in reason:
        return 'Excluded — redundant; higher-quality alternatives retained'
    if 'redundant for paper panel' in reason:
        return 'Excluded — redundant; another assembly retained'
    if reason.startswith('superseded by current ENA assembly'):
        matches = retained_accessions.get(row['species'], ())
        if len(matches) != 1 or not matches[0].startswith('GCA_'):
            raise ValueError(f'No unique retained GenBank replacement for {row["species"]}: {matches}')
        return f'Excluded — superseded by GenBank assembly {matches[0]}'
    if 'assembly too fragmented' in reason:
        return 'Excluded — assembly too fragmented (N50 2 kb; approximately 50% complete BUSCOs missing; inflated gene count)'
    if state == 'retained': return 'Retained'
    if state == 'transcriptome input': return 'Transcriptome input'
    if state == 'reference outgroup': return 'Reference outgroup'
    raise ValueError(f'Unmapped panel status: {state}: {reason}')

def accession(row):
    source, value = row['source_and_assembly'], row['accession_or_project']
    prefix = {'NCBI / genome assembly': 'NCBI: ',
              'ENA (assembled) / genome assembly': 'ENA reads: ',
              'NCBI GenBank / genome assembly': 'NCBI GenBank: ',
              'RNA-seq transcriptome assembly / Trinity / EviGene / TransDecoder': 'RNA-seq study/accession: ',
              'NCBI reference proteome / reference proteome/genome assembly': 'NCBI: ',
              'i5K OGS1.0 proteome / reference proteome/genome assembly': 'i5K: '}
    if source == 'USDA / genome assembly':
        assert value == 'USDA PBARC unpublished assembly'
        return 'USDA PBARC; reads: PRJNA1518226; assembly accession pending'
    if source == 'ENA (assembled) / genome assembly' and value == 'PRJNA647791':
        return 'SRA reads: PRJNA647791 (assembled in this study)'
    label = prefix[source] + value
    if source == 'ENA (assembled) / genome assembly':
        label += ' (assembled in this study)'
    return label

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--source-panel', type=Path, required=True)
    ap.add_argument('--catalog-pin', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    a = ap.parse_args()
    panel, catalog = read(a.source_panel), read(a.catalog_pin)
    if len(panel) != 99 or len(catalog) != 11:
        raise SystemExit('Unexpected panel or frozen catalog cohort')
    # The displayed Periclistus pair is sorted by displayed name; the source
    # panel puts the JH2016 assembly before P. caninae.
    panel[50], panel[51] = panel[51], panel[50]
    retained_accessions = {}
    for source in panel:
        if source['status'] == 'retained':
            retained_accessions.setdefault(source['species'], []).append(source['accession_or_project'])
    pins = {int(x['s1_row']): x for x in catalog}
    if len(pins) != 11: raise SystemExit('Duplicate frozen catalog row')
    rows = []
    for i, source in enumerate(panel, 1):
        old_accession = source['accession_or_project']
        score_a, score_p = source['assembly_busco_C'], source['proteome_busco_C']
        if i in pins:
            pin = pins[i]
            if pin['accession'] != old_accession or pin['mode'] not in ('assembly', 'proteome'):
                raise SystemExit(f'Frozen catalog identity/mode mismatch at row {i}')
            if pin['mode'] == 'assembly':
                if score_a: raise SystemExit(f'Catalog assembly value would overwrite panel at row {i}')
                score_a = pin['score']
            else:
                if score_p: raise SystemExit(f'Catalog proteome value would overwrite panel at row {i}')
                score_p = pin['score']
        rows.append(dict(zip(FIELDS, (CLASS[source['source_type']],
                  NAMES.get(source['species'], source['species']), accession(source),
                  PREP[source['source_and_assembly']], status(source, retained_accessions),
                  source['n50_bp'], score_a, score_p))))
    if any(not r['species_display'] or not r['source_accession'] for r in rows):
        raise SystemExit('Missing display identity')
    a.out.mkdir(parents=True, exist_ok=True)
    with (a.out / 'Table_S1.tsv').open('w', newline='', encoding='utf-8') as fh:
        writer = csv.DictWriter(fh, delimiter='\t', fieldnames=FIELDS, lineterminator='\n')
        writer.writeheader(); writer.writerows(rows)
    print(f'TABLE_S1_INVENTORY_PASS rows={len(rows)} catalog_cells={len(pins)}')

if __name__ == '__main__': main()
