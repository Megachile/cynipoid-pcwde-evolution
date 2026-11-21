#!/usr/bin/env python3
"""Build gene-resolved CAZy counts from the pinned full competitive scans."""

import argparse
import csv
import hashlib
from collections import defaultdict
from pathlib import Path

from cazy_gene_rule import map_genes
from cazy_data_root import resolve_input
from cazy_subfamily_rule import assign


def rows(path):
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def write(path, fields, values):
    with path.open('w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(fields)
        writer.writerows(values)


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--manifest', type=Path, required=True)
    p.add_argument('--data-root', '--source-root', dest='source_root', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--fig1-counts', type=Path, required=True,
                   help='pinned Figure 1 presence_counts input directory')
    args = p.parse_args()
    manifest = rows(args.manifest)
    if len(manifest) != 64 or len({r['species'] for r in manifest}) != 64:
        raise SystemExit('scan manifest must contain 64 unique Figure 1 species')
    figure_counts = {r['species']: int(r['GH1_e72'])
                     for r in rows(args.fig1_counts / 'gh1e72_gene_counts.tsv')}
    figure_counts.update({r['species']: int(r['GH1_e72'])
                          for r in rows(args.fig1_counts / 'internal_gh1e72_gene_counts.tsv')})
    source_to_species = {r['source_key']: r['species'] for r in manifest}
    for r in rows(args.fig1_counts / 'gh1e72_cynipoid_counts.tsv'):
        if r['gh1e72_locus_count']:
            figure_counts[source_to_species[r['accession']]] = int(r['gh1e72_locus_count'])
    if set(figure_counts) != {r['species'] for r in manifest}:
        raise ValueError('Figure 1 GH1_e72 roster differs from 64 pinned scans')
    reductions = Path(__file__).resolve().parents[2] / 'config/public_reference_reductions.tsv'
    pin_only = {r['species'] for r in rows(reductions) if r['removed_kind'] == 'reference_proteome'}
    groups = {}
    gene_rows, count_rows, component_rows, component_gene_rows = [], [], [], []
    gh1_rows, multi_rows, protein_rows, competing_rows = [], [], [], []
    for item in manifest:
        species = item['species']
        scan = resolve_input(args.source_root, item['scan_path'])
        if not item['proteome_path'] and species not in pin_only:
            raise ValueError(f'{species}: blank proteome path is allowed only for listed public references')
        proteome = resolve_input(args.source_root, item['proteome_path']) if item['proteome_path'] else None
        if sha(scan) != item['scan_sha256'] or (proteome and sha(proteome) != item['proteome_sha256']):
            raise ValueError(f'{species}: pinned scan or proteome hash mismatch')
        hits = rows(scan)
        proteins = {h['Target Name'] for h in hits}
        source = resolve_input(args.source_root, item['mapping_source']) if item['mapping_source'] else None
        if source and sha(source) != item['mapping_sha256']:
            raise ValueError(f'{species}: mapping source hash mismatch: {source}')
        genes = map_genes(item['mapping_method'], source, proteins)
        if set(genes) != proteins:
            raise ValueError(f'{species}: incomplete gene mapping')
        protein_rows.extend((species, protein, gene) for protein, gene in genes.items())
        units = defaultdict(set)
        components = defaultdict(set)
        gh1 = set()
        gene_family_hits = defaultdict(list)
        for hit in hits:
            group = hit['Subfam Name']
            composition = hit['Subfam Composition']
            unit, purity = assign(group, composition)
            record = (unit, round(purity, 6), composition)
            if group in groups and groups[group] != record:
                raise ValueError(f'{group}: inconsistent composition across scans')
            groups[group] = record
            gene = genes[hit['Target Name']]
            family = unit.split('_', 1)[0]
            gene_family_hits[(gene, family)].append((unit, float(hit['i-Evalue']), group))
            components[group].add(gene)
            if group == 'GH1_e72':
                gh1.add(gene)
        for (gene, family), candidates in gene_family_hits.items():
            named = [x for x in candidates if x[0] != family]
            if named:
                chosen = min(named, key=lambda x: (x[1], x[0], x[2]))[0]
                distinct = sorted({x[0] for x in named})
                if len(distinct) > 1:
                    competing_rows.append((species, gene, family, ';'.join(distinct), chosen))
            else:
                chosen = family
            units[chosen].add(gene)
        if len(gh1) != int(item['expected_GH1_e72']) or len(gh1) != figure_counts[species]:
            raise ValueError(f'{species}: GH1_e72 control {len(gh1)} != manifest '
                             f'{item["expected_GH1_e72"]} or Figure 1 {figure_counts[species]}')
        gh1_rows.append((species, len(gh1), figure_counts[species], 'match'))
        by_gene = defaultdict(set)
        for unit, members in units.items():
            count_rows.append((species, unit, len(members)))
            for gene in members:
                gene_rows.append((species, gene, unit))
                by_gene[gene].add(unit.split('_', 1)[0])
        for family, members in units.items():
            if '_' not in family:
                for unit, named_members in units.items():
                    if unit.startswith(family + '_') and members & named_members:
                        raise ValueError(f'{species}: gene in {family} and {unit}')
        for group, members in components.items():
            component_rows.append((species, group, len(members)))
            component_gene_rows.extend((species, gene, group) for gene in members)
        multi_rows.append((species, sum(len(families) > 1 for families in by_gene.values()), len(by_gene)))
    args.out.mkdir(parents=True, exist_ok=False)
    write(args.out / 'group_mapping.tsv', ('ecami_group','rank_unit','purity','composition'),
          ((group, *record) for group, record in sorted(groups.items())))
    write(args.out / 'gene_units.tsv', ('species','gene_id','rank_unit'), sorted(gene_rows))
    write(args.out / 'protein_genes.tsv', ('species','protein_id','gene_id'), sorted(protein_rows))
    write(args.out / 'unit_counts.tsv', ('species','rank_unit','gene_count'), sorted(count_rows))
    write(args.out / 'component_counts.tsv', ('species','ecami_group','gene_count'), sorted(component_rows))
    write(args.out / 'component_gene_units.tsv', ('species','gene_id','ecami_group'), sorted(component_gene_rows))
    write(args.out / 'gh1_control.tsv', ('species','derived_GH1_e72','expected_GH1_e72','status'), sorted(gh1_rows))
    write(args.out / 'derived_gh1e72_gene_counts.tsv', ('species', 'GH1_e72'),
          sorted((species, count) for species, count, _, _ in gh1_rows))
    write(args.out / 'multi_family_genes.tsv', ('species','genes_in_multiple_families','total_hit_genes'), sorted(multi_rows))
    write(args.out / 'competing_subfamilies.tsv', ('species','gene_id','family','named_subfamilies','selected_subfamily'), sorted(competing_rows))
    print(f'CAZY_SCAN_COUNTS_PASS genomes={len(manifest)} groups={len(groups)} gene_units={len(gene_rows)}')


if __name__ == '__main__':
    main()
