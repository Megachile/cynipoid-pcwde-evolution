#!/usr/bin/env python3
"""Score retained chromosome gains for accepted PCWDEs and native eCAMI groups.

Parsimony, pruning, node IDs and taxon sets come from the merged scorer.
"""
import argparse
import csv
import hashlib
import json
import statistics
from collections import Counter, defaultdict
from pathlib import Path

from score import (CHROMS, anchor_clade, descendants, score_character,
                   taxon_groups, write)
from tree import fitch, parse_newick, prune

FOCAL = ('GH5_2', 'PL1_8', 'PL4_4', 'B-lectin')
MODES = ('deltran', 'acctran', 'sankoff')


def rows(path):
    with path.open(newline='') as handle:
        yield from csv.DictReader(handle, delimiter='\t')


def empty_tip():
    return {'placed': set(), 'copies': set()}


def gain_lineage(event, cynipini):
    members = set(filter(None, event['descendant_tips'].split(',')))
    if not members:
        raise ValueError('gain branch has no descendant tips')
    return 'cynipini_gains' if members <= cynipini else 'elsewhere_gains'


def fixed_component_groups(components, mapping):
    sizes = defaultdict(lambda: defaultdict(set))
    for (species, group), genes in components.items():
        family = mapping.get(group)
        if family in FOCAL[:3]:
            sizes[family][group].update((species, gene) for gene in genes)
    selected = {family: min(groups, key=lambda group: (-len(groups[group]), group))
                for family, groups in sizes.items()}
    return sizes, selected


def main():
    p = argparse.ArgumentParser()
    for name in ('states', 'components', 'group-mapping', 'protein-genes',
                 'native-counts', 'crosswalk', 'primary-ranks', 'focal',
                 'panel', 'tree', 'figure-code', 'expected-focal', 'out'):
        p.add_argument('--' + name, type=Path, required=True)
    p.add_argument('--historical-regression', action='store_true')
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    panel = {r['species'].replace(' ', '_'): r['genome'] for r in rows(a.panel)
             if r['species'] != 'Telenomus remus'}
    expected_panel = 40 if a.historical_regression else 48
    if len(panel) != expected_panel:
        raise SystemExit(f'expected {expected_panel} cynipoid panel genomes')
    genome_species = {genome: species for species, genome in panel.items()}
    taxon = taxon_groups(a.figure_code)
    lineage = {species: side for side, members in taxon.items()
               for species in members}
    if set(panel) - set(lineage):
        raise SystemExit('panel species absent from Figure 2 taxonomy')
    tree = prune(parse_newick(a.tree.read_text()), set(panel))
    if len(descendants(tree)) != expected_panel:
        raise SystemExit('species tree/panel mismatch')
    cynipini = anchor_clade(tree, {'Dryocosmus_kuriphilus', 'Acraspis_erinacei'})
    if a.historical_regression and len(cynipini) != 10:
        raise SystemExit('unexpected Cynipini clade')

    # One placed state per gene; duplicated scan protein records remain one gene.
    gene_states = defaultdict(set)
    gene_proteins = defaultdict(set)
    state_rows = 0
    for r in rows(a.states):
        state_rows += 1
        key = r['species'], r['gene_id']
        gene_proteins[key].add(r['protein_id'])
        if r['state'] in CHROMS:
            gene_states[key].add(r['state'])
    expected_state_rows = sum(1 for r in rows(a.protein_genes) if r['species'] in panel)
    if state_rows != expected_state_rows or any(len(v) > 1 for v in gene_states.values()):
        raise SystemExit('gene collapse produced missing or multi-chromosome states')
    scan_gene = {(r['species'], r['protein_id']): r['gene_id']
                 for r in rows(a.protein_genes) if r['species'] in panel}
    accepted = defaultdict(set)
    for r in rows(a.crosswalk):
        if r['species'] in panel:
            accepted[r['species']].update(filter(None, r['scan_gene_ids'].split(',')))
    blectin_proteins = {(r['species'].replace(' ', '_'), protein)
                        for r in rows(a.focal) if r['family'] == 'B-lectin'
                        for protein in r['protein_ids'].split(',')}
    blectin_genes = {(sp, scan_gene[sp, protein]) for sp, protein in blectin_proteins
                     if (sp, protein) in scan_gene}
    accepted.update({sp: accepted[sp] | {gene for xsp, gene in blectin_genes if xsp == sp}
                     for sp in panel})
    mapping = {r['ecami_group']: r['rank_unit'] for r in rows(a.group_mapping)}
    components = defaultdict(set)
    blectin_groups = set()
    for r in rows(a.components):
        sp, gene, group = r['species'], r['gene_id'], r['ecami_group']
        if sp not in panel:
            continue
        components[sp, group].add(gene)
        if (sp, gene) in blectin_genes:
            blectin_groups.add(group)
    excluded_groups = {g for g, unit in mapping.items() if unit in FOCAL} | blectin_groups
    native = defaultdict(lambda: defaultdict(empty_tip))
    native_counts = Counter()
    for (sp, group), genes in components.items():
        if group in excluded_groups:
            continue
        genome = panel[sp]
        for gene in genes - accepted[sp]:
            tip = native[group][genome]
            tip['copies'].add(gene)
            tip['placed'].update(gene_states[sp, gene])
            native_counts[sp, group] += 1
    expected_native = {(r['species'], r['ecami_group']): int(r['gene_count'])
                       for r in rows(a.native_counts) if r['species'] in panel
                       and r['ecami_group'] not in excluded_groups}
    if native_counts != expected_native:
        diff = set(native_counts.items()) ^ set(expected_native.items())
        raise SystemExit(f'native component collapse differs from rank producer: {list(diff)[:5]}')

    focal = defaultdict(lambda: defaultdict(empty_tip))
    for r in rows(a.focal):
        sp = r['species'].replace(' ', '_')
        if sp not in panel or r['genome'] not in genome_species:
            continue
        tip = focal[r['family']][r['genome']]
        tip['copies'].add(r['locus_protein_id'])
        if r['bk_chromosome'] in CHROMS:
            tip['placed'].add(r['bk_chromosome'])
    primary_groups = {(r['species'], r['unit']): r['focal_ecami_group']
                      for r in rows(a.primary_ranks) if r['route'] == 'largest_group'
                      and r['species'] in panel and r['unit'] in FOCAL[:3]}
    component_sizes, fixed_groups = fixed_component_groups(components, mapping)
    if set(fixed_groups) != set(FOCAL[:3]):
        raise SystemExit('missing focal component group')
    write(a.out / 'focal_component_sizes.tsv', ['family', 'ecami_group', 'gene_count', 'selected'],
          [dict(family=family, ecami_group=group, gene_count=len(genes),
                selected=int(group == fixed_groups[family]))
           for family, groups in sorted(component_sizes.items()) for group, genes in sorted(groups.items())])
    focal_scan = defaultdict(lambda: defaultdict(empty_tip))
    for sp in panel:
        for family, group in fixed_groups.items():
            genes = components.get((sp, group), set())
            if not genes:
                continue
            tip = focal_scan[family][panel[sp]]
            for gene in genes:
                tip['copies'].add(gene)
                tip['placed'].update(gene_states[sp, gene])
    # Historical genome-wise maximum is retained only as a diagnostic trajectory.
    focal_scan_max = defaultdict(lambda: defaultdict(empty_tip))
    for (sp, family), group in primary_groups.items():
        genes = components.get((sp, group), set())
        if not genes:
            continue
        tip = focal_scan_max[family][panel[sp]]
        for gene in genes:
            tip['copies'].add(gene)
            tip['placed'].update(gene_states[sp, gene])

    eligibility = []
    selected = {}
    for group, tips in sorted(native.items()):
        placed = {g: v for g, v in tips.items() if v['placed']}
        sides = {lineage[genome_species[g]] for g in placed}
        if len(placed) < 5:
            reason = 'lt5_placed_tips'
        elif 'Cynipidae' not in sides or not ({'Diplolepididae', 'Figitidae'} & sides):
            reason = 'not_crown_present'
        else:
            reason = 'eligible'
            selected['NATIVE:' + group] = tips
        eligibility.append(dict(group=group, reason=reason, member_genomes=len(tips),
                                placed_tips=len(placed), placed_side_sets=','.join(sorted(sides))))
    eligible_count = sum(r['reason'] == 'eligible' for r in eligibility)
    write(a.out / 'background_eligibility.tsv',
          ['group', 'reason', 'member_genomes', 'placed_tips', 'placed_side_sets'], eligibility)
    if eligible_count < 10:
        raise SystemExit(f'only {eligible_count} native groups eligible; stop before comparison')
    for family in FOCAL:
        selected['FOCAL:' + family] = focal[family]
    for family in FOCAL[:3]:
        selected['SCAN:' + family] = focal_scan[family]
        selected['SCAN_MAX_DIAGNOSTIC:' + family] = focal_scan_max[family]
    scan_eligibility = []
    for family in FOCAL[:3]:
        tips = focal_scan[family]
        placed = {g for g, v in tips.items() if v['placed']}
        sides = {lineage[genome_species[g]] for g in placed}
        eligible = len(placed) >= 5 and 'Cynipidae' in sides and bool({'Diplolepididae', 'Figitidae'} & sides)
        scan_eligibility.append(dict(family=family, ecami_group=fixed_groups[family],
                                     member_genomes=len(tips), placed_tips=len(placed),
                                     crown_eligible=int(eligible), placed_side_sets=','.join(sorted(sides))))
    write(a.out / 'focal_scan_eligibility.tsv',
          ['family', 'ecami_group', 'member_genomes', 'placed_tips', 'crown_eligible', 'placed_side_sets'],
          scan_eligibility)

    records, changes, tiprows = [], [], []
    for ident, tips in sorted(selected.items()):
        placed = {g: v for g, v in tips.items() if v['placed']}
        if not placed:
            continue
        subtree = prune(tree, {genome_species[g] for g in placed})
        root_sets = {}
        totals = {mode: Counter() for mode in MODES}
        for genome, tip in tips.items():
            tiprows.append(dict(group=ident, genome=genome, species=genome_species[genome],
                                chromosomes=','.join(sorted(tip['placed'])),
                                placed=int(bool(tip['placed'])), copies=len(tip['copies'])))
        for chrom in CHROMS:
            ones = {genome_species[g] for g, v in placed.items() if chrom in v['placed']}
            possible, fitch_cost = fitch(subtree, ones)
            root_sets[chrom] = ''.join(map(str, sorted(possible[id(subtree)])))
            for mode in MODES:
                root_state, cost, gains, losses = score_character(subtree, ones, mode)
                if mode != 'sankoff' and cost != fitch_cost:
                    raise SystemExit(f'Fitch mismatch: {ident} {chrom} {mode}')
                totals[mode]['gains'] += len(gains)
                totals[mode]['losses'] += len(losses)
                for event in gains:
                    totals[mode][gain_lineage(event, cynipini)] += 1
                if ident.startswith(('FOCAL:', 'SCAN:')):
                    for kind, events in (('gain', gains), ('loss', losses)):
                        for event in events:
                            changes.append(dict(group=ident, mode=mode, chromosome=chrom,
                                                change=kind, root_state=root_state, **event))
        copies = [len(v['copies']) for v in tips.values()]
        rec = dict(group=ident, kind=ident.split(':', 1)[0].lower(),
                   placed_tips=len(placed), member_genomes=len(tips),
                   median_copies=statistics.median(copies),
                   single_copy_core=int(ident.startswith('NATIVE:') and all(n == 1 for n in copies)),
                   root_present=','.join(ch for ch, s in root_sets.items() if s == '1'),
                   root_ambiguous=','.join(ch for ch, s in root_sets.items() if s == '01'),
                   cynipini_placed_tips=sum(genome_species[g] in cynipini for g in placed),
                   elsewhere_placed_tips=sum(genome_species[g] not in cynipini for g in placed),
                   cynipini_gains=totals['deltran']['cynipini_gains'],
                   elsewhere_gains=totals['deltran']['elsewhere_gains'])
        for mode in MODES:
            rec[mode + '_gains'] = totals[mode]['gains']
            rec[mode + '_cynipini_gains'] = totals[mode]['cynipini_gains']
            rec[mode + '_elsewhere_gains'] = totals[mode]['elsewhere_gains']
            if totals[mode]['gains'] != totals[mode]['cynipini_gains'] + totals[mode]['elsewhere_gains']:
                raise SystemExit(f'lineage gain split mismatch: {ident} {mode}')
            rec[mode + '_losses'] = totals[mode]['losses']
            rec[mode + '_rate'] = round(totals[mode]['gains'] / len(placed), 8)
        records.append(rec)
    columns = ['group', 'kind', 'placed_tips', 'member_genomes', 'median_copies',
               'single_copy_core', 'root_present', 'root_ambiguous',
               'cynipini_placed_tips', 'elsewhere_placed_tips',
               'cynipini_gains', 'elsewhere_gains']
    for mode in MODES:
        columns += [mode + '_gains', mode + '_cynipini_gains', mode + '_elsewhere_gains', mode + '_losses', mode + '_rate']
    write(a.out / 'ecami_scores.tsv', columns, records)
    write(a.out / 'ecami_tip_states.tsv',
          ['group', 'genome', 'species', 'chromosomes', 'placed', 'copies'], tiprows)
    write(a.out / 'ecami_node_changes.tsv',
          ['group', 'mode', 'chromosome', 'change', 'root_state', 'parent',
           'child', 'depth', 'descendant_tips'], changes)
    focal_records = {r['group']: r for r in records if r['kind'] == 'focal'}
    if a.historical_regression:
        expected = {r['family']: r for r in rows(a.expected_focal)}
        if set(expected) != set(FOCAL):
            raise SystemExit('focal parity snapshot lacks a family')
        for family in FOCAL:
            current = focal_records['FOCAL:' + family]
            for key, value in expected[family].items():
                if key != 'family' and str(current[key]) != value:
                    raise SystemExit(f'focal gain parity failed: {family} {key}: {current[key]} != {value}')
        for family, chromosome in (('PL1_8', 'chr8'), ('PL4_4', 'chr7')):
            if chromosome not in focal_records['FOCAL:' + family]['root_present'].split(','):
                raise SystemExit(f'{family} root control failed')
    qc = dict(panel_genomes=len(panel), eligible_native_groups=eligible_count,
              excluded_groups=sorted(excluded_groups), blectin_mapped_genes=len(blectin_genes),
              scored_native_groups=sum(r['kind'] == 'native' for r in records),
              scored_focal_families=sum(r['kind'] == 'focal' for r in records),
              single_copy_native_groups=sum(r['single_copy_core'] for r in records))
    (a.out / 'ecami_score_qc.json').write_text(json.dumps(qc, indent=2) + '\n')
    def digest(path):
        return hashlib.sha256(path.read_bytes()).hexdigest()
    input_names = ('states', 'components', 'group_mapping', 'protein_genes',
                   'native_counts', 'crosswalk', 'primary_ranks', 'focal',
                   'panel', 'tree', 'figure_code', 'expected_focal')
    (a.out / 'ecami_score_sha256.json').write_text(json.dumps({
        'inputs': {name: digest(getattr(a, name)) for name in input_names},
        'outputs': {name: digest(a.out / name) for name in
                    ('ecami_scores.tsv', 'ecami_tip_states.tsv', 'ecami_node_changes.tsv', 'ecami_score_qc.json')},
    }, indent=2, sort_keys=True) + '\n')
    print(json.dumps(qc, sort_keys=True))


if __name__ == '__main__':
    main()
