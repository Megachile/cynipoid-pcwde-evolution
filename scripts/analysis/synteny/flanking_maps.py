#!/usr/bin/env python3
"""Build 12-gene flanks for assigned corrected synteny copies."""
import argparse
import csv
import hashlib
import importlib.util
from collections import defaultdict
from pathlib import Path

FLANKS = [f'U{i}' for i in range(12, 0, -1)] + [f'D{i}' for i in range(1, 13)]
FIELDS = ['family', 'block_id', 'accepted_locus_id', 'species', 'scaffold', *FLANKS]

def read(path):
    with Path(path).open(newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))

def checked(path, digest):
    actual = hashlib.sha256(Path(path).read_bytes()).hexdigest()
    if actual != digest:
        raise SystemExit(f'SHA256 mismatch for {path}: {actual}')

def produce(order_path, affected_path, member_path, old_path, labels_path, out_path, report_path):
    copy_table = Path(__file__).resolve().parents[2] / 'copy_table.py'
    spec = importlib.util.spec_from_file_location('paper_copy_table', copy_table)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    labels = {r['accession']: module.PAPER_NAMES.get(r['species'], r['species']).replace(' ', '_')
              for r in read(labels_path)}
    if len(labels) != len(read(labels_path)):
        raise SystemExit('Duplicate genome in species-label authority')
    order = defaultdict(list)
    for row in read(order_path):
        order[row['genome'], row['chrom']].append(row)
    # Reproduce 01_build_synteny.inject_reconstructed_order: affected scaffolds
    # replace their old rows, retaining CynOG labels by (genome, protein ID).
    affected = defaultdict(list)
    for row in read(affected_path):
        affected[row['genome'], row['scaffold']].append(row)
    gene_labels = {(g, gene['protein_id']): gene['cynog_id']
              for (g, _), genes in order.items() for gene in genes}
    for key, genes in affected.items():
        order[key] = [{'genome': key[0], 'chrom': key[1],
                       'position': gene['position'], 'protein_id': gene['protein_id'],
                       'cynog_id': gene_labels.get((key[0], gene['protein_id']), '')}
                      for gene in genes]
    positions = defaultdict(list)
    by_protein = defaultdict(list)
    for key, genes in order.items():
        genes.sort(key=lambda r: int(r['position']))
        for i, gene in enumerate(genes):
            positions[key[0], gene['protein_id']].append((key, i))
            by_protein[gene['protein_id']].append((key, i))

    members = read(member_path)
    assigned = [r for r in members if r['assignment_status'] == 'assigned']
    if len(members) != 675 or len({r['accepted_locus_id'] for r in members}) != 675 or not assigned:
        raise SystemExit('Accepted membership contract failed')
    by_locus = defaultdict(list)
    accepted_targets = {r['builder_target_id'] for r in assigned}
    for row in assigned:
        by_locus[row['builder_locus_id']].append(row)
    result = []
    for locus_id, copies in by_locus.items():
        genome_keys = set()
        occurrences = []
        for copy in copies:
            hits = by_protein.get(copy['builder_target_id'], [])
            occurrences.extend(hits)
            genome_keys.update(key[0] for key, _ in hits)
        if not occurrences or len(genome_keys) != 1:
            raise SystemExit(f'Cannot bind ordered locus {locus_id}')
        # Placement coordinate_subject can name the RagTag reference scaffold,
        # while the native locus uses its ordered ENA source scaffold.
        candidates = [(key, i) for key, i in occurrences if key[1] in locus_id]
        if not candidates:
            raise SystemExit(f'No ordered source scaffold for {locus_id}')
        key = candidates[0][0]
        if key[0] not in labels:
            raise SystemExit(f'Missing paper species label for {key[0]}')
        scaffold = key[1]
        indices = [i for k, i in candidates if k == key]
        genes = order[key]
        left, right = min(indices), max(indices)
        cells = {}
        for side, start, step in (('U', left-1, -1), ('D', right+1, 1)):
            neighbours = [gene for index, gene in enumerate(genes)
                          if (index <= start if side == 'U' else index >= start)
                          and gene['protein_id'] not in accepted_targets]
            if side == 'U': neighbours.reverse()
            for distance in range(1, 13):
                gene = neighbours[distance-1] if distance <= len(neighbours) else None
                cells[f'{side}{distance}'] = f"{gene['protein_id']}|{gene['cynog_id']}" if gene else ''
        for copy in copies:
            result.append({**{k: copy[k] for k in ('family', 'block_id', 'accepted_locus_id')},
                           'species': labels[key[0]]} |
                          {'scaffold': scaffold} | cells)
    if len(result) != len(assigned) or len({r['accepted_locus_id'] for r in result}) != len(assigned):
        raise SystemExit('Flank output is not one row per assigned copy')
    old = {r['accepted_locus_id']: r for r in read(old_path)}
    current = {r['accepted_locus_id']: r for r in result}
    shared = set(old) & set(current)
    if not shared:
        raise SystemExit('No historical shared IDs for flank regression')
    mismatch = []
    changed_order = {'Acraspis_erinacei', 'Andricus_quercuspetiolicola',
                     'Druon_quercuslanigerum'}
    duplicated_old = {'Amphibolips_quercusostensackenii', 'Dryocosmus_kuriphilus'}
    gate_fields = [field for field in FIELDS if field not in ('block_id', 'species')]
    for copy_id in sorted(shared):
        duplicate_case = current[copy_id]['species'] in duplicated_old
        for field in gate_fields:
            if duplicate_case and field in FLANKS:
                continue
            if current[copy_id][field] != old[copy_id][field]:
                mismatch.append((copy_id, current[copy_id]['species'], field, old[copy_id][field], current[copy_id][field]))
        if duplicate_case:
            for side in ('U', 'D'):
                keys = [f'{side}{i}' for i in range(1,13)]
                values = [old[copy_id][key] for key in keys]
                collapsed = [v for i,v in enumerate(values) if v and (i==0 or v!=values[i-1])]
                candidate = [current[copy_id][key] for key in keys]
                for i,(was,now) in enumerate(zip(collapsed,candidate),1):
                    if was != now:
                        mismatch.append((copy_id,current[copy_id]['species'],f'{side}{i}_deduplicated',was,now))
    # The block is a corrected-run assignment, never a historical flank cell.
    for row in result:
        member = next((r for r in assigned if r['accepted_locus_id'] == row['accepted_locus_id']), None)
        if member is None or row['block_id'] != member['block_id']:
            raise SystemExit(f"Corrected block membership mismatch: {row['accepted_locus_id']}")
    report_path = Path(report_path)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    with report_path.open('w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['accepted_locus_id', 'species', 'field', 'historical', 'corrected'])
        writer.writerows(mismatch)
    unaffected = [r for r in mismatch if r[1] not in changed_order]
    if unaffected:
        raise SystemExit(f'{len(unaffected)} unchanged-genome flank cells mismatch; see {report_path}')
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(sorted(result, key=lambda r: r['accepted_locus_id']))

def main():
    parser = argparse.ArgumentParser()
    for name in ('order', 'affected-order', 'membership', 'historical', 'species-labels', 'out', 'report'):
        parser.add_argument('--' + name, required=True)
    parser.add_argument('--order-sha256', required=True)
    parser.add_argument('--membership-sha256', required=True)
    parser.add_argument('--affected-order-sha256', required=True)
    args = parser.parse_args()
    checked(args.order, args.order_sha256)
    checked(args.membership, args.membership_sha256)
    checked(args.affected_order, args.affected_order_sha256)
    produce(args.order, args.affected_order, args.membership, args.historical, args.species_labels, args.out, args.report)

if __name__ == '__main__':
    main()
