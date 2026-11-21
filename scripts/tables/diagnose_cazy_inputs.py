#!/usr/bin/env python3
"""Check scan input identity and accepted-member sequence reconciliation.

This is a diagnostic gate. It reports unresolved evidence; it never repairs IDs
or promotes the inherited scan inventory to paper-proteome authority.
"""

import argparse
import csv
import hashlib
from collections import Counter, defaultdict
from pathlib import Path
from cazy_data_root import resolve_input


def table(path):
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def write(path, columns, rows):
    with path.open('w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(columns)
        writer.writerows(rows)


def file_hash(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def fasta(path):
    result = defaultdict(set)
    name, pieces = None, []
    with path.open() as handle:
        for line in handle:
            if line.startswith('>'):
                if name is not None:
                    result[name].add(hashlib.sha256(''.join(pieces).encode()).hexdigest())
                name, pieces = line[1:].split()[0], []
            else:
                pieces.append(line.strip().upper())
        if name is not None:
            result[name].add(hashlib.sha256(''.join(pieces).encode()).hexdigest())
    return result


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--manifest', type=Path, required=True)
    p.add_argument('--paper-sources', type=Path, required=True)
    p.add_argument('--member-table', type=Path, required=True)
    p.add_argument('--data-root', '--source-root', dest='source_root', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    args = p.parse_args()
    manifest = table(args.manifest)
    if len(manifest) != 64 or len({r['species'] for r in manifest}) != 64:
        raise ValueError('expected 64 distinct Figure 1 rows')
    paper = {r['species']: r for r in table(args.paper_sources)}
    if len(paper) != len(table(args.paper_sources)) or not set(paper) <= {r['species'] for r in manifest}:
        raise ValueError('invalid paper-source roster')
    members = defaultdict(list)
    for row in table(args.member_table):
        if row['family'] in {'GH5_2', 'PL1_8', 'PL4_4'}:
            species = row['species'].replace(' ', '_')
            if species not in {r['species'] for r in manifest}:
                raise ValueError(f'accepted member outside roster: {species}')
            members[species].append(row)
    coverage, reconciliation, unresolved = [], [], []
    for row in manifest:
        sp = row['species']
        scan = resolve_input(args.source_root, row['scan_path'])
        if not row['proteome_path']:
            # Public reference proteome (Table S1 accession), kept only as a SHA-256 pin.
            if members[sp] or file_hash(scan) != row['scan_sha256']:
                raise ValueError(f'{sp}: unredistributed proteome with members or scan hash mismatch')
            coverage.append((sp, row['source_key'], '', row['proteome_sha256'], row['scan_path'],
                             row['scan_sha256'], '', '', 'public_reference_not_redistributed'))
            continue
        scan_input = resolve_input(args.source_root, row['proteome_path'])
        if file_hash(scan_input) != row['proteome_sha256'] or file_hash(scan) != row['scan_sha256']:
            raise ValueError(f'{sp}: scan manifest hash mismatch')
        candidate = paper.get(sp)
        paper_hash = ''
        relation = 'paper_source_unverified'
        paper_path = ''
        if candidate:
            paper_path = candidate['proteome_path']
            source = resolve_input(args.source_root, paper_path)
            paper_hash = file_hash(source)
            if paper_hash != candidate['proteome_sha256']:
                raise ValueError(f'{sp}: paper source hash mismatch')
            relation = 'identical' if paper_hash == row['proteome_sha256'] else 'different'
        coverage.append((sp, row['source_key'], row['proteome_path'], row['proteome_sha256'],
                         row['scan_path'], row['scan_sha256'], paper_path, paper_hash, relation))
        if not members[sp]:
            continue
        scan_seqs = fasta(scan_input)
        paper_seqs = fasta(resolve_input(args.source_root, paper_path)) if paper_path else {}
        scan_by_seq = defaultdict(set)
        for ident, seqs in scan_seqs.items():
            for seq in seqs:
                scan_by_seq[seq].add(ident)
        for member in members[sp]:
            ident = member['protein_id']
            seqs = paper_seqs.get(ident, set())
            if ident in scan_seqs:
                if seqs and scan_seqs[ident] != seqs:
                    outcome, matches = 'same_id_changed_sequence', ident
                else:
                    outcome, matches = 'exact_scan_id', ident
            elif len(seqs) == 1:
                matches_set = set().union(*(scan_by_seq[seq] for seq in seqs))
                matches = ','.join(sorted(matches_set))
                outcome = ('sequence_rename_one_to_one' if len(matches_set) == 1
                           else 'sequence_ambiguous' if matches_set else 'changed_model_or_missing_scan')
            elif len(seqs) > 1:
                outcome, matches = 'duplicate_paper_id_conflicting_sequences', ''
            else:
                outcome, matches = ('paper_id_absent' if paper_path else 'paper_source_unverified'), ''
            item = (sp, member['family'], ident, paper_path, outcome, matches)
            reconciliation.append(item)
            if outcome != 'exact_scan_id' and outcome != 'sequence_rename_one_to_one':
                unresolved.append(item)
    args.out.mkdir(parents=True, exist_ok=False)
    write(args.out / 'coverage.tsv', ('species', 'source_key', 'scan_proteome_path',
          'scan_proteome_sha256', 'scan_path', 'scan_sha256', 'paper_proteome_path',
          'paper_proteome_sha256', 'paper_scan_relation'), coverage)
    columns = ('species', 'family', 'accepted_protein_id', 'paper_proteome_path', 'outcome', 'scan_matches')
    write(args.out / 'member_reconciliation.tsv', columns, reconciliation)
    write(args.out / 'unresolved.tsv', columns, unresolved)
    counts = Counter(r[4] for r in reconciliation)
    print(f'CAZY_DIAGNOSTIC rows={len(coverage)} members={len(reconciliation)} '
          f'unresolved={len(unresolved)} outcomes={dict(sorted(counts.items()))}')


if __name__ == '__main__':
    main()
