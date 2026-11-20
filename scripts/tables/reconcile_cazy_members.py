#!/usr/bin/env python3
"""Map accepted member sequences to processed scan proteins without ID heuristics."""

import argparse
import csv
import hashlib
from collections import defaultdict
from pathlib import Path

from cazy_gene_rule import map_genes
from cazy_data_root import resolve_input

FAMILIES = ('GH5_2', 'PL1_8', 'PL4_4')


def rows(path):
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def sequences(path):
    found = defaultdict(set)
    name, parts = None, []
    with path.open() as handle:
        for line in handle:
            if line.startswith('>'):
                if name is not None:
                    found[name].add(hashlib.sha256(''.join(parts).encode()).hexdigest())
                name, parts = line[1:].split()[0], []
            else:
                parts.append(line.strip().upper())
        if name is not None:
            found[name].add(hashlib.sha256(''.join(parts).encode()).hexdigest())
    return found


def sequence_text(path):
    found = defaultdict(set)
    name, parts = None, []
    with path.open() as handle:
        for line in handle:
            if line.startswith('>'):
                if name is not None:
                    found[name].add(''.join(parts))
                name, parts = line[1:].split()[0], []
            else:
                parts.append(line.strip().upper())
        if name is not None:
            found[name].add(''.join(parts))
    return found


def relation(accepted, processed):
    if not accepted:
        return 'unavailable_sequence', '', ''
    if accepted == processed:
        return 'exact_sequence', '1', str(len(processed))
    start = processed.find(accepted)
    if start >= 0:
        return 'contiguous_accepted_fragment', str(start + 1), str(start + len(accepted))
    # Only the documented scan-input replacement of invalid dot residues is allowed.
    if accepted.replace('.', 'X') == processed:
        return 'normalization_only', '1', str(len(processed))
    return 'other_sequence_difference', '', ''


def write(path, fields, values):
    with path.open('w', newline='') as handle:
        w = csv.writer(handle, delimiter='\t', lineterminator='\n')
        w.writerow(fields)
        w.writerows(values)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--manifest', type=Path, required=True)
    p.add_argument('--processed-map', type=Path, required=True)
    p.add_argument('--member-table', type=Path, required=True)
    p.add_argument('--sequences-dir', type=Path, required=True)
    p.add_argument('--data-root', '--source-root', dest='source_root', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    manifest = {r['species']: r for r in rows(a.manifest)}
    processed = {r['species']: r for r in rows(a.processed_map)}
    if len(manifest) != 64 or set(processed) != set(manifest):
        raise ValueError('expected 64 unique processed scan inputs')
    accepted = {fam: sequence_text(a.sequences_dir / f'complete_locus_{fam}.faa')
                for fam in FAMILIES}
    members = defaultdict(list)
    for r in rows(a.member_table):
        if r['family'] in FAMILIES:
            sp = r['species'].replace(' ', '_')
            if sp not in manifest:
                raise ValueError(f'accepted member outside Figure 1 roster: {sp}')
            members[sp].append(r)
    output = []
    for sp in sorted(manifest):
        item = manifest[sp]
        pin = processed[sp]
        path = resolve_input(a.source_root, pin['processed_path'])
        if sha(path) != pin['processed_sha256']:
            raise ValueError(f'{sp}: processed FASTA hash mismatch')
        scan = sequence_text(path)
        if any(len(v) != 1 for v in scan.values()):
            raise ValueError(f'{sp}: duplicate processed ID with conflicting sequence')
        reverse = defaultdict(set)
        for ident, seqs in scan.items():
            for seq in seqs:
                reverse[seq].add(ident)
        scan_file = resolve_input(a.source_root, item['scan_path'])
        if sha(scan_file) != item['scan_sha256']:
            raise ValueError(f'{sp}: scan output hash mismatch')
        hits = {r['Target Name'] for r in rows(scan_file)}
        for member in members[sp]:
            ident, fam = member['protein_id'], member['family']
            accepted_seqs = accepted[fam].get(ident, set())
            accepted_hashes = {hashlib.sha256(s.encode()).hexdigest() for s in accepted_seqs}
            if len(accepted_seqs) > 1:
                outcome, matches = 'ambiguous_accepted_sequence', set()
            elif ident in scan:
                matches = {ident}
                outcome = ('exact_id_same_sequence' if accepted_seqs == scan[ident]
                           else 'exact_id_sequence_unavailable' if not accepted_seqs
                           else 'exact_id_changed_sequence')
            elif len(accepted_seqs) == 1:
                matches = set().union(*(reverse[s] for s in accepted_seqs))
                outcome = ('unique_sequence_rename' if len(matches) == 1
                           else 'ambiguous_sequence_rename' if matches
                           else 'accepted_sequence_absent_from_scan')
            else:
                outcome, matches = 'accepted_sequence_unavailable', set()
            genes = map_genes(item['mapping_method'],
                              resolve_input(a.source_root, item['mapping_source']) if item['mapping_source'] else None,
                              matches) if matches else {}
            accepted_seq = next(iter(accepted_seqs)) if len(accepted_seqs) == 1 else ''
            processed_seq = next(iter(scan[next(iter(matches))])) if len(matches) == 1 else ''
            rel, start, end = relation(accepted_seq, processed_seq) if processed_seq else ('unavailable_sequence', '', '')
            if len(matches) > 1 or len(accepted_seqs) > 1:
                rel, start, end = 'ambiguous_sequence', '', ''
            if outcome == 'exact_id_changed_sequence':
                outcome = ('exact_id_contiguous_fragment' if rel == 'contiguous_accepted_fragment'
                           else 'exact_id_normalization_only' if rel == 'normalization_only'
                           else 'exact_id_other_sequence_difference')
            output.append((sp, fam, ident, outcome, ','.join(sorted(matches)),
                           ','.join(sorted(set(genes.values()))),
                           'yes' if matches & hits else 'no',
                           next(iter(accepted_hashes)) if len(accepted_hashes) == 1 else '',
                           rel, len(accepted_seq) if accepted_seq else '',
                           len(processed_seq) if processed_seq else '', start, end,
                           pin['processed_sha256'],
                           hashlib.sha256(processed_seq.encode()).hexdigest() if processed_seq else '',
                           'unique_id' if ident in scan and len(accepted_seqs) <= 1 else
                           'unique_sequence' if len(matches) == 1 else
                           'ambiguous' if matches or len(accepted_seqs) > 1 else 'unresolved'))
    a.out.mkdir(parents=True, exist_ok=False)
    fields = ('species', 'family', 'accepted_protein_id', 'outcome', 'scan_protein_ids',
              'scan_gene_ids', 'dbcan_sub_hit', 'accepted_sequence_sha256',
              'sequence_relation', 'accepted_length', 'processed_length',
              'fragment_start_1based', 'fragment_end_1based', 'processed_fasta_sha256',
              'processed_sequence_sha256', 'identity_status')
    write(a.out / 'member_crosswalk.tsv', fields, output)
    write(a.out / 'unresolved.tsv', fields,
          (r for r in output if r[3] not in ('exact_id_same_sequence', 'unique_sequence_rename',
                                            'exact_id_contiguous_fragment',
                                            'exact_id_normalization_only')))
    print(f'CAZY_MEMBER_CROSSWALK members={len(output)} '
          f'unresolved={sum(r[3] not in ("exact_id_same_sequence", "unique_sequence_rename", "exact_id_contiguous_fragment", "exact_id_normalization_only") for r in output)}')


if __name__ == '__main__':
    main()
