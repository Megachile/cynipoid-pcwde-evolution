#!/usr/bin/env python3
"""Normalize invalid protein symbols for a pinned competitive dbCAN-sub scan.

The source FASTA is preserved. Headers, order, and sequence lengths are unchanged;
the producer's existing rule replaces non-letter amino-acid placeholders with X.
dbCAN itself assigns suffixes to duplicate FASTA identifiers on import.
"""

import argparse
import csv
import hashlib
from collections import Counter
from pathlib import Path

ALLOWED = set('ABCDEFGHIKLMNPQRSTUVWXYZ')


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--source', type=Path, required=True)
    p.add_argument('--source-sha256', required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--report', type=Path, required=True)
    args = p.parse_args()
    if sha(args.source) != args.source_sha256:
        raise ValueError('paper source FASTA hash mismatch')
    if args.out.exists() or args.report.exists():
        raise FileExistsError('scan input or report already exists')
    ids, replacement = Counter(), Counter()
    residues = records = 0
    with args.source.open() as source, args.out.open('w') as dest:
        for line in source:
            if line.startswith('>'):
                records += 1
                ids[line[1:].split()[0]] += 1
                dest.write(line)
            else:
                clean = []
                for char in line.rstrip('\n'):
                    if char.upper() in ALLOWED:
                        clean.append(char.upper())
                    else:
                        clean.append('X')
                        replacement[char] += 1
                    residues += 1
                dest.write(''.join(clean) + '\n')
    with args.report.open('w', newline='') as handle:
        w = csv.writer(handle, delimiter='\t', lineterminator='\n')
        w.writerow(('metric', 'value'))
        w.writerows((('source_sha256', args.source_sha256),
                     ('scan_input_sha256', sha(args.out)),
                     ('records', records), ('unique_ids', len(ids)),
                     ('duplicate_ids', sum(v > 1 for v in ids.values())),
                     ('residues', residues),
                     ('total_replacements', sum(replacement.values()))))
        for char, n in sorted(replacement.items()):
            w.writerow(('replacement_' + repr(char), n))
    print(f'CAZY_INPUT_PREPARED records={records} unique_ids={len(ids)} '
          f'replacements={sum(replacement.values())} sha256={sha(args.out)}')


if __name__ == '__main__':
    main()
