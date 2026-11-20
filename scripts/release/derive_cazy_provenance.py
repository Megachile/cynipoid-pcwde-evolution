#!/usr/bin/env python3
"""Make a portable record for the provenance-only Mayetiola mapping source."""
import argparse
import hashlib
from pathlib import Path

ORIGINAL_SHA = '115956861ea001a889e49a6555f04c8b336ff7150b3f7c380ff1ecf7c3bcc989'
PROTEOME_SHA = '183fe79cfd45f4c9899bbf2135e0059896dc5da1c5ce12c91652a9bf9828a933'


def derive(source):
    raw = source.read_bytes()
    if hashlib.sha256(raw).hexdigest() != ORIGINAL_SHA:
        raise ValueError('Mayetiola source identity changed')
    text = raw.decode('utf-8')
    required = ('Ma_destructor', '2026-05-24', 'i5K OGS1.0',
                '22,635 protein records', PROTEOME_SHA,
                'C 92.5, S 79.6, D 12.9, F 2.8, M 4.8')
    if any(item not in text for item in required):
        raise ValueError('Mayetiola provenance facts changed')
    return (f'# Mayetiola destructor — portable mapping provenance\n\n'
            f'Original record SHA256: `{ORIGINAL_SHA}`. Original record: '
            '`good_data/references/genomes/Ma_destructor/PROVENANCE.md` '
            '(retained privately).\n\n'
            'The canonical layout record was created on 2026-05-24 by the '
            '`genome_intake_for_cynogv4` goal. It records the i5K OGS1.0 '
            'peptide set (`hf_OGS1.0.pep.fasta`) as a proteome-only outgroup '
            'annotation: 22,635 protein records, SHA256 '
            f'`{PROTEOME_SHA}`. BUSCO 6.0.0 protein mode, `insecta_odb10`: '
            'C 92.5, S 79.6, D 12.9, F 2.8, M 4.8. No GFF, transcript '
            'FASTA or annotation TSV is claimed.\n\n'
            'Gene mapping method: `i5k_ogs_transcript_suffix`; remove the '
            'terminal `-R[A-Z][0-9]*` transcript suffix from each protein ID. '
            'The mapping values come from protein IDs; this record supplies '
            'provenance and a hash pin only.\n').encode('utf-8')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('source', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    content = derive(args.source)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    if not args.output.exists() or args.output.read_bytes() != content:
        args.output.write_bytes(content)
    print(f'PORTABLE_PROVENANCE_PASS original={ORIGINAL_SHA} derived={hashlib.sha256(content).hexdigest()} bytes={len(content)}')
