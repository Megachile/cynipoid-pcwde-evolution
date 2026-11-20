#!/usr/bin/env python3
"""Freeze the eleven catalog BUSCO values used by the Table S1 source producer.

Invoke with a host-local snapshot returned by pipeline_db_request.py snapshot.
The snapshot is opened read-only; the export records its exact byte hash.
"""
import csv
import hashlib
import sqlite3
import sys
from pathlib import Path

TARGETS = [(15, 157, "proteome", "42.7"), (17, 152, "proteome", "45.4"),
           (22, 150, "proteome", "73.9"), (24, 155, "proteome", "70.8"),
           (27, 205, "proteome", "76.2"), (49, 156, "proteome", "39.9"),
           (51, 149, "proteome", "91.1"), (54, 154, "proteome", "34.0"),
           (66, 151, "proteome", "52.5"), (68, 153, "proteome", "66.4"),
           (94, 114, "assembly", "95.7")]

def main():
    db, out = map(Path, sys.argv[1:3])
    if not str(db).startswith('/tmp/'):
        raise SystemExit('Snapshot must be host-local under /tmp')
    digest = hashlib.sha256(db.read_bytes()).hexdigest()
    con = sqlite3.connect(f'file:{db}?mode=ro', uri=True)
    con.row_factory = sqlite3.Row
    with out.open('w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['s1_row', 'accession', 'assembly_id', 'annotation_id', 'lineage',
                    'mode', 'score', 'catalog_record_id', 'snapshot_digest'])
        for row, ident, mode, expected in TARGETS:
            if mode == 'assembly':
                q = con.execute('SELECT assembly_id, accession, busco_lineage, busco_complete FROM assemblies WHERE assembly_id=?', (ident,)).fetchone()
                accession, assembly, annotation, lineage, score = q['accession'], ident, '', q['busco_lineage'], q['busco_complete']
            else:
                q = con.execute('SELECT a.annotation_id, a.assembly_id, a.proteome_busco_lineage, a.proteome_busco_C, b.accession FROM annotations a JOIN assemblies b ON a.assembly_id=b.assembly_id WHERE a.annotation_id=?', (ident,)).fetchone()
                accession, assembly, annotation, lineage, score = q['accession'], q['assembly_id'], ident, q['proteome_busco_lineage'], q['proteome_busco_C']
            if f'{score:.1f}' != expected:
                raise SystemExit(f'Catalog value changed for S1 row {row}: {score} != {expected}')
            w.writerow([row, accession, assembly, annotation, lineage, mode, f'{score:.1f}', ident, digest])
    print(f'{out}: {digest}')

if __name__ == '__main__':
    main()
