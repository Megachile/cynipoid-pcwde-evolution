#!/usr/bin/env python3
"""Bounded contract check for scan gene IDs through annotation objects."""
import csv
import hashlib
import subprocess
import sys
import tempfile
from pathlib import Path


def write_tsv(path, header, rows):
    with path.open('w', newline='') as stream:
        writer = csv.writer(stream, delimiter='\t', lineterminator='\n')
        writer.writerow(header)
        writer.writerows(rows)


with tempfile.TemporaryDirectory() as temp:
    root = Path(temp)
    species = ['Acraspis_erinacei', 'AHIE', 'AMIN', 'CQUE', 'DRUB',
               'PCAN', 'PLAM', 'SHIS', 'TPHE']
    sources, agps, genes = [], [], []
    for name in species:
        gff = root / (name + '.gff3')
        agp = root / (name + '.agp')
        lines = []
        agp_lines = []
        for category, object_name in [('native', 'chr1_RagTag'),
                                      ('focal', 'chr2_RagTag')]:
            gene = name + '_' + category
            genes.append((name, gene + '.1', gene))
            source_object = object_name if name == 'Acraspis_erinacei' else gene + '_source'
            lines.append(f'{source_object}\tHelixer\tmRNA\t1\t9\t.\t+\t.\tID={gene}.1;Parent={gene}\n')
            agp_lines.append(f'{object_name}\t1\t9\t1\tW\t{source_object}\t1\t9\t+\n')
        # Multiple source objects must remain visible, even when both place.
        gene = name + '_ambiguous'
        genes.append((name, gene + '.1', gene))
        for obj in ['source_a', 'source_b']:
            lines.append(f'{obj}\tHelixer\tmRNA\t1\t9\t.\t+\t.\tID={gene}.1;Parent={gene}\n')
        gff.write_text(''.join(lines))
        agp.write_text(''.join(agp_lines))
        sources.append((name, gff.name, gff.stat().st_size,
                        hashlib.sha256(gff.read_bytes()).hexdigest()))
        if name != 'Acraspis_erinacei':
            agps.append((name, agp.name, agp.stat().st_size,
                         hashlib.sha256(agp.read_bytes()).hexdigest()))
    source_path, agp_path, gene_path, output = (root / name for name in
                                      ['sources.tsv', 'agps.tsv', 'genes.tsv', 'bridge.tsv'])
    write_tsv(source_path, ['species', 'project_path', 'bytes', 'sha256'], sources)
    write_tsv(agp_path, ['species', 'project_path', 'bytes', 'sha256'], agps)
    write_tsv(gene_path, ['species', 'protein_id', 'gene_id'], genes)
    subprocess.run([sys.executable, str(Path(__file__).with_name('ecami_bridge.py')),
                    '--sources', str(source_path), '--agp-sources', str(agp_path),
                    '--project-root', str(root),
                    '--protein-genes', str(gene_path), '--out', str(output)], check=True)
    with output.open() as stream:
        rows = list(csv.DictReader(stream, delimiter='\t'))
    assert len(rows) == len(species) * 3
    for row in rows:
        expected = 'multiple' if row['gene_id'].endswith('_ambiguous') else 'unique'
        assert row['object_status'] == expected, row
        assert len(row['ragtag_objects'].split(',')) == (2 if expected == 'multiple' else 1)
        if expected == 'unique':
            assert row['ragtag_objects'] in {'chr1_RagTag','chr2_RagTag'}
    print(f'BRIDGE_FIXTURE_PASS species={len(species)} genes={len(rows)} native_and_focal=yes ambiguous_unplaced=yes')
