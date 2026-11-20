#!/usr/bin/env python3
"""Resolve corrected dbCAN-sub gene IDs to source RagTag objects.

The scan producer already collapses proteins to genes.  A gene can span more
than one source scaffold; retain all objects so placement can fail closed if
they disagree.  This producer never assigns BK chromosomes.
"""
import argparse
import csv
import hashlib
from collections import defaultdict
from pathlib import Path


def rows(path):
    with path.open(newline='') as handle:
        yield from csv.DictReader(handle, delimiter='\t')


def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def attributes(text):
    return {key: value for part in text.split(';') if '=' in part
            for key, value in [part.split('=', 1)]}


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--sources', type=Path, required=True)
    p.add_argument('--agp-sources', type=Path)
    p.add_argument('--project-root', type=Path, required=True)
    p.add_argument('--protein-genes', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    wanted = defaultdict(set)
    for row in rows(a.protein_genes):
        wanted[row['species']].add(row['gene_id'])
    found = defaultdict(lambda: defaultdict(set))
    sources = list(rows(a.sources))
    if len(sources) != len({r['species'] for r in sources}):
        raise SystemExit('duplicate bridge species')
    for source in sources:
        species = source['species']
        relative = Path(source['project_path'])
        if relative.is_absolute() or '..' in relative.parts:
            raise SystemExit('unsafe source path')
        path = a.project_root / relative
        if path.stat().st_size != int(source['bytes']) or sha256(path) != source['sha256']:
            raise SystemExit(f'corrected GFF source drift: {species}')
        # Restrict lookups to the scanned genes; a feature ID may have a
        # transcript, split-protein or exon suffix after the stable gene ID.
        by_prefix = wanted[species]
        with path.open() as handle:
            for line in handle:
                if line.startswith('#'):
                    continue
                fields = line.rstrip('\n').split('\t')
                if len(fields) != 9 or fields[2] not in {'gene', 'mRNA', 'CDS', 'exon'}:
                    continue
                attrs = attributes(fields[8])
                identifiers = [attrs.get('ID', ''), *attrs.get('Parent', '').split(',')]
                for identifier in identifiers:
                    candidate = identifier
                    while candidate not in by_prefix and '.' in candidate:
                        candidate = candidate.rsplit('.', 1)[0]
                    if candidate in by_prefix and (identifier == candidate or identifier.startswith(candidate + '.')):
                        found[species][candidate].add(fields[0])
    agp_objects = {}
    for source in rows(a.agp_sources) if a.agp_sources else []:
        species = source['species']
        relative = Path(source['project_path'])
        if relative.is_absolute() or '..' in relative.parts:
            raise SystemExit('unsafe AGP source path')
        path = a.project_root / relative
        if path.stat().st_size != int(source['bytes']) or sha256(path) != source['sha256']:
            raise SystemExit(f'AGP source drift: {species}')
        mapping = defaultdict(set)
        with path.open() as handle:
            for line in handle:
                if line.startswith('#'):
                    continue
                fields = line.rstrip('\n').split('\t')
                if len(fields) > 5 and fields[4] == 'W':
                    mapping[fields[5]].add(fields[0])
        if species in agp_objects:
            raise SystemExit(f'duplicate AGP source: {species}')
        agp_objects[species] = mapping
    if agp_objects and set(agp_objects) - {s['species'] for s in sources}:
        raise SystemExit('AGP species absent from annotation sources')
    output = []
    missing_annotation = 0
    for source in sources:
        species = source['species']
        for gene in sorted(wanted[species]):
            objects = sorted(found[species][gene])
            if not objects:
                missing_annotation += 1
            if species in agp_objects and len(objects) == 1:
                objects = sorted(agp_objects[species].get(objects[0], set()))
                status = 'unique' if len(objects) == 1 else ('multiple' if objects else 'missing_agp')
            else:
                status = 'unique' if len(objects) == 1 else ('multiple' if objects else 'missing')
            output.append((species, gene, ','.join(objects), status))
    a.out.parent.mkdir(parents=True, exist_ok=True)
    with a.out.open('w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(('species', 'gene_id', 'ragtag_objects', 'object_status'))
        writer.writerows(output)
    missing = sum(row[3] == 'missing' for row in output)
    print(f'ECAMI_BRIDGE rows={len(output)} missing={missing} multiple={sum(row[3] == "multiple" for row in output)}')
    if missing_annotation:
        raise SystemExit('scan genes missing from pinned annotation GFFs')


if __name__ == '__main__':
    main()
