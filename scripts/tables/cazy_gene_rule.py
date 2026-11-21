"""Collapse dbCAN-sub protein hits to the gene units used for GH1_e72.

The mapping is selected by the pinned scan manifest, never by the hit label.
All CAZy groups therefore receive the same isoform treatment within a genome.
"""

import csv
import gzip
import re
from collections import defaultdict
from pathlib import Path
from urllib.parse import unquote


def _attrs(text):
    return {k: unquote(v) for part in text.split(';') if '=' in part
            for k, v in [part.split('=', 1)]}


def _gff_index(path, bk=False):
    """Protein key -> candidate gene keys, exactly as the GFF rule resolves them."""
    opener = gzip.open if path.suffix == '.gz' else open
    features = {}
    parents = defaultdict(set)
    direct = defaultdict(set)
    with opener(path, 'rt') as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 9:
                continue
            typ = fields[2]
            a = _attrs(fields[8])
            fid = a.get('ID', '')
            ps = set(filter(None, a.get('Parent', '').split(',')))
            if fid:
                features[fid] = (typ, ps)
            if typ != 'CDS':
                continue
            ids = set(filter(None, (a.get('protein_id', '') + ',' + a.get('proteinId', '')).split(',')))
            for x in a.get('Dbxref', '').split(','):
                if x.startswith('Genbank:') and not ids:
                    ids.add(x.split(':', 1)[1])
            for protein in ids:
                if bk:
                    if a.get('gene'):
                        direct[re.sub(r'\.\d+$', '', protein)].add(a['gene'])
                else:
                    parents[protein].update(ps)
                    if a.get('gene'):
                        direct[protein].add('gene_name:' + a['gene'])
                    for x in a.get('Dbxref', '').split(','):
                        if x.startswith('GeneID:'):
                            direct[protein].add(x)
    if bk:
        return dict(direct)

    def ancestors(starts):
        stack = list(starts)
        seen, genes = set(), set()
        while stack:
            x = stack.pop()
            if x in seen:
                continue
            seen.add(x)
            item = features.get(x)
            if item is None:
                continue
            typ, ps = item
            if typ == 'gene':
                genes.add(x)
            else:
                stack.extend(ps)
        return genes
    return {p: ancestors(ps) or direct[p] for p, ps in parents.items()}


def _gff_genes(path, wanted, bk=False):
    return _resolve(_gff_index(path, bk=bk), wanted)


def write_protein_gene_index(index, out):
    """Write a GFF-derived index; one row per (protein key, gene key), empty gene kept."""
    with Path(out).open('w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(('protein_key', 'gene_key'))
        for key in sorted(index):
            for gene in sorted(index[key]) or ['']:
                writer.writerow((key, gene))


def read_protein_gene_index(path):
    index = defaultdict(set)
    with Path(path).open(newline='') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            index[row['protein_key']]
            if row['gene_key']:
                index[row['protein_key']].add(row['gene_key'])
    return dict(index)


def _resolve(index, wanted):
    by_stem = defaultdict(set)
    for key, values in index.items():
        by_stem[re.sub(r'\.\d+$', '', key)].update(values)
    result = {}
    for protein in wanted:
        values = index.get(protein) or by_stem.get(re.sub(r'\.\d+$', '', protein), set())
        if len(values) != 1:
            raise ValueError(f'ambiguous or missing gene for {protein}: {sorted(values)}')
        result[protein] = next(iter(values))
    return result


def map_genes(method, mapping_source, proteins):
    """Return exactly one stable gene key for each hit protein."""
    proteins = set(proteins)
    if method in ('ncbi_gff', 'bk_ncbi_gff'):
        return _gff_genes(Path(mapping_source), proteins, bk=method == 'bk_ncbi_gff')
    if method == 'protein_gene_index':
        # Index derived from a public NCBI annotation GFF (Table S1 accession) by
        # scripts/deposit_inputs/reduce_public_references.py; same resolution rule.
        return _resolve(read_protein_gene_index(mapping_source), proteins)
    if method == 'canonical_annotation':
        with Path(mapping_source).open(newline='') as handle:
            index = {r['protein_id']: r.get('gene') or r.get('transcript_id') or
                     re.sub(r'\.\d+$', '', r['protein_id'])
                     for r in csv.DictReader(handle, delimiter='\t')}
        missing = proteins - set(index)
        if missing:
            raise ValueError(f'canonical annotation lacks {len(missing)} proteins: {sorted(missing)[:5]}')
        return {p: index[p] for p in proteins}
    if method == 'i5k_ogs_transcript_suffix':
        result = {p: re.sub(r'-R[A-Z][0-9]*$', '', p) for p in proteins}
        if any(result[p] == p for p in proteins):
            raise ValueError('i5k protein lacks expected transcript suffix')
        return result
    if method == 'frozen_model_stem':
        result = {p: re.sub(r'\.\d+$', '', re.sub(r'\.\d+\.p\d+$', '', p, flags=re.I)) for p in proteins}
        return {p: re.sub(r'(DKUR_\d+)-P\d+$', r'\1', v) for p, v in result.items()}
    if method == 'dbcan_duplicate_suffix':
        # dbCAN assigns _1, _2, ... to repeated FASTA headers. The suffix is
        # distinct from the paper model version after the final dot.
        return {p: re.sub(r'\.\d+$', '', re.sub(r'\.\d+_\d+$', '', p)) for p in proteins}
    raise ValueError(f'unsupported gene mapping method {method}')
