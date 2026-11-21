#!/usr/bin/env python3
"""Replace redistributed public reference files in the input packet.

NCBI annotation GFFs used only for protein->gene collapse become compact
protein->gene index tables (same resolution rule, scripts/tables/cazy_gene_rule.py),
and outgroup reference proteomes used only as hash pins are dropped. Both are
public downloads identified by Table S1 accessions; their original SHA-256 values
stay pinned in config/public_reference_reductions.tsv so a downloaded copy can be
checked and the index regenerated with:

    python3 scripts/deposit_inputs/reduce_public_references.py --derive-index GFF OUT [--bk]

Usage (one-time packet transform):
    reduce_public_references.py --repo REPO --inputs DATA_ROOT/repository_inputs
"""
import argparse
import csv
import hashlib
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'tables'))
from cazy_gene_rule import _gff_index, write_protein_gene_index  # noqa: E402

# Source accessions for the public references, as listed in Table S1.
ACCESSIONS = {
    'Apis_mellifera': 'NCBI GCF_003254395.2', 'Athalia_rosae': 'NCBI GCF_917208135.1',
    'Belonocnema_kinseyi': 'NCBI GCF_010883055.1', 'Bombus_terrestris': 'NCBI GCF_910591885.1',
    'Bradysia_coprophila': 'NCBI GCF_014529535.1', 'Cephus_cinctus': 'NCBI GCF_000341935.2',
    'Ceratina_calcarata': 'NCBI GCF_001652005.2', 'Contarinia_nasturtii': 'NCBI GCF_009176525.2',
    'Drosophila_melanogaster': 'NCBI GCA_000001215.4', 'Habropoda_laboriosa': 'NCBI GCF_001263275.1',
    'Mayetiola_destructor': 'i5K GCA_000149185.1 (OGS1.0)', 'Microplitis_demolitor': 'NCBI GCF_026212275.2',
    'Nasonia_vitripennis': 'NCBI GCF_009193385.2', 'Neodiprion_lecontei': 'NCBI GCF_021901455.1',
    'Orussus_abietinus': 'NCBI GCF_000612105.2', 'Solenopsis_invicta': 'NCBI GCF_016802725.1',
}
# Kept: BK's NCBI GFF is also an upstream input for intron restaging and synteny member
# identity (config/intron_sources.tsv, intron_snapshot_sources.tsv, intron_manifest_bridge.json).
# BK's RefSeq proteome is also kept: accepted BK members are matched to it by sequence.
KEEP_GFF = {'Belonocnema_kinseyi'}
CONFIGS = ('config/logical_inputs.tsv', 'config/source_aliases.tsv',
           'config/cazy_stage_sources.tsv', 'config/upstream_locked_inputs.tsv')


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def read(path):
    with Path(path).open(newline='') as f:
        r = csv.reader(f, delimiter='\t')
        return next(r), list(r)


def write(path, header, rows):
    with Path(path).open('w', newline='') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n')
        w.writerow(header)
        w.writerows(rows)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--repo', type=Path)
    p.add_argument('--inputs', type=Path)
    p.add_argument('--derive-index', nargs=2, type=Path, metavar=('GFF', 'OUT'))
    p.add_argument('--bk', action='store_true')
    a = p.parse_args()
    if a.derive_index:
        write_protein_gene_index(_gff_index(a.derive_index[0], bk=a.bk), a.derive_index[1])
        return
    repo, inputs = a.repo, a.inputs
    head, scans = read(repo / 'config/cazy_scans.tsv')
    col = {k: i for i, k in enumerate(head)}
    removed = {}      # repository_inputs-relative path -> (sha, bytes)
    replaced = {}     # old rel -> new rel
    record = []
    for row in scans:
        sp, method, source = row[col['species']], row[col['mapping_method']], row[col['mapping_source']]
        if method in ('ncbi_gff', 'bk_ncbi_gff') and sp not in KEEP_GFF:
            assert source.startswith('repository_inputs/cazy_inputs/mapping/'), source
            old = inputs / source.split('/', 1)[1]
            old_sha, old_bytes = sha(old), old.stat().st_size
            assert old_sha == row[col['mapping_sha256']], sp
            name = old.name.split('_mapping.gff')[0] + '_protein_gene.tsv'
            new = old.with_name(name)
            write_protein_gene_index(_gff_index(old, bk=method == 'bk_ncbi_gff'), new)
            new_rel = 'repository_inputs/cazy_inputs/mapping/' + name
            removed[source] = (old_sha, old_bytes)
            replaced[source] = (new_rel, sha(new), new.stat().st_size)
            record.append((sp, ACCESSIONS[sp], 'annotation_gff', source, old_sha, old_bytes,
                           method, new_rel, replaced[source][1]))
            row[col['mapping_method']] = 'protein_gene_index'
            row[col['mapping_source']] = new_rel
            row[col['mapping_sha256']] = replaced[source][1]
        if row[col['source']] == 'outgroup_15':
            prot = row[col['proteome_path']]
            f = inputs / prot.split('/', 1)[1]
            assert sha(f) == row[col['proteome_sha256']], sp
            removed[prot] = (row[col['proteome_sha256']], f.stat().st_size)
            record.append((sp, ACCESSIONS[sp], 'reference_proteome', prot, removed[prot][0],
                           removed[prot][1], 'hash_pin_only', '', ''))
            row[col['proteome_path']] = ''
    write(repo / 'config/cazy_scans.tsv', head, scans)
    for rel, field in (('config/cazy_processed_proteomes.tsv', 'processed_path'),
                       ('config/cazy_paper_proteomes.tsv', 'proteome_path')):
        h, rows = read(repo / rel)
        i = h.index(field)
        for r in rows:
            if r[i] in removed:
                r[i] = ''
        write(repo / rel, h, rows)

    def tail(path):
        return path.split('/', 1)[1] if path.startswith(('repository_inputs/', 'source_inputs/')) else path

    names = {Path(k).name: k for k in removed}

    def hit(r):
        return next((names[Path(x).name] for x in r if Path(x).name in names), None)

    def clone(r, old_rel):
        new_rel, new_sha, new_bytes = replaced[old_rel]
        old_sha, old_bytes = removed[old_rel]
        oname, nname = Path(old_rel).name, Path(new_rel).name
        out = []
        for x in r:
            x = x.replace(oname, nname).replace(old_sha, new_sha)
            if x == str(old_bytes):
                x = str(new_bytes)
            x = x.replace(new_sha + '.gff.gz', new_sha + '.tsv').replace(new_sha + '.gff', new_sha + '.tsv')
            x = x.replace('annotation GFF', 'protein-to-gene index from annotation GFF')
            x = x.replace('gene mapping source', 'gene mapping index derived from NCBI GFF ' + old_sha)
            x = x.replace('(computed_in_place_2026-11-02)', '(derived by scripts/deposit_inputs/reduce_public_references.py)')
            x = x.replace('computed in place for r41; exact input', 'derived from the pinned public GFF; exact input')
            out.append(x)
        return out
    manifests = [repo / c for c in CONFIGS] + [inputs / 'INPUT_MANIFEST.tsv']
    for m in manifests:
        h, rows = read(m)
        out, n_drop, n_add = [], 0, 0
        for r in rows:
            old = hit(r)
            if old is None:
                out.append(r)
                continue
            n_drop += 1
            if old in replaced:
                out.append(clone(r, old))
                n_add += 1
        write(m, h, out)
        print(f'{m.name}: dropped {n_drop}, replaced {n_add}')
    for rel in removed:
        (inputs / rel.split('/', 1)[1]).unlink()
    write(repo / 'config/public_reference_reductions.tsv',
          ('species', 'table_s1_accession', 'removed_kind', 'removed_path', 'removed_sha256',
           'removed_bytes', 'original_method', 'derived_path', 'derived_sha256'), sorted(record))
    print(f'REDUCED removed={len(removed)} bytes={sum(b for _, b in removed.values())} derived={len(replaced)}')


if __name__ == '__main__':
    main()
