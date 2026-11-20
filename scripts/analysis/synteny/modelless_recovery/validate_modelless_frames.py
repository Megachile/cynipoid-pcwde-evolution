#!/usr/bin/env python3
"""Put every accepted model-less locus into the coordinate frame of the frozen
synteny gene order and validate it. Read-only, fai random access, no search.

Frames:
  native        scaffold already in the gene order for that genome; check bounds
                against the searched assembly.
  agp_liftover  recorded interval is on a pre-RagTag contig; lift through the
                RagTag AGP and require exact sequence identity native vs lifted.
  exact_match   pre-RagTag contig without an AGP: locate the exact interval
                sequence (either strand) in the accepted-chromosome RagTag
                scaffold; require one occurrence.
  no_order_scaffold  scaffold absent from gene order and no liftover: the locus
                cannot lie inside any detected search region (stays as is).

Every input (gene order, the three recovery tables, each .fai, AGP and FASTA
read) is checked against RECOVERY_PINS.tsv before derivation (recovery_pins.py);
a missing or changed pin stops the run before any output is written.
"""

import sys as _source_sys
from pathlib import Path as _SourcePath
_source_repo = next(p for p in _SourcePath(__file__).resolve().parents if (p / 'config/logical_inputs.tsv').is_file())
_source_sys.path.insert(0, str(_source_repo / 'scripts/release'))
from source_aliases import source_path
import csv, json, sys
from collections import defaultdict
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import recovery_pins

# usage: validate_modelless_frames.py PROJECT_ROOT RECOVERY_DIR
P = Path(sys.argv[1]).resolve()
O = Path(sys.argv[2]).resolve()
GO_KEY = 'source_inputs/synteny_corrected_placement/reconstruct/stage/final_gene_order/gene_order_final.tsv'
GO = source_path(P, GO_KEY)
BKCHR = {'chr3': 'NC_046659.1', 'chr6': 'NC_046662.1'}
GB = {'NC_046659.1': 'CM021340.1', 'NC_046662.1': 'CM021343.1'}
# genome -> (native contig fasta, RagTag AGP or None, RagTag fasta)
LIFT = {
    'GCA_035041075.1': ('good_data/references/genomes/Ac_erinacei/assembly.fasta',
                        'analysis/pcwde_effectors/pcwde_paper/analysis/ae_bk_ragtag_helixer_dbcan/ragtag_output_bk_asm10/ragtag.scaffold.agp',
                        'analysis/pcwde_effectors/pcwde_paper/analysis/ae_bk_ragtag_helixer_dbcan/ragtag_output_bk_asm10/ragtag.scaffold.fasta'),
    'GCA_010014905.2': ('genomes/druon/druon_genome.fna',
                        'analysis/pcwde_effectors/pcwde_paper/analysis/dl_bk_ragtag_helixer_validation_redo/ragtag_output_bk_asm10/ragtag.scaffold.agp',
                        'analysis/pcwde_effectors/pcwde_paper/analysis/dl_bk_ragtag_helixer_validation_redo/ragtag_output_bk_asm10/ragtag.scaffold.fasta'),
    'GCA_035222095.1': ('good_data/references/genomes/AQP/assembly.fasta', None,
                        'good_data/references/genomes/AQP/assembly.ragtag.secondary.fasta'),
}
COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
def rc(s): return s.translate(COMP)[::-1]
V = recovery_pins.Verifier()
def pkey(path):
    if Path(path) == GO:
        return GO_KEY
    return Path(path).relative_to(P).as_posix()
TABLES = ('COMPLETE_COORDINATE_RECOVERY.tsv', 'OVERLAP_TOLERANT_RECHAIN.tsv', 'RAW_LOCUS_RESCUE.tsv')
V.require([(pkey(GO), GO)] + [(recovery_pins.RECOVERY_PREFIX + n, O / n) for n in TABLES])

def fai(path):
    return {r[0]: (int(r[1]), int(r[2]), int(r[3]), int(r[4])) for r in csv.reader(open(V.check(Path(str(path) + '.fai'))), delimiter='\t')}
def fetch(path, name, start, end):
    ln, off, lb, lw = fai(path)[name]
    assert 1 <= start <= end <= ln, (name, start, end, ln)
    with open(V.check(path), 'rb') as f:
        a = off + (start - 1) // lb * lw + (start - 1) % lb
        f.seek(a)
        raw = f.read((end - start + 1) + ((end - start) // lb + 2) * (lw - lb))
    s = raw.replace(b'\n', b'').replace(b'\r', b'')[:end - start + 1].decode().upper()
    assert len(s) == end - start + 1
    return s

order = defaultdict(set)
for r in csv.DictReader(open(V.check(GO)), delimiter='\t'):
    order[r['genome']].add(r['chrom'])

rows = list(csv.DictReader(open(V.check(O / 'COMPLETE_COORDINATE_RECOVERY.tsv')), delimiter='\t'))
rechain = {r['accepted_id']: r for r in csv.DictReader(open(V.check(O / 'OVERLAP_TOLERANT_RECHAIN.tsv')), delimiter='\t')}
loci = []
for r in rows:
    x = dict(genome=r['genome'], family=r['family'], accepted_id=r['accepted_id'], accepted_bk_chromosome=r['accepted_bk_chromosome'],
             detection=('tblastn_only' if r['method'] else ''), method=r['method'], scaffold=r['scaffold'], start=r['start'], end=r['end'],
             strand=r['strand'], identity=r['identity'], query_coverage=r['query_coverage'], status=r['status'],
             assembly=r['assembly_project_path'], status_strict_chain=r['status'])
    if r['status'] != 'resolved' and r['accepted_id'] in rechain and rechain[r['accepted_id']]['status'] == 'resolved':
        q = rechain[r['accepted_id']]
        x.update(scaffold=q['scaffold'], start=q['start'], end=q['end'], strand=q['strand'], identity=q['identity'],
                 query_coverage=q['query_coverage'], status='resolved',
                 method='serial_subject_tblastn_collinear_global_bitscore_overlap_tolerant')
    loci.append(x)
for r in csv.DictReader(open(V.check(O / 'RAW_LOCUS_RESCUE.tsv')), delimiter='\t'):
    g = 'Isocolus_centaureae' if r['accession'] == 'Isocolus_centaureae' else r['accession']
    lo, hi = sorted((int(r['sstart']), int(r['send'])))
    loci.append(dict(genome=g, family=r['pcwde_family'], accepted_id='raw_rescue:' + r['best_contig'], accepted_bk_chromosome=r['bk_equivalent_chr'],
                     detection='tblastn_only_raw_rescue', method='recorded_tblastn_hsp', scaffold=r['best_contig'], start=lo, end=hi,
                     strand='+' if int(r['send']) >= int(r['sstart']) else '-', identity=r['best_percent_identity'], query_coverage='',
                     status='resolved', assembly='', status_strict_chain='resolved'))

# Pin every assembly index, AGP and FASTA the frame derivation below will read.
need = {}
for x in loci:
    if x['status'] != 'resolved':
        continue
    if x['scaffold'] in order.get(x['genome'], set()):
        if x['assembly']:
            f = source_path(P, x['assembly'] + '.fai'); need[pkey(f)] = f
    elif x['genome'] in LIFT:
        for v in LIFT[x['genome']]:
            if v:
                need[v] = source_path(P, v)
                if not v.endswith('.agp'):
                    need[v + '.fai'] = source_path(P, v + '.fai')
V.require(sorted(need.items()))
print('INPUT_PINS_VERIFIED', len(V.verified))

out = []
for x in loci:
    g, s = x['genome'], x['scaffold']
    y = dict(x, order_scaffold='', order_start='', order_end='', frame='', frame_check='', chromosome_concordance='')
    if x['status'] != 'resolved':
        y['frame'] = 'unresolved'
    elif s in order.get(g, set()):
        st, en = int(x['start']), int(x['end'])
        if x['assembly']:
            ln = fai(source_path(P, x['assembly']))[s][0]
            assert en <= ln
        y.update(order_scaffold=s, order_start=st, order_end=en, frame='native', frame_check='scaffold_in_gene_order;bounds_ok')
    elif g in LIFT:
        native, agp, rag = (source_path(P, v) if v else None for v in LIFT[g])
        st, en = int(x['start']), int(x['end'])
        seq = fetch(native, s, st, en)
        if agp:
            comp = [a for a in csv.reader(open(V.check(agp)), delimiter='\t') if len(a) > 8 and a[4] == 'W' and a[5] == s]
            assert len(comp) == 1, (s, len(comp))
            a = comp[0]; obj, os_, cs, ce, ori = a[0], int(a[1]), int(a[6]), int(a[7]), a[8]
            assert cs <= st and en <= ce
            if ori == '+': ls, le = os_ + st - cs, os_ + en - cs
            else: oe = int(a[2]); ls, le = oe - (en - cs), oe - (st - cs)
            got = fetch(rag, obj, ls, le)
            ok = got == (seq if ori == '+' else rc(seq))
            y.update(order_scaffold=obj, order_start=ls, order_end=le, frame='agp_liftover',
                     frame_check=f'agp_component_{ori};exact_sequence_{"match" if ok else "MISMATCH"}')
            assert ok
        else:
            obj = BKCHR[x['accepted_bk_chromosome']] + '_RagTag'
            ln = fai(rag)[obj][0]
            chrom = fetch(rag, obj, 1, ln)
            hits = []
            for strand, q in (('+', seq), ('-', rc(seq))):
                i = chrom.find(q)
                while i >= 0:
                    hits.append((strand, i + 1)); i = chrom.find(q, i + 1)
            assert len(hits) == 1, hits
            strand, i = hits[0]
            y.update(order_scaffold=obj, order_start=i, order_end=i + len(seq) - 1, frame='exact_match',
                     frame_check=f'unique_exact_{strand}_strand_{len(seq)}bp')
        assert y['order_scaffold'] in order[g]
    else:
        y.update(frame='no_order_scaffold', frame_check='scaffold_absent_from_gene_order_no_liftover')
    sc = y['order_scaffold']
    if sc:
        acc = x['accepted_bk_chromosome']
        names = {BKCHR.get(acc, '-'), GB.get(BKCHR.get(acc, ''), '-')}
        y['chromosome_concordance'] = ('concordant_bk_name' if sc.replace('_RagTag', '') in names else 'non_bk_name_scaffold')
    out.append(y)
w = csv.DictWriter(open(O / 'MODELLESS_LOCI_ORDER_FRAME.tsv', 'w'), fieldnames=list(out[0]), delimiter='\t')
w.writeheader(); w.writerows(out)
for y in out:
    print(y['accepted_id'], y['genome'], y['status'], y['frame'], y['order_scaffold'], y['order_start'], y['order_end'], y['frame_check'], y['chromosome_concordance'])

# Producer input: exactly the accepted copies (no raw rescue row, which is not
# one of the accepted copies without an order identity).
prod = [y for y in out if not y['accepted_id'].startswith('raw_rescue:')]
cols = ['genome', 'family', 'accepted_id', 'detection', 'method', 'status', 'frame', 'frame_check',
        'order_scaffold', 'order_start', 'order_end', 'identity', 'query_coverage', 'status_strict_chain', 'chromosome_concordance']
w = csv.DictWriter(open(O / 'MODELLESS_LOCI_DISPOSITION.tsv', 'w'), fieldnames=cols, delimiter='\t', extrasaction='ignore', lineterminator='\n')
w.writeheader(); w.writerows(prod)
print('DISPOSITION_ROWS', len(prod))
