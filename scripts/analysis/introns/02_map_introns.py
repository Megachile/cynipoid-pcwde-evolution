#!/usr/bin/env python3
"""Map current-roster introns to accepted BK reference protein coordinates."""
import argparse,csv,json
from pathlib import Path
from Bio import SeqIO
from Bio.Align import PairwiseAligner,substitution_matrices

FAMILIES=('GH5_2','PL1_8','PL4_4')
REFS={'GH5_2':'XP_033220919','PL1_8':'XP_033225824','PL4_4':'XP_033222742'}

def norm(x): return x[:-2] if x.endswith('.1') else x
def chrom(x):
    x=str(x or '').removeprefix('chr')
    return 'unplaced' if x in ('','NA','nan','None','unknown','unresolved') else x
def load_map(legacy_map,synteny_map):
    with legacy_map.open() as h: old={norm(x['protein_id']):chrom(x['bk_chr']) for x in csv.DictReader(h,delimiter='\t')}
    with synteny_map.open() as h: new={norm(x['protein_id']):chrom(x['bk_chromosome']) for x in csv.DictReader(h,delimiter='\t')}
    conflicts={p:(old[p],new[p]) for p in old.keys()&new.keys() if old[p]!=new[p]}
    if conflicts: raise AssertionError(f'chromosome-map conflicts: {list(conflicts.items())[:5]}')
    return old,new
def map_positions(query,ref):
    a=PairwiseAligner();a.substitution_matrix=substitution_matrices.load('BLOSUM62');a.open_gap_score=-10;a.extend_gap_score=-0.5;a.mode='global'
    aln=a.align(ref,query)[0]; out={}
    for (rs,re),(qs,qe) in zip(aln.aligned[0],aln.aligned[1]):
        for ri,qi in zip(range(rs,re),range(qs,qe)):out[qi+1]=ri+1
    return out
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--input-dir',type=Path,required=True)
    ap.add_argument('--sequence-dir',type=Path,required=True,help='Built proteins selected from pinned annotation proteomes.')
    ap.add_argument('--member-table',type=Path,required=True,help='Accepted member-layer protein and BK chromosome classes.')
    ap.add_argument('--output-dir',type=Path,required=True)
    a=ap.parse_args()
    with a.member_table.open() as h: members=list(csv.DictReader(h,delimiter='\t'))
    member_by_id={r['protein_id']:r for r in members if r['family'] in FAMILIES}
    new={norm(r['protein_id']):chrom(r['bk_chromosome']) for r in members if r['family'] in FAMILIES}
    if len(new)!=629:raise AssertionError(f'expected 629 three-family member chromosomes, got {len(new)}')
    old={};a.output_dir.mkdir(parents=True,exist_ok=False);qc={}
    for fam in FAMILIES:
        seqs={norm(r.id):str(r.seq) for r in SeqIO.parse(a.sequence_dir/f'{fam}_current_genome_proteins.faa','fasta')}
        ref=seqs[REFS[fam]]
        with (a.input_dir/f'{fam}_intron_metadata.tsv').open() as h: raw=list(csv.DictReader(h,delimiter='\t'))
        bypid={}
        for r in raw:bypid.setdefault(r['protein_id'],[]).append(r)
        rows=[]; missing_seq=[]
        for pid,pr in sorted(bypid.items()):
            key=norm(pid); locus=new.get(key,old.get(key,'unplaced')); locus='unplaced' if locus=='unplaced' else 'chr'+locus
            if int(pr[0]['intron_count'])==0:
                rows.append(dict(protein_id=pid,species=member_by_id[pid]['species'],intron_number=0,aa_position=0,ref_aa_pos=0,phase_before=-1,phase_after=-1,donor_dinuc='',acceptor_dinuc='',locus_class=locus));continue
            if key not in seqs:missing_seq.append(pid);continue
            mp=map_positions(seqs[key],ref)
            for r in pr:
                pos=int(r['aa_position']); rp=mp.get(pos)
                if rp is None:
                    for d in (1,-1,2,-2,3,-3):
                        if pos+d in mp:rp=mp[pos+d];break
                rows.append(dict(protein_id=pid,species=member_by_id[pid]['species'],intron_number=r['intron_number'],aa_position=pos,ref_aa_pos=rp or -1,phase_before=r['phase_before'],phase_after=r['phase_after'],donor_dinuc=r['donor_dinuc'],acceptor_dinuc=r['acceptor_dinuc'],locus_class=locus))
        if missing_seq:raise AssertionError(f'{fam} sequences missing: {missing_seq}')
        cols=list(rows[0])
        with (a.output_dir/f'{fam}_intron_positions.tsv').open('w',newline='') as h:
            w=csv.DictWriter(h,fieldnames=cols,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
        qc[fam]={'eligible_proteins':len(bypid),'with_introns':sum(int(x[0]['intron_count'])>0 for x in bypid.values()),'mapped_intron_records':sum(int(r['ref_aa_pos'])>0 for r in rows),'unmapped_intron_records':sum(int(r['intron_number'])>0 and int(r['ref_aa_pos'])<1 for r in rows),'locus_counts':{}}
        for pid,pr in bypid.items():
            key=norm(pid);v=new.get(key,old.get(key,'unplaced'));v='unplaced' if v=='unplaced' else 'chr'+v
            qc[fam]['locus_counts'][v]=qc[fam]['locus_counts'].get(v,0)+1
    (a.output_dir/'QC.json').write_text(json.dumps(qc,indent=2,sort_keys=True)+'\n');print(json.dumps(qc,indent=2,sort_keys=True))
if __name__=='__main__':main()
