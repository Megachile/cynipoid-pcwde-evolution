#!/usr/bin/env python3
import argparse,json
from collections import Counter
from pathlib import Path
import pandas as pd
FAMILIES=('GH5_2','PL1_8','PL4_4'); MAX_DIST=3
FLANK=12
def readfa(p):
    d={};k=None;a=[]
    for line in open(p):
        if line.startswith('>'):
            if k is not None:d[k]=''.join(a)
            k=line[1:].split()[0];a=[]
        else:a.append(line.strip())
    if k is not None:d[k]=''.join(a)
    return d
def window(seq,pos):
    i=pos-1
    return seq[max(0,i-FLANK):i].rjust(FLANK,'.')+seq[i:i+FLANK].ljust(FLANK,'.')
def mode_seek(items):
    alive=list(range(len(items)));out=[]
    while alive:
        best=max(alive,key=lambda i:sum(abs(items[j][0]-items[i][0])<=MAX_DIST for j in alive))
        members=[j for j in alive if abs(items[j][0]-items[best][0])<=MAX_DIST]
        ps=sorted(items[j][0] for j in members);out.append((ps[len(ps)//2],members));ms=set(members);alive=[j for j in alive if j not in ms]
    return sorted(out)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--input-dir',type=Path,required=True);ap.add_argument('--sequence-dir',type=Path,required=True);ap.add_argument('--output-dir',type=Path,required=True);a=ap.parse_args();a.output_dir.mkdir(parents=True,exist_ok=False);qc={};summary=[]
    for fam in FAMILIES:
        d=pd.read_csv(a.input_dir/f'{fam}_intron_positions.tsv',sep='\t').reset_index(drop=True);d['group_id']=''
        idx=d.index[(d.ref_aa_pos>0)&(d.intron_number>0)].tolist();items=[(int(d.loc[i,'ref_aa_pos']),i) for i in idx];groups=[]
        for center,members in mode_seek(items):
            src=[items[j][1] for j in members];kept=[];seen=set()
            for i in src:
                p=d.loc[i,'protein_id']
                if p not in seen:seen.add(p);kept.append(i);d.loc[i,'group_id']=f'{fam}_g{center}'
            f=d.loc[kept];ph=Counter(int(x) for x in f.phase_before if int(x)>=0);pm=ph.most_common(1)[0][0] if ph else -1
            loci=Counter(f.locus_class.fillna('unplaced'));species=Counter(f.species)
            row=dict(family=fam,group_id=f'{fam}_g{center}',ref_aa_pos=center,n_genes=len(f),n_species=len(species),phase_mode=pm,phase_mode_n=ph[pm] if ph else 0,phase_mode_pct=round(100*ph[pm]/sum(ph.values()),1) if ph else 0,gt_ag_count=int(((f.donor_dinuc=='GT')&(f.acceptor_dinuc=='AG')).sum()),gt_ag_pct=round(100*((f.donor_dinuc=='GT')&(f.acceptor_dinuc=='AG')).sum()/len(f),1),locus_class_breakdown=json.dumps(dict(sorted(loci.items()))),major_group=str(len(f)>=5).lower())
            groups.append(row)
            if len(f)>=5:summary.append(row)
        d.to_csv(a.output_dir/f'{fam}_intron_positions.tsv',sep='\t',index=False);d.to_csv(a.output_dir/f'{fam}_positions_grouped.tsv',sep='\t',index=False);pd.DataFrame(groups).to_csv(a.output_dir/f'{fam}_groups.tsv',sep='\t',index=False)
        seqs=readfa(a.sequence_dir/f'{fam}_current_genome_proteins.faa');rich=[];ctx=[]
        for g in groups:
            f=d[d.group_id==g['group_id']].drop_duplicates('protein_id');wins=[]
            for _,r in f.iterrows():
                if r.protein_id in seqs:wins.append(window(seqs[r.protein_id],int(r.aa_position)))
            consensus=[]
            for col in range(2*FLANK):
                chars=[w[col] for w in wins if w[col]!='.'];top,n=(Counter(chars).most_common(1)[0] if chars else ('_',0));pct=100*n/len(chars) if chars else 0
                consensus.append(top);rich.append(dict(family=fam,group_id=g['group_id'],ref_aa_pos=g['ref_aa_pos'],col_offset=col-FLANK,top_residue=top,top_pct=round(pct,1),n_observed=len(chars)))
            c=''.join(consensus);ctx.append(dict(family=fam,group_id=g['group_id'],ref_aa_pos=g['ref_aa_pos'],n_genes=g['n_genes'],n_species=g['n_species'],n_members_with_context=len(wins),consensus_pretty=c[:FLANK]+'|'+c[FLANK:],n_species_cross_locus=0,cross_locus_species=''))
        pd.DataFrame(rich).to_csv(a.output_dir/f'{fam}_context_rich.tsv',sep='\t',index=False);pd.DataFrame(ctx).to_csv(a.output_dir/f'{fam}_context_summary.tsv',sep='\t',index=False)
        qc[fam]={'eligible_proteins':int(d.protein_id.nunique()),'mapped_intron_records':len(idx),'groups':len(groups),'recurrent_groups':[{k:r[k] for k in ('ref_aa_pos','n_genes','n_species','locus_class_breakdown')} for r in groups if r['n_genes']>=5]}
    pd.DataFrame(summary).to_csv(a.output_dir/'recurrent_position_summary.tsv',sep='\t',index=False);(a.output_dir/'QC.json').write_text(json.dumps(qc,indent=2)+'\n');print(json.dumps(qc,indent=2))
if __name__=='__main__':main()
