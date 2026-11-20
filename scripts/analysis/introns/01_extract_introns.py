#!/usr/bin/env python3
"""Rebuild intron metadata for the complete current Table S2 genome roster.

The accepted Table S2 metadata is the membership authority.  Each source group is
bound to one proteome/GFF/assembly triplet containing every accepted protein.  No
protein is silently omitted: the build fails unless all accepted IDs resolve to
both a protein sequence and CDS features.
"""
from __future__ import annotations
import argparse, csv, json, re
from collections import defaultdict
from pathlib import Path
from pyfaidx import Fasta

FAMILIES=('GH5_2','PL1_8','PL4_4')

def norm(s: str) -> str:
    s=s.strip()
    for prefix in ('transcript:','rna-','cds-'):
        if s.startswith(prefix): s=s[len(prefix):]
    return s[:-2] if s.endswith('.1') else s

def read_fasta(path):
    seqs={}; key=None; parts=[]
    with open(path,errors='replace') as h:
        for line in h:
            if line.startswith('>'):
                if key is not None: seqs[key]=''.join(parts)
                key=line[1:].split()[0]; parts=[]
            else: parts.append(line.strip())
    if key is not None: seqs[key]=''.join(parts)
    return seqs

def source_alias(source, pid):
    if source!='DL_bkragtag': return pid
    x=pid.replace('DL_bkragtag_','D_quercuslanigerum_',1)
    x=re.sub(r'(\.1)_[12]$',r'\1',x)
    return x

def load_source_manifest(path: Path):
    """Load the explicitly accepted protein/GFF/assembly triplet per source."""
    def resolve(value: str) -> Path:
        candidate=Path(value)
        return candidate if candidate.is_absolute() else (path.parent/candidate).resolve()
    with path.open() as h:
        rows=list(csv.DictReader(h,delimiter='\t'))
    required={'source_file','proteome','gff','assembly'}
    if not rows or not required.issubset(rows[0]):
        raise ValueError(f'{path}: expected columns {sorted(required)}')
    result={}
    for row in rows:
        source=row['source_file']
        if source in result: raise ValueError(f'{path}: duplicate source_file {source}')
        result[source]=(resolve(row['proteome']),resolve(row['gff']),resolve(row['assembly']))
    return result

def parse_cds(gff, alias_to_pid):
    out={}; invalid=set()
    with open(gff,errors='replace') as h:
        for line in h:
            if line.startswith('#'): continue
            p=line.rstrip().split('\t')
            if len(p)!=9 or p[2]!='CDS': continue
            attrs={}
            for bit in p[8].split(';'):
                if '=' in bit:
                    k,v=bit.split('=',1); attrs[k]=v
            matched=None
            for key in ('protein_id','Parent','ID'):
                for value in attrs.get(key,'').split(','):
                    if norm(value) in alias_to_pid:
                        matched=alias_to_pid[norm(value)]; break
                if matched: break
            if matched:
                if matched in invalid: continue
                phase=int(p[7]) if p[7] in ('0','1','2') else 0
                entry=out.setdefault(matched,{'scaffold':p[0],'strand':p[6],'cds':[]})
                if entry['scaffold']!=p[0] or entry['strand']!=p[6]:
                    invalid.add(matched); out.pop(matched,None); continue
                entry['cds'].append((int(p[3]),int(p[4]),phase))
    for x in out.values(): x['cds'].sort()
    return out,invalid

def reverse_complement(seq):
    return seq.upper().translate(str.maketrans('ACGTN','TGCAN'))[::-1]

def introns(data):
    cds=data['cds']; strand=data['strand']; total=sum(e-s+1 for s,e,_ in cds)
    ans=[]; cumulative=0
    for i,((s,e,phase),(ns,ne,nphase)) in enumerate(zip(cds,cds[1:]),1):
        cumulative += e-s+1
        if ns<=e+1: continue
        rel=cumulative/total
        ans.append(dict(intron_number_genomic=i,start=e+1,end=ns-1,length=ns-e-1,
                        phase_before=(phase+e-s+1)%3,phase_after=nphase,
                        relative_cds_position=rel))
    if strand=='-':
        n=len(ans)
        for x in ans:
            x['intron_number']=n-x['intron_number_genomic']+1
            x['relative_cds_position']=1-x['relative_cds_position']
    else:
        for x in ans:x['intron_number']=x['intron_number_genomic']
    aa_len=total//3
    for x in ans:
        x['relative_cds_position']=round(x['relative_cds_position'],4)
        x['aa_position']=round(x['relative_cds_position']*aa_len)
    return ans

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--metadata',type=Path,required=True,help='Accepted Table S2 protein metadata TSV.')
    ap.add_argument('--source-manifest',type=Path,required=True,help='Accepted source_file-to-proteome/GFF/assembly manifest.')
    ap.add_argument('--output-dir',type=Path,required=True)
    a=ap.parse_args()
    a.output_dir.mkdir(parents=True,exist_ok=False)
    with a.metadata.open() as h:
        rows=[r for r in csv.DictReader(h,delimiter='\t') if r['source']=='genome' and r['family'] in FAMILIES]
    accepted_sources=load_source_manifest(a.source_manifest)
    groups=defaultdict(list)
    for r in rows: groups[r['source_file']].append(r)
    missing_sources=sorted(set(groups)-set(accepted_sources))
    unused_sources=sorted(set(accepted_sources)-set(groups))
    if missing_sources or unused_sources:
        raise RuntimeError(f'source manifest mismatch: missing={missing_sources}; unused={unused_sources}')
    results=[]; logs=[]; seq_by_family=defaultdict(dict); sources=[]
    for source,rs in sorted(groups.items()):
        ids=[r['protein_id'] for r in rs]
        faa,gff,assembly=accepted_sources[source]
        for p in (faa,gff,assembly):
            if not p.is_file(): raise FileNotFoundError(p)
        seqs=read_fasta(faa)
        alias_to_pid={norm(source_alias(source,pid)):pid for pid in ids}
        selected={pid:seqs.get(source_alias(source,pid)) or seqs.get(norm(source_alias(source,pid))) for pid in ids}
        # Some FASTAs retain .1 while the accepted ID does not (and conversely).
        nseq={norm(k):v for k,v in seqs.items()}
        selected={pid:(seq or nseq.get(norm(source_alias(source,pid)))) for pid,seq in selected.items()}
        miss_seq=[pid for pid,s in selected.items() if not s]
        if miss_seq: raise RuntimeError(f'{source}: sequences unresolved: {miss_seq}')
        cds,invalid=parse_cds(gff,alias_to_pid)
        miss_cds=sorted(set(ids)-set(cds)-invalid)
        if miss_cds: raise RuntimeError(f'{source}: CDS unresolved ({len(miss_cds)}): {miss_cds[:10]}')
        # Inputs are immutable. Require the companion FASTA index instead of
        # creating files beside a released or shared genome assembly.
        genome=Fasta(str(assembly),rebuild=False)
        for r in rs:
            pid=r['protein_id']; seq_by_family[r['family']][pid]=selected[pid]
            if pid in invalid:
                logs.append(dict(protein_id=pid,family=r['family'],species=r['species'],source_file=source,status='INELIGIBLE',reason='CDS features span multiple scaffolds or strands'))
                continue
            data=cds[pid]; xs=introns(data)
            base=dict(protein_id=pid,family=r['family'],species=r['species'],source='genome',scaffold=data['scaffold'],strand=data['strand'])
            if not xs:
                results.append(base|dict(intron_count=0,intron_number=0,start=0,end=0,length=0,phase_before=-1,phase_after=-1,relative_cds_position=0.0,donor_dinuc='',acceptor_dinuc='',aa_position=0))
            for x in xs:
                seq=str(genome[data['scaffold']][x['start']-1:x['end']])
                if data['strand']=='-': seq=reverse_complement(seq)
                results.append(base|dict(intron_count=len(xs),intron_number=x['intron_number'],start=x['start'],end=x['end'],length=x['length'],phase_before=x['phase_before'],phase_after=x['phase_after'],relative_cds_position=x['relative_cds_position'],donor_dinuc=seq[:2].upper(),acceptor_dinuc=seq[-2:].upper(),aa_position=x['aa_position']))
            logs.append(dict(protein_id=pid,family=r['family'],species=r['species'],source_file=source,status='OK',reason=f'{len(xs)} introns from {len(data["cds"])} CDS exons'))
        sources.append(dict(source_file=source,n_proteins=len(rs),proteome=str(faa),gff=str(gff),assembly=str(assembly)))
    cols=['protein_id','family','species','source','scaffold','intron_count','intron_number','start','end','length','strand','phase_before','phase_after','relative_cds_position','donor_dinuc','acceptor_dinuc','aa_position']
    for fam in FAMILIES:
        rr=[x for x in results if x['family']==fam]
        with (a.output_dir/f'{fam}_intron_metadata.tsv').open('w',newline='') as h:
            w=csv.DictWriter(h,fieldnames=cols,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rr)
        with (a.output_dir/f'{fam}_current_genome_proteins.faa').open('w') as h:
            for pid,seq in sorted(seq_by_family[fam].items()): h.write(f'>{pid}\n{seq}\n')
    for name,data,cols2 in [('extraction_log.tsv',logs,list(logs[0])),('source_manifest.tsv',sources,list(sources[0]))]:
        with (a.output_dir/name).open('w',newline='') as h:
            w=csv.DictWriter(h,fieldnames=cols2,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(data)
    qc={'accepted_proteins':len(rows),'accounted_proteins':len(logs),'intron_eligible':sum(r['status']=='OK' for r in logs),'intron_ineligible':sum(r['status']=='INELIGIBLE' for r in logs),'families':{f:{'accepted':sum(r['family']==f for r in rows),'sequence_resolved':len(seq_by_family[f]),'intron_eligible':sum(r['family']==f and r['status']=='OK' for r in logs),'intron_ineligible':sum(r['family']==f and r['status']=='INELIGIBLE' for r in logs),'with_introns':len({r['protein_id'] for r in results if r['family']==f and r['intron_count']>0})} for f in FAMILIES}}
    if len(logs)!=len(rows): raise AssertionError(qc)
    (a.output_dir/'QC.json').write_text(json.dumps(qc,indent=2)+'\n'); print(json.dumps(qc,indent=2))

if __name__=='__main__': main()
