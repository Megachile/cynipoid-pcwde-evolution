#!/usr/bin/env python3
"""Locate frozen accepted copies only; no block construction or state changes."""

import sys as _source_sys
from pathlib import Path as _SourcePath
_source_repo = next(p for p in _SourcePath(__file__).resolve().parents if (p / 'config/logical_inputs.tsv').is_file())
_source_sys.path.insert(0, str(_source_repo / 'scripts/release'))
from source_aliases import source_path
import csv, hashlib, json, subprocess, sys
from pathlib import Path
from collections import defaultdict
import os
R=Path(os.environ['MODELLESS_RUN_ROOT']).resolve()
P=Path(os.environ['MODELLESS_PROJECT_ROOT']).resolve()
C=R/'build/rebase_repo/config'
O=R/'00_REVIEW/modelless_recovery'
B=O/'pre_rule'
MINI=Path(os.environ['MODELLESS_TOOL_BIN'])/'miniprot'
def rows(p):
 with Path(p).open() as f:return list(csv.DictReader(f,delimiter='\t'))
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for x in iter(lambda:f.read(1048576),b''):h.update(x)
 return h.hexdigest()
def fasta(p):
 d={};k=None
 with Path(p).open() as f:
  for l in f:
   if l.startswith('>'):k=l[1:].split()[0];d[k]=''
   elif k:d[k]+=l.strip()
 return d
def write(p,rs):
 with Path(p).open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rs[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rs)
sources={x['source_file']:x for x in rows(C/'intron_sources.tsv')}
gmap={x['genome']:sources[x['source_file']] for x in rows(C/'intron_genome_sources.tsv')}
mem=rows(B/'member_table.tsv');ids={x['accepted_id']:x for x in rows(B/'synteny_existing_identity.tsv')}
old={x['accepted_id']:x for x in rows(R/'00_REVIEW/UNRESOLVED_18_R105.tsv')}
seqpath=source_path(P, 'analysis/pcwde_effectors/pcwde_paper/analysis/beta_uniform_transcript_membership/results/outputs/accepted_sequences.faa')
assert sha(seqpath)=='274eb8e8839d2edd89cddbecdb6467784dfa580eecfe4a47bfd5c8632255f761'
seq=fasta(seqpath)
# Enumerate exact source-annotation identity for every accepted member, all families.
bygen=defaultdict(list)
for m in mem:bygen[m['genome']].append(m)
annot={};audit=[]
for genome,ms in sorted(bygen.items()):
 s=gmap[genome]; gff=source_path(P, s['gff_project_path']); assert sha(gff)==s['gff_sha256']
 wanted={m['protein_id'] for m in ms}; hits=defaultdict(list)
 with gff.open() as f:
  for l in f:
   if l.startswith('#'):continue
   a=l.rstrip().split('\t')
   if len(a)!=9 or a[2] not in ('gene','mRNA','transcript','CDS'):continue
   attrs=dict(x.split('=',1) for x in a[8].split(';') if '=' in x)
   for key in ('ID','protein_id'):
    v=attrs.get(key,'')
    for pid in (v,v[:-2] if v.endswith('.1') else ''):
     if pid in wanted:hits[pid].append((a[0],a[3],a[4],a[2],v))
 for m in ms:
  pid=m['protein_id'];annot[pid]=hits[pid]
  audit.append(dict(genome=genome,family=m['family'],accepted_id=pid,order_id=ids[pid]['order_id'],annotation_exact_features=json.dumps(hits[pid]),annotation_gff=s['gff_project_path'],annotation_sha256=s['gff_sha256'],issue=ids[pid]['issue']))
write(O/'ALL_ACCEPTED_ANNOTATION_AUDIT.tsv',audit)
recovery=[];aligned=set()
for m in mem:
 pid=m['protein_id']; ident=ids[pid]
 if ident['order_id']:continue
 s=gmap[m['genome']]
 item=dict(genome=m['genome'],family=m['family'],accepted_id=pid,unit_id=m['unit_id'],accepted_bk_chromosome=m['bk_chromosome'],recorded_subject=m['unit_coordinate_subject'],recorded_start=m['unit_coordinate_start'],recorded_end=m['unit_coordinate_end'],assembly_project_path=s['assembly_project_path'],assembly_sha256=s['assembly_sha256'],source_gff_exact_features=json.dumps(annot[pid]),scaffold='',start='',end='',strand='',method='',identity='',query_coverage='',best_score='',n_best_loci='',status='')
 # The nine explicitly specified by the signed amendment require independent alignment.
 align=pid in old and old[pid]['source_interval']=='missing'
 if not align:
  item.update(scaffold=m['unit_coordinate_subject'],start=m['unit_coordinate_start'],end=m['unit_coordinate_end'],method='recorded_source_interval',status='resolved' if m['unit_coordinate_subject'] and m['unit_coordinate_start'] and m['unit_coordinate_end'] else 'unresolved_missing_interval')
 else:
  asm=source_path(P, s['assembly_project_path']); assert sha(asm)==s['assembly_sha256'];aligned.add(str(asm));assert len(aligned)<=18
  q=O/(pid+'.faa');q.write_text('>'+pid+'\n'+seq[m['unit_id']]+'\n')
  paf=O/(pid+'.paf');log=O/(pid+'.miniprot.log')
  with paf.open('w') as out, log.open('w') as err:
   subprocess.run([str(MINI),'-t','4','-N','100','--outn=100','--outs=0','--outc=0.1',str(asm),str(q)],stdout=out,stderr=err,check=True)
  hs=[]
  for l in paf.read_text().splitlines():
   if l.startswith('#'):continue
   a=l.split('\t')
   if len(a)<12:raise ValueError('invalid PAF '+l)
   tags={x.split(':',2)[0]:x.split(':',2)[2] for x in a[12:]}
   cov=(int(a[3])-int(a[2]))/int(a[1]);identity=int(a[9])/int(a[10])
   hs.append(dict(scaffold=a[5],start=int(a[7])+1,end=int(a[8]),strand=a[4],identity=identity,query_coverage=cov,best_score=float(tags['AS'])))
  best=max((h['best_score'] for h in hs),default=None)
  top=[h for h in hs if h['best_score']==best]
  unique={(h['scaffold'],h['start'],h['end'],h['strand']) for h in top}
  item['n_best_loci']=len(unique)
  if len(unique)==1 and top[0]['identity']>=.95 and top[0]['query_coverage']>=.90:
   item.update(top[0]);item.update(method='miniprot_best_AS',status='resolved')
  else:item.update(method='miniprot_best_AS',status='unresolved_tie_or_below_threshold',best_score=best if best is not None else '')
 recovery.append(item)
 write(O/'COORDINATE_RECOVERY.tsv',recovery)
 print(pid,item['status'],item['scaffold'],flush=True)
# The extra raw PL1_8 rescue is a fixed accepted presence copy outside member_table.
rescue=source_path(P, 'analysis/pcwde_effectors/pcwde_paper/analysis/raw_genome_tblastn_presence_loop/outputs/fig5_validated_raw_locus_rescue.tsv')
write(O/'RAW_LOCUS_RESCUE.tsv',rows(rescue))
# Read-only overlapping annotation evidence for the resolved loci, including BK.
overlap=[]
for genome in sorted({x['genome'] for x in recovery if x['status']=='resolved'}):
 wanted=[x for x in recovery if x['genome']==genome and x['status']=='resolved']
 with (source_path(P, gmap[genome]['gff_project_path'])).open() as f:
  for l in f:
   if l.startswith('#'):continue
   a=l.rstrip().split('\t')
   if len(a)!=9 or a[2] not in ('gene','mRNA','transcript','CDS'):continue
   for x in wanted:
    lo,hi=sorted((int(x['start']),int(x['end'])))
    if a[0]==x['scaffold'] and int(a[3])<=hi and int(a[4])>=lo:
     overlap.append(dict(accepted_id=x['accepted_id'],genome=genome,feature=a[2],scaffold=a[0],start=a[3],end=a[4],attributes=a[8]))
if overlap:write(O/'OVERLAPPING_ANNOTATION.tsv',overlap)
aa_hits=[]
for genome in sorted({x['genome'] for x in recovery}):
 s=gmap[genome];protein=source_path(P, s['proteome_project_path']);assert sha(protein)==s['proteome_sha256']
 proteins=fasta(protein)
 for x in recovery:
  if x['genome']!=genome:continue
  query=seq.get(x['unit_id'],'').rstrip('*')
  for pid,aa in proteins.items():
   if query and aa.rstrip('*')==query:aa_hits.append(dict(genome=genome,accepted_id=x['accepted_id'],exact_aa_model=pid))
if aa_hits:write(O/'EXACT_AA_ANNOTATION.tsv',aa_hits)
(O/'SUMMARY.json').write_text(json.dumps(dict(accepted_member_rows=len(mem),no_order_rows=len(recovery),aligned_assemblies=len(aligned),unresolved=sum(x['status']!='resolved' for x in recovery),raw_rescue_rows=len(rows(rescue))),indent=2)+'\n')
print('COORDINATE_RECOVERY_COMPLETE',flush=True)
