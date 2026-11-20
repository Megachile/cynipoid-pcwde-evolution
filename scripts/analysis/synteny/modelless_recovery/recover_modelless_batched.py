#!/usr/bin/env python3
"""Resume only five unfinished accepted queries; serial bounded BLAST subjects."""

import sys as _source_sys
from pathlib import Path as _SourcePath
_source_repo = next(p for p in _SourcePath(__file__).resolve().parents if (p / 'config/logical_inputs.tsv').is_file())
_source_sys.path.insert(0, str(_source_repo / 'scripts/release'))
from source_aliases import source_path
import csv, hashlib, json, os, subprocess, time
from pathlib import Path
from collections import defaultdict
from Bio import SeqIO

import os
R=Path(os.environ['MODELLESS_RUN_ROOT']).resolve()
P=Path(os.environ['MODELLESS_PROJECT_ROOT']).resolve()
O=R/'00_REVIEW/modelless_recovery'
C=R/'build/rebase_repo/config'
BIN=Path(os.environ['MODELLESS_TOOL_BIN'])
T=O/'batch_scratch'
T.mkdir(exist_ok=True)
def rows(p):
 with Path(p).open() as f:return list(csv.DictReader(f,delimiter='\t'))
def write(p,rs):
 with Path(p).open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rs[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rs)
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
sources={x['source_file']:x for x in rows(C/'intron_sources.tsv')}
gmap={x['genome']:sources[x['source_file']] for x in rows(C/'intron_genome_sources.tsv')}
mem=rows(O/'pre_rule/member_table.tsv')
old={x['accepted_id']:x for x in rows(R/'00_REVIEW/UNRESOLVED_18_R105.tsv')}
retained=rows(O/'COORDINATE_RECOVERY.tsv')
done={x['accepted_id']:x for x in retained}
pending=[m for m in mem if m['protein_id'] in old and old[m['protein_id']]['source_interval']=='missing' and m['protein_id'] not in done]
assert len(pending)==5
pending.sort(key=lambda m:(source_path(P, gmap[m['genome']]['assembly_project_path'])).stat().st_size)
seqpath=source_path(P, 'analysis/pcwde_effectors/pcwde_paper/analysis/beta_uniform_transcript_membership/results/outputs/accepted_sequences.faa')
assert sha(seqpath)=='274eb8e8839d2edd89cddbecdb6467784dfa580eecfe4a47bfd5c8632255f761'
seq={r.id:str(r.seq) for r in SeqIO.parse(seqpath,'fasta')}
receipt=[]; new=[]; pilot_rss=None
for m in pending:
 pid=m['protein_id']; s=gmap[m['genome']]; asm=source_path(P, s['assembly_project_path'])
 assert sha(asm)==s['assembly_sha256']
 query=seq[m['unit_id']]; q=O/(pid+'.faa');q.write_text('>'+pid+'\n'+query+'\n')
 hits={}; mappings={}; batch=[]; bases=0; n=0; scanned=0
 raw=O/(pid+'.tblastn.tsv')
 raw.write_text('')
 def search():
  global batch,bases,n,pilot_rss
  if not batch:return
  subject=T/'subject.fna'; output=T/'hits.tsv'; metrics=T/'time.txt'
  subject.write_text(''.join('>'+k+'\n'+v+'\n' for k,v in batch))
  t=time.time()
  # Subject mode has no database/index files. All HSPs are retained for global ranking.
  cmd=[str(BIN/'tblastn'),'-query',str(q),'-subject',str(subject),'-seg','no','-evalue','1e-5','-num_threads','4','-max_target_seqs','1000000','-outfmt','6 sseqid qstart qend sstart send nident length bitscore evalue','-out',str(output)]
  subprocess.run(['/usr/bin/time','-f','%M','-o',str(metrics)]+cmd,check=True)
  rss=int(metrics.read_text().strip()); assert rss<12*1024*1024
  if pilot_rss is None:pilot_rss=rss
  for line in output.read_text().splitlines():
   a=line.split('\t'); scaffold,offset=mappings[a[0]]
   qs,qe,ss,se,ni,al=int(a[1]),int(a[2]),int(a[3])+offset,int(a[4])+offset,int(a[5]),int(a[6])
   h=(scaffold,qs,qe,ss,se,ni,al,float(a[7]),float(a[8]))
   hits[h[:5]]=h
   with raw.open('a') as f:f.write('\t'.join(map(str,h))+'\n')
  receipt.append(dict(copy=pid,batch=n,bases=bases,subjects=len(batch),maxrss_kib=rss,seconds=round(time.time()-t,3),database_bytes=0,index_bytes=0))
  write(O/'BATCH_RESOURCES.tsv',receipt)
  print(pid,'batch',n,'bases',bases,'rss_kib',rss,flush=True)
  batch=[];bases=0;n+=1
 # Split with 200 kb overlap, exceeding the 100 kb permitted HSP-chain gap.
 # Pilot is one <=2 Mb window from the smallest remaining assembly. Later
 # subject batches <=10 Mb only after measured pilot memory has ample headroom.
 for rec in SeqIO.parse(asm,'fasta'):
  scanned+=len(rec)
  for start in range(0,len(rec),1800000):
   window=str(rec.seq[start:start+2000000]); key='w'+str(len(mappings))
   mappings[key]=(rec.id,start)
   limit=2000000 if pilot_rss is None else min(10000000,max(2000000,int(2000000*(8*1024*1024)/max(pilot_rss,1))))
   if batch and bases+len(window)>limit:search()
   batch.append((key,window));bases+=len(window)
   if pilot_rss is None:search()
 search()
 assert scanned>0
 # Collinear, nonoverlapping query HSP chains, <=100 kb genomic gap.
 # DP maximizes summed bit score, with no identity/coverage filtering until
 # global ranking. Alternative loci and ties are retained, not batch winners.
 groups=defaultdict(list)
 for h in hits.values():groups[(h[0],1 if h[4]>=h[3] else -1)].append(h)
 candidates=[]
 for (scaffold,strand),hs in groups.items():
  hs.sort(key=lambda h:(h[1],h[2],h[3]))
  dp=[]
  for i,h in enumerate(hs):
   best=(h[7],[h])
   for j,prev in enumerate(hs[:i]):
    gap=(h[3]-prev[4])*strand
    if prev[2]<h[1] and 0<gap<=100000:
     option=(dp[j][0]+h[7],dp[j][1]+[h])
     if option[0]>best[0]:best=option
   dp.append(best)
  for score,chain in dp:
   candidates.append(dict(scaffold=scaffold,start=min(min(h[3:5]) for h in chain),end=max(max(h[3:5]) for h in chain),strand='+' if strand==1 else '-',identity=sum(h[5] for h in chain)/sum(h[6] for h in chain),query_coverage=sum(h[2]-h[1]+1 for h in chain)/len(query),best_score=score,hsps=len(chain)))
 candidates.sort(key=lambda h:(-h['best_score'],h['scaffold'],h['start'],h['end']))
 if candidates:write(O/(pid+'.candidates.tsv'),candidates)
 top=[h for h in candidates if h['best_score']==candidates[0]['best_score']]
 # Duplicate overlapping windows have already been deduplicated by exact HSP.
 unique={(h['scaffold'],h['start'],h['end'],h['strand']) for h in top}
 item=dict(genome=m['genome'],family=m['family'],accepted_id=pid,unit_id=m['unit_id'],accepted_bk_chromosome=m['bk_chromosome'],recorded_subject=m['unit_coordinate_subject'],recorded_start=m['unit_coordinate_start'],recorded_end=m['unit_coordinate_end'],assembly_project_path=s['assembly_project_path'],assembly_sha256=s['assembly_sha256'],source_gff_exact_features='',scaffold='',start='',end='',strand='',method='serial_subject_tblastn_collinear_global_bitscore',identity='',query_coverage='',best_score='',n_best_loci=len(unique),status='unresolved_tie_or_below_threshold')
 if len(unique)==1 and top[0]['identity']>=.95 and top[0]['query_coverage']>=.90:
  item.update({k:v for k,v in top[0].items() if k!='hsps'});item['status']='resolved'
 elif top:item['best_score']=top[0]['best_score']
 new.append(item);write(O/'BATCHED_COORDINATE_RECOVERY.tsv',new)
 print(pid,item['status'],flush=True)
 # Record completion before any interrupted annotation audit; never rerun blindly.
 (O/(pid+'.complete.json')).write_text(json.dumps(dict(assembly_sha256=s['assembly_sha256'],assembly_bases=scanned,unique_hsps=len(hits),batches=n,result=item),indent=2)+'\n')
print('BATCHED_SEARCHES_COMPLETE',flush=True)
# Finish the interrupted read-only annotation audit, retaining original methods.
ids={x['accepted_id']:x for x in rows(O/'pre_rule/synteny_existing_identity.tsv')}
all_done={x['accepted_id']:x for x in retained+new}
recovery=[]
for m in mem:
 pid=m['protein_id']
 if ids[pid]['order_id']:continue
 if pid in all_done:recovery.append(all_done[pid]);continue
 s=gmap[m['genome']]
 item=dict(retained[0])
 for k in item:item[k]=''
 item.update(genome=m['genome'],family=m['family'],accepted_id=pid,unit_id=m['unit_id'],accepted_bk_chromosome=m['bk_chromosome'],recorded_subject=m['unit_coordinate_subject'],recorded_start=m['unit_coordinate_start'],recorded_end=m['unit_coordinate_end'],assembly_project_path=s['assembly_project_path'],assembly_sha256=s['assembly_sha256'],scaffold=m['unit_coordinate_subject'],start=m['unit_coordinate_start'],end=m['unit_coordinate_end'],method='recorded_source_interval',status='resolved' if all(m[k] for k in ('unit_coordinate_subject','unit_coordinate_start','unit_coordinate_end')) else 'unresolved_missing_interval')
 recovery.append(item)
write(O/'COMPLETE_COORDINATE_RECOVERY.tsv',recovery)
def fasta(p):return {r.id:str(r.seq) for r in SeqIO.parse(p,'fasta')}
aligned={x['assembly_project_path'] for x in recovery if x['method']!='recorded_source_interval'}
original=(Path(__file__).resolve().parent/'recover_modelless.py').read_text()
exec(compile(original[original.index('# The extra raw PL1_8 rescue'):],str(Path(__file__).resolve().parent/'recover_modelless.py'),'exec'))
