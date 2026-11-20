#!/usr/bin/env python3
"""Build uniform, tip-complete Figure S4 full-protein architecture from Pfam-A hmmscan.

Scan is fixed at HMMER 3.4, Pfam-A SHA256 6dd42235..., --domE 1
--incdomE 1. The historical S3 PF01453 tree-usable recovery requires
independent E <= 1, HMM coverage >= 0.5 and length >= 40 aa. Other Pfam
domains use conditional E < 1e-3. Overlaps > half the shorter aligned region
are resolved by conditional E value.
"""
import argparse,csv,hashlib
from collections import defaultdict
from pathlib import Path

def digest(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--mapping',required=True);ap.add_argument('--domtbl',required=True);ap.add_argument('--fasta',required=True);ap.add_argument('--out',required=True);a=ap.parse_args()
 mapping=list(csv.DictReader(open(a.mapping),delimiter='\t'))
 assert len(mapping)==246 and len({r['tip_id'] for r in mapping})==246
 assert all(r['state']=='mapped' and (not r['domain_start'] or 1<=int(r['domain_start'])<=int(r['domain_end'])<=int(r['full_length_aa'])) for r in mapping)
 seqs={};name=None
 for line in open(a.fasta):
  if line.startswith('>'):name=line[1:].strip();seqs[name]=''
  else:seqs[name]+=line.strip()
 assert len(seqs)==239
 for r in mapping:
  name='S4_'+r['sequence_sha256'][:20]
  assert name in seqs and digest_text(seqs[name])==r['sequence_sha256'] and len(seqs[name])==int(r['full_length_aa'])
 hits=defaultdict(list)
 for line in open(a.domtbl):
  if line.startswith('#') or not line.strip():continue
  f=line.split(maxsplit=22)
  if len(f)<22:raise ValueError('short domtbl row')
  name,pfam,qlen,query=f[0],f[1],int(f[5]),f[3]
  ce=float(f[11]); ie=float(f[12]); start,end=int(f[19]),int(f[20]); ali_start,ali_end=int(f[17]),int(f[18])
  if query not in seqs or qlen!=len(seqs[query]):raise ValueError('unknown query or length')
  if not 1<=start<=end<=qlen:raise ValueError('invalid protein coordinates')
  coverage=(int(f[16])-int(f[15])+1)/int(f[2])
  b_recovered=name=='B_lectin' and ie<=1 and coverage>=0.5 and end-start+1>=40
  if ce>=1e-3 and not b_recovered:continue
  if name=='B_lectin' and not b_recovered:continue
  hits[query].append(dict(domain=name,pfam=pfam,start=start,end=end,ali_start=ali_start,ali_end=ali_end,conditional_evalue=ce,independent_evalue=ie,weak_b='yes' if name=='B_lectin' and ie>1e-5 else 'no'))
 kept={}
 for query,raw in hits.items():
  accepted=[]
  for h in sorted(raw,key=lambda h:(h['conditional_evalue'],h['start'],h['end'],h['domain'])):
   bad=False
   for k in accepted:
    overlap=max(0,min(h['ali_end'],k['ali_end'])-max(h['ali_start'],k['ali_start'])+1)
    shorter=min(h['ali_end']-h['ali_start']+1,k['ali_end']-k['ali_start']+1)
    if overlap/shorter>0.5:bad=True;break
   if not bad:accepted.append(h)
  kept[query]=sorted(accepted,key=lambda h:(h['start'],h['end'],h['domain']))
 out=[]
 for r in mapping:
  query='S4_'+r['sequence_sha256'][:20]
  if not kept.get(query):raise ValueError('no domains for '+r['tip_id'])
  for h in kept[query]:out.append(dict(family=r['family'],tip_id=r['tip_id'],protein_length_aa=r['full_length_aa'],sequence_sha256=r['sequence_sha256'],**h))
 fields=list(out[0]);p=Path(a.out);p.parent.mkdir(parents=True,exist_ok=True)
 with p.open('w') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(out)
 print('tips',len({r['tip_id'] for r in out}),'proteins',len(seqs),'blocks',len(out),'weak B',sum(r['weak_b']=='yes' for r in out),'sha256',digest(p))

def digest_text(s):return hashlib.sha256(s.encode()).hexdigest()
if __name__=='__main__':main()
