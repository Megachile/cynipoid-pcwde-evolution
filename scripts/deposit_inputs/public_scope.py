#!/usr/bin/env python3
"""Project the public input layer; preserve retained cells and FASTA bytes."""
import csv,hashlib,json
from pathlib import Path
REPO=Path(__file__).resolve().parents[2]
EXCLUDED={'BIP_adult_OVX','BIP_adult_VGX'}
ALIASES={'MC_0006_1_s1_domains.faa':'CynOGv3_000234_domains.faa','MC_0030_1_s1_domains.faa':'CynOGv3_000148_domains.faa','MC_0047_4_s1_domains.faa':'CynOGv3_000168_domains.faa'}
LEGACY={'mc_v1','cynogv4_id'}
def read(p):
 with p.open(newline='') as f:
  rd=csv.reader(f,delimiter='\t');return next(rd),list(rd)
def write(p,h,rs):
 with p.open('w',newline='') as f:
  w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(h);w.writerows(rs)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def project(inputs):
 inputs=Path(inputs)
 # The explicitly approved retirements are pinned in the committed method list.
 retired=[s for s in (REPO/'config/public_retired_inputs.txt').read_text().splitlines() if s and not s.startswith('#')]
 if len(retired)!=31 or len(set(retired))!=31:raise ValueError('public retirement inventory differs')
 receipt=[]
 for rel in retired:
  p=inputs/rel
  if not p.is_file() or p.is_symlink():raise ValueError(f'missing retirement source: {rel}')
  receipt.append({'path':rel,'operation':'retire_unused_input','old_sha256':sha(p),'new_sha256':None});p.unlink()
 roster=inputs/'tables/authority/tables_s5_s6_current/baseline/source_sample_roster.tsv'
 h,rs=read(roster); idx=h.index('full_sample_id');ids={row[idx]:f'EXPSAMPLE-{i:03d}' for i,row in enumerate(sorted(rs,key=lambda x:x[idx]),1)}
 if len(ids)!=160:raise ValueError('historical sample assignment inventory differs')
 excluded_ids={v for k,v in ids.items() if k.split(':')[-1] in EXCLUDED}
 if len(excluded_ids)!=2:raise ValueError('excluded sample assignment inventory differs')
 for p in sorted(inputs.rglob('*.tsv')):
  if p.name=='INPUT_MANIFEST.tsv':continue
  h,rs=read(p);keep=[i for i,k in enumerate(h) if k not in LEGACY];new=[];removed=0
  for row in rs:
   if len(row)!=len(h):raise ValueError(f'nonrectangular input: {p}')
   bad=any(any(token in c for token in EXCLUDED) for c in row)
   bad=bad or any(row[i] in excluded_ids for i,k in enumerate(h) if k=='sample_id')
   if bad:removed+=1;continue
   new.append([row[i] for i in keep])
  nh=[h[i] for i in keep]
  if p==roster:
   nh.append('sample_id');new=[row+[ids[row[idx]]] for row in new]
  if nh!=h or removed:
   old=sha(p);write(p,nh,new);receipt.append({'path':p.relative_to(inputs).as_posix(),'operation':'project_fields_and_samples','old_sha256':old,'new_sha256':sha(p),'old_rows':len(rs),'new_rows':len(new),'removed_fields':sorted(set(h)-set(nh))})
 for old,new in ALIASES.items():
  src=inputs/'figures/origin_trees'/old;dst=src.with_name(new)
  if not src.is_file() or dst.exists():raise ValueError('missing or conflicting public FASTA alias')
  digest=sha(src);src.rename(dst)
  if sha(dst)!=digest:raise ValueError('FASTA alias changed sequence bytes')
  receipt.append({'path':src.relative_to(inputs).as_posix(),'operation':'sequence_alias','public_path':dst.relative_to(inputs).as_posix(),'old_sha256':digest,'new_sha256':digest})
 # Preserve order and roles of every retained inventory row.
 for rel in ['figures/INPUT_MANIFEST.tsv','INPUT_MANIFEST.tsv']:
  p=inputs/rel;h,rs=read(p);pi=h.index('relative_path');si=h.index('size_bytes');prefix=inputs/'figures' if rel.startswith('figures/') else inputs;new=[]
  for row in rs:
   path=Path(row[pi]);name=ALIASES.get(path.name,path.name);row[pi]=str(path.with_name(name));target=prefix/row[pi]
   if not target.is_file():continue
   row[si]=str(target.stat().st_size);new.append(row)
  write(p,h,new)
 return receipt
if __name__=='__main__':
 import argparse
 a=argparse.ArgumentParser();a.add_argument('inputs',type=Path);a.add_argument('--receipt',type=Path,required=True);args=a.parse_args();args.receipt.write_text(json.dumps(project(args.inputs),indent=2)+'\n')
