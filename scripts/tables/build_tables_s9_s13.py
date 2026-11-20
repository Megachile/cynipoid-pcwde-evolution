#!/usr/bin/env python3
"""Prepare lossless or RESULT-reconciled summary data; no biological recomputation."""
import csv,json,statistics,hashlib,re
from pathlib import Path
import argparse
parser=argparse.ArgumentParser();parser.add_argument('--data-root',type=Path,required=True);parser.add_argument('--out',type=Path,required=True);parser.add_argument('--analysis',type=Path,required=True);parser.add_argument('--catalytic-summary',type=Path,required=True);args=parser.parse_args()
ROOT=args.data_root/'repository_inputs/tables/authority/supplement_new'
OUT=args.out; OUT.mkdir(parents=True,exist_ok=True)
FAMILIES=['GH5_2','PL1_8','PL4_4','B-lectin']
NAMES={'CATALYTIC_CLUSTERING.tsv':9,'LOCUS_CONFLICT.tsv':12,'DI_SORT_LIMITS.tsv':13}
def get(run,file):return ROOT/run/file
def tsv(p):return list(csv.DictReader(p.open(),delimiter='\t'))
PREVIEW_SHA={9:'3d7762f0eb0606d47c7a54bcc7a6158bd511bd8c8408e593d4a8a9fc106bbb57',12:'6ad54be1fee012896fb393120c142f296244bfbdeeebf16ffe9cc3824b229a86',13:'f795008378c0199b2f3af7d5c0b42b24c1684562b9d23e6d8a02916fe5d0c779'}
CAPTIONS={
9:'Catalytic-residue clustering tests. Retained gives mapped residues over callable residues; exceptions are substitutions; observed changes are Fitch counts. Null mean and clustering P use 10,000 permutations.',
12:'Locus-tree conflict tests. Delta log likelihood and AU P compare stated alternatives; site concordance and informative-site ranges are shown for observed and alternative branches; long-branch sensitivity is separate.',
13:'DI inquiline sorting limits by family. Own-bin, confirmed-leak and unresolved are read-rule outcomes; possible leakage intervals include every unresolved tip. The lineage-qualified upper bound on possible misassignment is {accepted_bound} ({report_bound} before one PL4_4 tip assigned own-bin by the rule, but with only nearest-lineage support, was excluded as organismal evidence).'}
def write(name,rows,fields=None):
 if not rows:raise ValueError(name)
 fields=fields or list(rows[0]);number=NAMES[name];p=OUT/f'Table_S{number}.tsv'
 with p.open('w') as h:
  w=csv.DictWriter(h,fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
 digest=hashlib.sha256(p.read_bytes()).hexdigest()
 if digest!=PREVIEW_SHA[number]:raise SystemExit(f'Table S{number} differs from 2026-10-28 field regression: {digest}')
 (OUT/f'Table_S{number}_caption.txt').write_text(CAPTIONS[number]+'\n')
 print('Table_S'+str(number)+'.tsv',len(rows),'PREVIEW_PASS')
# Catalytic change clustering: every family/site test from accepted JSON.
a=json.loads(args.analysis.read_text())
summary={(r['family'],r['site_label']):r for r in tsv(args.catalytic_summary) if r['alignment_set']=='de_novo' and r['site_role']=='catalytic'}
for fam in FAMILIES[:3]:
 for site,x in a['families'][fam]['sites'].items():
  y=summary[(fam,site)]
  got=(x['mapped_retained'],x['mapped_total']); expected=(int(y['n_expected_residue']),int(y['n_site_observed']))
  if got!=expected: raise SystemExit(f'S9/S10 count disagreement {fam}/{site}: {got} != {expected}')
print('S9_S10_COUNT_AGREEMENT_PASS',len(summary))
rows=[]
for fam in FAMILIES[:3]:
 for site,x in a['families'][fam]['sites'].items():
  rows.append(dict(family=fam,site=site,retained=f"{x['mapped_retained']}/{x['mapped_total']}",exceptions=x['exceptions'],observed_changes=x['observed_changes'],null_mean=x['null_mean'],permutations=x['permutations'],p_clustered=x['p_clustered']))
write('CATALYTIC_CLUSTERING.tsv',rows)
# Locus conflict: all 17 accepted AU alternatives, sCF range and L1 in matching units.
run='locus_conflict';au=tsv(get(run,'AU_RESULTS.tsv'));scf=tsv(get(run,'SCF_RESULTS.tsv'));lb=tsv(get(run,'LONG_BRANCH_SENSITIVITY.tsv'))
lbmap={(x['case_id'],x['family']):x for x in lb};rows=[]
for x in au:
 if x['test'] in ('control','positive_control'):continue
 key=(x['test'],x['family']); observed=[y for y in scf if (y['test'],y['family'])==key and y['branch_role']=='observed']; alternative=[y for y in scf if (y['test'],y['family'])==key and y['branch_role']=='alternative']
 def rng(ys,col):
  vv=[float(y[col]) for y in ys if y[col] not in ('','NA')];return '' if not vv else str(min(vv)) if min(vv)==max(vv) else f'{min(vv)}–{max(vv)}'
 l=lbmap.get(key,{})
 rows.append(dict(test=x['test'],family=x['family'],delta_logL=x['deltaL'],p_AU=x['p_AU'],observed_sCF_pct=rng(observed,'sCF'),observed_sN=rng(observed,'sN'),alternative_sCF_pct=rng(alternative,'sCF'),alternative_sN=rng(alternative,'sN'),long_branch_result=l.get('l1_verdict','not run'),long_branch_UF=l.get('max_violating_l1_ufboot','')))
assert len(rows)==17,len(rows)
write('LOCUS_CONFLICT.tsv',rows)
# DI audit: parse the accepted report's primary-rule table, and require the
# coordinator's lineage-bound correction before reporting either bound.
report=get('di_sort_limits','SORT_LIMITS.md').read_text()
note=get('di_sort_limits','COORDINATOR_ACCEPTANCE_NOTE.md').read_text()
# Lineage-qualified bound: one shared intermediate (DI_LINEAGE_BOUND.tsv) feeds
# the caption and numbers.json. The report states (tips - attributed)/tips; the
# acceptance note lowers the strongly attributed count and states the result.
m=re.search(r'lineage-qualified \((\d+) − (\d+)\)/(\d+) = \*\*(\d+)/(\d+)\*\*',report)
if not m:raise SystemExit('DI report lineage-qualified bound statement missing')
tips,report_attr,den,report_num,den2=map(int,m.groups())
if tips!=den or den!=den2 or tips-report_attr!=report_num:raise SystemExit('DI report lineage-qualified bound arithmetic disagrees')
m=re.search(r'Treat strongly attributed tips as (\d+)/(\d+); the lineage-qualified upper bound is then (\d+)/(\d+), not (\d+)/(\d+)',note)
if not m:raise SystemExit('DI acceptance-note lineage bound statement missing')
note_attr,note_den,note_num,note_den2,old_num,old_den=map(int,m.groups())
if len({tips,note_den,note_den2,old_den})!=1 or note_den-note_attr!=note_num or (old_num,old_den)!=(report_num,den):
 raise SystemExit('DI acceptance-note bound disagrees with report or its own arithmetic')
bound=dict(tips=tips,report_strongly_attributed=report_attr,report_bound_numerator=report_num,accepted_strongly_attributed=note_attr,accepted_bound_numerator=note_num,denominator=tips,authority='COORDINATOR_ACCEPTANCE_NOTE.md',authority_sha256=hashlib.sha256(get('di_sort_limits','COORDINATOR_ACCEPTANCE_NOTE.md').read_bytes()).hexdigest())
with (OUT/'DI_LINEAGE_BOUND.tsv').open('w') as h:
 w=csv.DictWriter(h,list(bound),delimiter='\t',lineterminator='\n');w.writeheader();w.writerow(bound)
CAPTIONS[13]=CAPTIONS[13].format(accepted_bound=f'{note_num}/{tips}',report_bound=f'{report_num}/{tips}')
match=re.search(r'\| Family \| Tips \| Own bin supported \| Leaked \| Unresolved \|\n\|[-|: ]+\|\n((?:\|[^\n]+\|\n){4})',report)
if not match:raise SystemExit('DI primary-rule table missing')
rows=[]
for line in match.group(1).splitlines():
 fam,n,s,l,u=[v.strip().replace('**','') for v in line.strip('|').split('|')]
 n,s,l,u=map(int,(n,s,l,u))
 rows.append(dict(family=fam,tips=n,own_bin=s,confirmed_leak=l,unresolved=u,possible_leak_interval=f'0–{u}/{n}'))
if [(r['family'],r['tips'],r['own_bin'],r['confirmed_leak'],r['unresolved']) for r in rows] != [('GH5_2',7,4,0,3),('PL1_8',12,5,0,7),('PL4_4',15,5,0,10),('Total',34,14,0,20)]:
 raise SystemExit('DI primary-rule source counts changed')
if rows[-1]['tips']!=tips:raise SystemExit('DI lineage bound denominator differs from read-rule total')
write('DI_SORT_LIMITS.tsv',rows)
