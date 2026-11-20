#!/usr/bin/env python3
"""Select ten branch changes with explicit tip evidence for manual review."""
from __future__ import annotations
import argparse,csv
from collections import defaultdict
from pathlib import Path
from tree import parse_newick,prune
from score import CHROMS,score_character
def rows(path):
    with open(path) as handle:yield from csv.DictReader(handle,delimiter='\t')
p=argparse.ArgumentParser()
for name in ('scores','tips','changes','tree','out'):p.add_argument('--'+name,type=Path,required=True)
p.add_argument('--ecami',action='store_true')
a=p.parse_args();scores=list(rows(a.scores));states=defaultdict(dict)
unit='group' if a.ecami else 'orthogroup'
for r in rows(a.tips):
    if r['placed']=='1':states[r[unit]][r['species']]=set(r['chromosomes'].split(','))
focal=[r for r in rows(a.changes) if r['mode']=='deltran' and r['change']=='gain']
selected=[]
for name,count in [('FOCAL:GH5_2',2),('FOCAL:PL1_8',1),('FOCAL:PL4_4',2)]:
    subset=[r for r in focal if r[unit]==name]
    if len(subset)<count:raise ValueError(f'insufficient focal gains: {name}')
    selected.extend(subset[:count])
tree=parse_newick(a.tree.read_text());background=[r for r in scores if r['kind']==('native' if a.ecami else 'background') and int(r['deltran_gains'])>0]
for r in background[:5]:
    og=r[unit];sub=prune(tree,set(states[og]))
    for chrom in CHROMS:
        ones={tip for tip,labels in states[og].items() if chrom in labels}
        root,cost,gains,losses=score_character(sub,ones,'deltran')
        if gains:
            selected.append({unit:og,'chromosome':chrom,**gains[0]});break
assert len(selected)==10
out=[]
for r in selected:
    og=r[unit];chrom=r['chromosome'];inside=set(r['descendant_tips'].split(','))
    if not inside<=set(states[og]):raise ValueError((og,inside-set(states[og])))
    occupied={tip for tip,labels in states[og].items() if chrom in labels}
    out.append({unit:og,'chromosome':chrom,'parent':r['parent'],'child':r['child'],
                'descendants':len(inside),'inside_occupied':','.join(sorted(inside&occupied)),
                'inside_absent':','.join(sorted(inside-occupied)),
                'outside_occupied':','.join(sorted(occupied-inside)),
                'outside_absent_count':len(set(states[og])-inside-occupied)})
with a.out.open('w',newline='') as handle:
    cols=list(out[0]);w=csv.DictWriter(handle,cols,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(out)
print('HAND_CHECK_ROWS',len(out))
