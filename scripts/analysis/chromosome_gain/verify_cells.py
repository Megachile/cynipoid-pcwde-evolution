#!/usr/bin/env python3
"""Require focal protein counts to equal every Figure 1 chromosome cell."""
import argparse,csv,json
from collections import Counter
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--focal',type=Path,required=True)
p.add_argument('--figure-cells',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
a=p.parse_args()
rename={'Ceroptres_kovalevi':'Ceroptres_masudai','Eschatocerus_acaciae':'Eschatocerus_acacia','Eubothrus_centaureae':'Isocolus_centaureae'}
count=Counter()
with a.focal.open() as handle:
    for r in csv.DictReader(handle,delimiter='\t'):
        sp=r['species'].replace(' ','_');sp=rename.get(sp,sp)
        fam='CBM13' if r['family']=='B-lectin' else r['family']
        count[(sp,fam,r['bk_chromosome'])]+=len(r['protein_ids'].split(','))
tested=0;errors=[]
with a.figure_cells.open() as handle:
    for r in csv.DictReader(handle,delimiter='\t'):
        if not r['position'].startswith('chr'):continue
        tested+=1;key=r['species'],r['family'],r['position']
        if count[key]!=int(r['value']):errors.append((key,int(r['value']),count[key]))
if tested!=1248 or errors:raise ValueError(f'Figure 1 parity failed: {tested} cells, {errors[:5]}')
a.out.write_text(json.dumps({'figure1_chromosome_cells':tested,'mismatches':0},indent=2)+'\n')
print('FIGURE1_CELLS_PASS',tested)
