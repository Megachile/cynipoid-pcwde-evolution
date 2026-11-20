#!/usr/bin/env python3
"""Record the 41-genome CynOGv3 panel and source-tree coverage."""
import argparse,csv,json,re
from collections import Counter,defaultdict
from pathlib import Path
p=argparse.ArgumentParser()
for name in ('membership','catalog','tree','panel','out'):p.add_argument('--'+name,type=Path,required=True)
a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
def rows(path):
    with path.open() as handle:yield from csv.DictReader(handle,delimiter='\t')
tips=set(re.findall(r'(?<=[(,])([^():,;]+)(?=[:,)])',a.tree.read_text()))
panel={r['genome']:r['species'].replace(' ','_') for r in rows(a.panel)}
catalog=defaultdict(Counter)
for r in rows(a.catalog):
    catalog[r['genome']]['total']+=1
    if r['bk_chr'] not in ('','unplaced'):catalog[r['genome']]['catalog_placed']+=1
members=defaultdict(Counter);groups=set()
for r in rows(a.membership):members[r['genome']]['total']+=1;groups.add(r['cynog_id'])
assert len(members)==len(catalog)==len(panel)==41
assert set(members)==set(panel) and set(panel.values())<=tips
with (a.out/'input_panel.tsv').open('w',newline='') as handle:
    w=csv.writer(handle,delimiter='\t',lineterminator='\n')
    w.writerow(['genome','species','member_proteins','catalog_proteins','catalog_placed'])
    for g in sorted(panel):w.writerow([g,panel[g],members[g]['total'],catalog[g]['total'],catalog[g]['catalog_placed']])
summary={'member_genomes':41,'source_tree_tips':len(tips),'cynog_groups':len(groups),
         'member_proteins':sum(v['total'] for v in members.values()),'panel_species':sorted(panel.values())}
(a.out/'input_preflight.json').write_text(json.dumps(summary,indent=2)+'\n')
print('INPUT_PREFLIGHT_PASS',summary['member_proteins'],summary['cynog_groups'])
