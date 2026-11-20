#!/usr/bin/env python3
"""Exercise every projection class and reject missing sources and re-entry."""
import csv,sys,tempfile
from pathlib import Path
from public_scope import project,ALIASES,REPO,read,write
sys.path.insert(0,str(REPO/'scripts/release'))
from public_scope_policy import validate
retired=[s for s in (REPO/'config/public_retired_inputs.txt').read_text().splitlines() if s and not s.startswith('#')]
def fixture(root):
 for rel in retired:
  p=root/rel;p.parent.mkdir(parents=True,exist_ok=True);p.write_text('id\nold\n')
 p=root/'tables/authority/tables_s5_s6_current/baseline/source_sample_roster.tsv';p.parent.mkdir(parents=True,exist_ok=True)
 rows=[[f'A{i:03}',f'A{i:03}'] for i in range(18)]+[[f'Z{i:03}',f'Z{i:03}'] for i in range(140)]+[['BIP__unsorted:BIP_adult_OVX','BIP_adult_OVX'],['BIP__unsorted:BIP_adult_VGX','BIP_adult_VGX']]
 write(p,['full_sample_id','original_sample_id'],rows)
 for name in ALIASES:
  p=root/'figures/origin_trees'/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_text('>protein\nAMINO\n')
 for name,header,rows in [('legacy.tsv',['mc_v1','cynogv4_id','cynog_id','value'],[['MC_1','CynOGv4_1','CynOGv3_1','1.23456789']]),('matrix.tsv',['sample','tpm'],[['retained','1.23456789'],['BIP_adult_OVX','9']])]:write(root/name,header,rows)
 write(root/'figures/INPUT_MANIFEST.tsv',['component','relative_path','role','size_bytes'],[['origin_trees','origin_trees/'+a,'alignment','0'] for a in ALIASES])
 write(root/'INPUT_MANIFEST.tsv',['relative_path','size_bytes','source_class'],[[str(p.relative_to(root)),str(p.stat().st_size),'fixture'] for p in root.rglob('*') if p.is_file()])
with tempfile.TemporaryDirectory() as d:
 root=Path(d);fixture(root);receipt=project(root);validate(root,REPO,check_pins=False)
 assert read(root/'legacy.tsv')==(['cynog_id','value'],[['CynOGv3_1','1.23456789']])
 assert read(root/'matrix.tsv')[1]==[['retained','1.23456789']]
 assert all((root/'figures/origin_trees'/b).read_bytes()==b'>protein\nAMINO\n' for b in ALIASES.values())
 assert len(read(root/'tables/authority/tables_s5_s6_current/baseline/source_sample_roster.tsv')[1])==158
 for header,row in [(['mc_v1'],['MC_1']),(['sample'],['BIP_adult_OVX']),(['sample_id'],['EXPSAMPLE-019'])]:
  p=root/'reentry.tsv';write(p,header,[row])
  try:validate(root,REPO,check_pins=False)
  except ValueError:pass
  else:raise AssertionError('excluded input accepted')
  p.unlink()
for absent in [retired[0],'figures/origin_trees/'+next(iter(ALIASES))]:
 with tempfile.TemporaryDirectory() as d:
  root=Path(d);fixture(root);(root/absent).unlink()
  try:project(root)
  except ValueError:pass
  else:raise AssertionError('missing input accepted')
print('PUBLIC_SCOPE_FIXTURE_PASS')
