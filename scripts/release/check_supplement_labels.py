#!/usr/bin/env python3
"""Check stable-ID label references in the final supplement render inputs."""
import argparse,csv,re
from pathlib import Path
from blectin_domain_counts import load_counts,verify_counts,verify_table_s3_caption
p=argparse.ArgumentParser();p.add_argument('--repo',type=Path,required=True);p.add_argument('--build',type=Path,required=True);a=p.parse_args()
def labels(path,key):
 with path.open(newline='') as h:return {r[key]:r['label'] for r in csv.DictReader(h,delimiter='\t')}
f=labels(a.repo/'config/figure_numbers.tsv','figure_id')
t=labels(a.repo/'config/table_numbers.tsv','table_id')
assert len(set(f.values()))==len(f)
assert len(set(t.values()))==len(t)
assert set(t.values())=={f'S{i}' for i in range(1,11)}
caption=(a.build/'table_s2/candidate/supplementary_tables/Table_S2_caption.txt').read_text()
assert caption.startswith(f"Table {t['outgroup']}.")
verify_table_s3_caption(caption,verify_counts(a.build/'figures/b_lectin_domain_tree/primary/trees/Figure_S4_B_lectin_domains.rooted.treefile',a.build/'figures/b_lectin_domain_tree/primary/sidecars/Figure_S4_domain_manifest.tsv',load_counts(a.build/'blectin_domain_tree_counts.json')),f['b_lectin_domain_tree'])
assert '69 Figure S3 tips' not in caption
for stable,old in [('loci',3),('synteny',4),('expression',6)]:
 candidates=list(a.build.glob(f'**/Table_S{old}_caption.txt'))
 assert candidates,stable
 assert any(p.read_text().startswith(f'Table {t[stable]}' + ('b.' if stable == 'expression' else '.')) for p in candidates),stable
assert (a.build/'table_s6/supplementary_tables/Table_S6a_caption.txt').read_text().startswith('Table S2a.')
assert (a.build/'tables/Table_S2a.pdf').is_file()
assert (a.build/'tables/Table_S2b.pdf').is_file()
assert f"\\section*{{Table {t['intron']}}}" in (a.build/f"tables/Table_{t['intron']}.tex").read_text()
tex=(a.build/'supplement/supplement.tex').read_text()
for stable in ('euceroptres_genome','tree_counts_full','origin_trees_full','accessory_domain_trees','complete_locus_trees','b_lectin_domain_tree'):
 assert f'figlabel{stable}' in tex or f'Figure {f[stable]}' in tex,stable
assert f'../tables/Table_{t["cazy_rank"]}.pdf' in tex
assert '../tables/Table_S2a.pdf' in tex and '../tables/Table_S2b.pdf' in tex
assert f'../tables/Table_{t["chromosome_gain_summary"]}.pdf' in tex
assert f'../tables/Table_{t["gc3"]}.pdf' in tex
assert not re.search(r'pending',tex,re.IGNORECASE)
assert 'figlabelchromosome_gain_distribution' not in tex and 'chromosome_gain/chromosome_gain.pdf' not in tex
assert 'figlabelchromosome_gain_trees' not in tex
for path in (a.repo/'config').glob('supplement_*.tex'):
 assert not re.search(r'\b(?:Table|Figure) S\d+\b',path.read_text()),path
for path in (a.build/'tables').glob('Table_S*.tsv'):
 with path.open(newline='') as handle:
  reader=csv.DictReader(handle,delimiter='\t')
  fields=reader.fieldnames or []
  assert not any(field in {'species_id','species_code','internal_species_code'} for field in fields),path
  for row in reader:
   for field in fields:
    if 'species' in field and field not in {'species_id','species_code','internal_species_code'}:
     assert '_' not in row[field],(path,field,row[field])
print('SUPPLEMENT_LABELS_PASS')
