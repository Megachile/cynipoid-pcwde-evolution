#!/usr/bin/env python3
"""Build Table S5 from the closed current-roster intron release."""
import argparse,csv,json
from collections import defaultdict
from pathlib import Path
FAMS=('GH5_2','PL1_8','PL4_4')
FIELDS=['part','family_id','eligible_proteins','proteins_with_introns','total_introns','protein_prevalence_pct','reference_position_aa','observation_count','observation_proportion_pct','phase','taxonomic_breadth','locus_distribution','taxonomic_families']
FAMILY_ORDER=('Diplolepididae','Figitidae','Cynipidae')
def read(p):
    with open(p) as h:return list(csv.DictReader(h,delimiter='\t'))
def fmt(x):return f'{x:.3f}'.rstrip('0').rstrip('.')
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--extraction-dir',type=Path,required=True)
    ap.add_argument('--summary-dir',type=Path,required=True)
    ap.add_argument('--b-lectin-extraction-dir',type=Path,required=True,help='Built 28-copy whole-gene B-lectin intron metadata.')
    ap.add_argument('--output-dir',type=Path,required=True)
    a=ap.parse_args();a.output_dir.mkdir(parents=True,exist_ok=False)
    taxonomy = read(Path(__file__).resolve().parents[3]/'config/cazy_lineage_groups.tsv')
    species_family = {r['species'].replace('_',' '): r['source_family'] for r in taxonomy}
    rows=[];den={}
    for fam in FAMS:
        raw=read(a.extraction_dir/f'{fam}_intron_metadata.tsv');ids={r['protein_id'] for r in raw};positive={r['protein_id'] for r in raw if int(r['intron_count'])>0};nobs=sum(int(r['intron_count'])>0 for r in raw);den[fam]=len(ids)
        rows.append(dict(part='A',family_id=fam,eligible_proteins=len(ids),proteins_with_introns=len(positive),total_introns=nobs,protein_prevalence_pct=fmt(100*len(positive)/len(ids)),reference_position_aa='',observation_count='',observation_proportion_pct='',phase='',taxonomic_breadth='',locus_distribution=''))
    braw=read(a.b_lectin_extraction_dir/'B_lectin_intron_metadata.tsv')
    bids={r['protein_id'] for r in braw};bpositive={r['protein_id'] for r in braw if int(r['intron_count'])>0}
    if len(bids)!=28:raise AssertionError(f'expected 28 B-lectin models, got {len(bids)}')
    bnobs=sum(int(r['intron_count'])>0 for r in braw)
    rows.append(dict(part='A',family_id='B_lectin',eligible_proteins=len(bids),proteins_with_introns=len(bpositive),total_introns=bnobs,protein_prevalence_pct=fmt(100*len(bpositive)/len(bids)),reference_position_aa='',observation_count='',observation_proportion_pct='',phase='',taxonomic_breadth='',locus_distribution=''))
    for fam in FAMS:
        groups=read(a.summary_dir/f'{fam}_groups.tsv')
        members=defaultdict(set)
        for position in read(a.summary_dir/f'{fam}_intron_positions.tsv'):
            if position['group_id']:
                species=position['species']
                family=species_family.get(species)
                if family not in FAMILY_ORDER:
                    raise AssertionError(f'no cynipoid family assignment for {species!r}')
                members[position['group_id']].add((position['protein_id'],species,family))
        for g in groups:
            n=int(g['n_genes'])
            if n<5:continue
            loci=json.loads(g['locus_class_breakdown'])
            dist='; '.join(f'{k}:{v}' for k,v in sorted(loci.items()))
            group_members=members[g['group_id']]
            if len({p for p,_,_ in group_members}) != n or len({s for _,s,_ in group_members}) != int(g['n_species']):
                raise AssertionError(f'intron group membership differs for {g["group_id"]}')
            families=', '.join(f for f in FAMILY_ORDER if f in {x[2] for x in group_members})
            if not families: raise AssertionError(f'empty family list for {g["group_id"]}')
            rows.append(dict(part='B',family_id=fam,eligible_proteins='',proteins_with_introns='',total_introns='',protein_prevalence_pct='',reference_position_aa=int(g['ref_aa_pos']),observation_count=n,observation_proportion_pct=fmt(100*n/den[fam]),phase=f"{g['phase_mode']} ({int(g['phase_mode_n'])}/{n})",taxonomic_breadth=f"{g['n_species']} species",locus_distribution=dist,taxonomic_families=families))
    with (a.output_dir/'Table_S5.tsv').open('w',newline='') as h:
        w=csv.DictWriter(h,fieldnames=FIELDS,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    caption=('Table S5. Intron prevalence and recurrent positions in current family members with coherent single-scaffold genomic CDS models. Part A includes all eligible models, including intronless copies, and reports proteins with at least one intron and individual intron observations. Part B reports positions represented in at least five distinct proteins after mapping full-length proteins to the family reference and clustering coordinates within three amino acids; each percentage uses all eligible family models as its denominator. A copy can contribute to more than one position. Phase is the modal standard intron phase (0 between codons, 1 after the first base, 2 after the second), with supporting observations shown as n/N.')
    (a.output_dir/'Table_S5_caption.txt').write_text(caption+'\n')
    print((a.output_dir/'Table_S5.tsv').read_text())
if __name__=='__main__':main()
