#!/usr/bin/env python3
"""Assign paper proteins with the shared per-protein placement rule."""
import argparse,csv,importlib.util,json
from collections import defaultdict,Counter
from pathlib import Path
parser=argparse.ArgumentParser()
parser.add_argument('--data-root',type=Path,required=True)
parser.add_argument('--out',type=Path,required=True)
parser.add_argument('--protein-genes',type=Path)
parser.add_argument('--scan-manifest',type=Path)
parser.add_argument('--panel',type=Path)
parser.add_argument('--bridge',type=Path)
args=parser.parse_args()
R=args.data_root/'repository_inputs/chromosome_gain'
spec=importlib.util.spec_from_file_location('paper_placement',Path(__file__).resolve().parents[2]/'analysis/synteny/placement.py')
placement=importlib.util.module_from_spec(spec);spec.loader.exec_module(placement)
def tsv(path):
    with path.open() as f:yield from csv.DictReader(f,delimiter='\t')
def bed(path):
    out={}
    with path.open() as f:
        for line in f:
            s,a,b,key=line.rstrip('\n').split('\t')[:4]
            out[key]=(s,int(a),int(b))
    return out
def segments(simple,left,right,left_ix,right_ix):
    # Each JCVI row has two gene endpoints per genome, then score and orientation.
    # The paired endpoints delimit one syntenic span in each genome.
    seg=[];partner=defaultdict(set)
    with simple.open() as f:
        for row in csv.reader(f,delimiter='\t'):
            x,y=left[row[left_ix]],left[row[left_ix+1]]
            z,w=right[row[right_ix]],right[row[right_ix+1]]
            assert x[0]==y[0] and z[0]==w[0]
            chrom=int(x[0].split('CM0213')[1].split('.')[0])-37
            seg.append((z[0],min(z[1],w[1]),max(z[2],w[2]),'chr'+str(chrom)))
            partner[z[0]].add('chr'+str(chrom))
    return seg,partner
def jcvi_evidence(coord,seg,partners):
    chrom,start,end=coord
    contained={r[3] for r in seg if r[0]==chrom and r[1]<=start and end<=r[2]}
    return contained,partners.get(chrom,set())
ep_bed=bed(R/'ep_usda.bed');bk_ep=bed(R/'belonocnema_kinseyi_gcf.bed')
epseg,eppartner=segments(R/'ep_usda.belonocnema_kinseyi_gcf.lifted.anchors.simple',bk_ep,ep_bed,2,0)
epbridge=defaultdict(set)
for r in tsv(R/'helixer_id_map.tsv'):
    epbridge[r['pseudo_chr_id'].removesuffix('.1')].add(r['chromosome_id'])
dk_bed=bed(R/'gca_038277395_1.bed');bk_dk=bed(R/'belonocnema_kinseyi_gcf.bed')
dkseg,dkpartner=segments(R/'DK_JCVI_SIMPLE.tsv',bk_dk,dk_bed,0,2)
# DK source annotation and JCVI annotation use the same native coordinates.
if args.protein_genes:
    if not (args.scan_manifest and args.panel):raise SystemExit('eCAMI placement needs scan manifest and panel')
    scanned={r['species'] for r in tsv(args.scan_manifest)}
    species_genome={r['species'].replace(' ','_'):r['genome'] for r in tsv(args.panel) if r['species']!='Telenomus remus'}
    if set(species_genome)-scanned:raise SystemExit('panel genome lacks pinned scan')
    panel=set(species_genome.values())
    roster=[dict(cynog_id='',species=r['species'],gene_id=r['gene_id'],genome=species_genome[r['species']],protein_id=r['protein_id'])
            for r in tsv(args.protein_genes) if r['species'] in species_genome]
    if len(species_genome) not in {40,48} or not roster:raise SystemExit(f'unexpected eCAMI scan panel: {len(species_genome)} genomes, {len(roster)} proteins')
    bridge={(r['species'],r['gene_id']):r['ragtag_objects'].split(',') if r['ragtag_objects'] else [] for r in tsv(args.bridge or R/'corrected_gene_objects.tsv')}
else:
    roster=list(tsv(R/'membership_table.tsv'))
    bridge={}
needed={r['protein_id'] for r in roster if r['genome']=='GCA_038277395.1'}
old={};new=defaultdict(list)
for path,wanted,out in [(R/'D_kuriphilus_helixer.gff3',needed,old),(R/'GCA_038277395.1_helixer.gff3',None,new)]:
    with path.open() as f:
        for line in f:
            if line.startswith('#'):continue
            a=line.rstrip('\n').split('\t')
            if len(a)<9 or a[2]!='mRNA':continue
            ident=a[8].split('ID=',1)[1].split(';',1)[0]
            if wanted is not None:
                key=ident.removesuffix('.1')
                if key in wanted:out[key]=(a[0],int(a[3]),int(a[4]),a[6])
            else:new[a[0]].append((int(a[3]),int(a[4]),a[6],ident))
for chrom in new:new[chrom].sort()
dkbridge={}
by_chrom=defaultdict(list)
for key,(chrom,start,end,strand) in old.items():by_chrom[chrom].append((start,end,strand,key))
for chrom,queries in by_chrom.items():
    entries=new.get(chrom,[]); j=0; active=[]
    for start,end,strand,key in sorted(queries):
        while j<len(entries) and entries[j][0]<=end:
            active.append(entries[j]);j+=1
        active=[x for x in active if x[1]>=start]
        hits=[(min(end,b)-max(start,a)+1,ident) for a,b,st,ident in active if a<=end and st==strand]
        if hits:
            top=max(score for score,_ in hits)
            winners={ident for score,ident in hits if score==top}
            if len(winners)==1:
                ident=next(iter(winners))
                if ident in dk_bed:dkbridge[key]=dk_bed[ident]
# Corrected curated DK B-lectin unit is an accepted source with its own bridge.
for r in tsv(R/'DK_BLECTIN_JCVI_BRIDGE.tsv'):
    dkbridge[r['biological_unit_id']]=(r['jcvi_scaffold'],int(r['bed_start']),int(r['bed_end']))
# Exact catalog records for every CynOGv3 member.
keys={(r['genome'],r['protein_id']) for r in roster}
catalog={}
for r in tsv(R/'protein_to_chromosome.tsv'):
    key=r['genome'],r['protein_id']
    if key in keys:catalog[key]=r
if args.protein_genes:
    missing=keys-set(catalog)
    expected={(r['genome'],r['protein_id']) for r in roster if (r['species'],r['gene_id']) in bridge}
    if missing!=expected:raise SystemExit(f'eCAMI catalog/bridge mismatch: missing={len(missing)} expected={len(expected)}')
else:
    assert len(catalog)==len(keys),(len(catalog),len(keys))
special={'GCA_010014905.2':'agp_DL','GCA_035041075.1':'agp_AE','GCA_035222095.1':'agp_AQP'}
nine={'AHIE','AMIN','CQUE','DRUB','PCAN','PLAM','SHIS','TPHE','Pe_aceris'}
ragobj={}
for genome in {g for g,p in keys}-{'Belonocnema_kinseyi_GCF','Euceroptres_primus_USDA','GCA_038277395.1'}:
    role=special.get(genome,'agp_'+genome)
    path=R/role
    assert path.is_file(),path
    wanted={catalog[(g,p)]['scaffold'] for g,p in keys if g==genome and (g,p) in catalog}
    objects=defaultdict(set)
    with path.open() as f:
        for line in f:
            if line.startswith('#'):continue
            a=line.rstrip('\n').split('\t')
            if len(a)>5 and a[4]=='W' and a[5] in wanted:objects[a[5]].add(a[0])
    ragobj[genome]=objects
args.out.parent.mkdir(parents=True,exist_ok=True)
counts=Counter();coverage=defaultdict(Counter)
with args.out.open('w') as f:
 writer=csv.DictWriter(f,['cynog_id','species','gene_id','genome','protein_id','route','source_scaffold','state','reason'],delimiter='\t',lineterminator='\n')
 writer.writeheader()
 for r in roster:
    genome,protein=r['genome'],r['protein_id']
    bridged=(r.get('species',''),r.get('gene_id','')) in bridge
    scaffold=','.join(bridge[r['species'],r['gene_id']]) if bridged else catalog[(genome,protein)]['scaffold']
    if genome=='Belonocnema_kinseyi_GCF':
        chrom,reason=placement.place('NCBI-BK',seqid=scaffold);route='NCBI-BK'
    elif genome=='Euceroptres_primus_USDA':
        ids=epbridge.get(protein,set());coord=ep_bed[next(iter(ids))] if len(ids)==1 and next(iter(ids)) in ep_bed else None
        contained,partners=jcvi_evidence(coord,epseg,eppartner) if coord else (set(), eppartner.get(scaffold,set()) if len(eppartner.get(scaffold,set()))==1 else set())
        chrom,reason=placement.place('JCVI',contained_segment_chromosomes=contained,partner_chromosomes=partners);route='JCVI'
        if coord is None and not partners:reason='missing_unique_ep_bridge'
    elif genome=='GCA_038277395.1':
        coord=dkbridge.get(protein)
        contained,partners=jcvi_evidence(coord,dkseg,dkpartner) if coord else (set(),set())
        chrom,reason=placement.place('JCVI',contained_segment_chromosomes=contained,partner_chromosomes=partners);route='JCVI'
        if coord is None:reason='missing_unique_dk_bridge'
    else:
        route='RagTag'
        if bridged:
            objects=bridge[r['species'],r['gene_id']]
            chrom,reason=placement.place('RagTag',ragtag_object=objects[0]) if len(objects)==1 else ('unplaced','ambiguous_annotation_objects' if objects else 'missing_annotation_object')
        elif scaffold.endswith('_RagTag'):obj=scaffold;status='named_object'
        else:
            hits=ragobj[genome].get(scaffold,set())
            obj=next(iter(hits)) if len(hits)==1 else ''
            status='unique_agp_object' if len(hits)==1 else ('ambiguous_agp_objects' if hits else 'no_agp_membership')
        if not bridged:
            chrom,reason=placement.place('RagTag',ragtag_object=obj)
            if not obj:reason=status
    writer.writerow(dict(cynog_id=r['cynog_id'],species=r.get('species',''),gene_id=r.get('gene_id',''),genome=genome,protein_id=protein,route=route,source_scaffold=scaffold,state=chrom if chrom!='unplaced' else '?',reason=reason))
    counts[(genome,reason)]+=1
    coverage[genome]['total']+=1;coverage[genome]['placed']+=int(chrom!='unplaced');coverage[genome][route]+=1
qc={'rows':sum(counts.values()),'placed':sum(v['placed'] for v in coverage.values()),
    'unplaced':sum(v['total']-v['placed'] for v in coverage.values()),
    'reason_counts':{f'{g}|{reason}':n for (g,reason),n in sorted(counts.items())}}
(args.out.parent/'state_qc.json').write_text(json.dumps(qc,indent=2)+'\n')
with (args.out.parent/'placement_coverage.tsv').open('w',newline='') as f:
 w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['genome','total','placed','placed_fraction','route'])
 for genome,c in sorted(coverage.items()):
  route=next(r for r in placement.ROUTES if c[r])
  w.writerow([genome,c['total'],c['placed'],f'{c["placed"]/c["total"]:.6f}',route])
print('states',qc['rows'],'placed',qc['placed'],flush=True)
