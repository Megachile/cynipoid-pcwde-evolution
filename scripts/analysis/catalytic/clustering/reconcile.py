import argparse, collections, csv, json, pathlib, random, re
parser=argparse.ArgumentParser()
parser.add_argument("--calls", type=pathlib.Path, required=True)
parser.add_argument("--authority-root", type=pathlib.Path, required=True)
parser.add_argument("--out", type=pathlib.Path, required=True)
parser.add_argument("--crosswalk", type=pathlib.Path, required=True)
args=parser.parse_args()
ROOT=args.authority_root
rows=list(csv.DictReader(args.calls.open(),delimiter="\t"))

sites={family:{site:(sum(r['conserved']=='True' for r in rows if r['alignment_set']=='de_novo' and r['family']==family and r['site_label']==site and r['site_status']=='observed'),sum(1 for r in rows if r['alignment_set']=='de_novo' and r['family']==family and r['site_label']==site and r['site_status']=='observed')) for site in labels} for family,labels in {'GH5_2':['acid_base','nucleophile'],'PL1_8':['catalytic_arginine'],'PL4_4':['proton_abstractor','proton_donor']}.items()}

def read_tree(path):
 s=path.read_text().strip(); toks=re.findall(r'[^(),:;\s]+|[(),:;]',s); nodes=[]; stack=[]; children=[]; names=[]; i=0
 def new(kids,name):
  n=len(children);children.append(kids);names.append(name);return n
 def parse():
  nonlocal i
  if toks[i]=='(':
   i+=1;kids=[parse()]
   while toks[i]==',':i+=1;kids.append(parse())
   assert toks[i]==')';i+=1
   label=None
   if toks[i] not in (':',',',')',';'):label=toks[i];i+=1
   n=new(kids,label)
  else:n=new([],toks[i]);i+=1
  if toks[i]==':':i+=2
  return n
 root=parse();return children,names,root

N=10000
seed=20261014
metadata={(r['family'],r['protein_id']):r for r in csv.DictReader(open(ROOT/'figure_s1/metadata.tsv'),delimiter='\t')}
crosswalk={(r['family'],r['source_id']):r for r in csv.DictReader(open(args.crosswalk),delimiter='\t')}
coverage=[]; results={}; summary=collections.Counter()
for family, accepted_sites in sites.items():
    children,names,root=read_tree(ROOT/f'figure3_s1_origin_current/{family}_tree.treefile')
    raw_tips={name:i for i,name in enumerate(names) if not children[i]}
    tips={sid:raw_tips[row['tree_id']] for (fam,sid),row in crosswalk.items() if fam==family}
    assert len(tips)==len(raw_tips) and len(set(tips.values()))==len(raw_tips)
    site_rows=[r for r in rows if r['alignment_set']=='de_novo' and r['family']==family]
    by_id=collections.defaultdict(dict)
    for r in site_rows:
        sid=r['sequence_id']; site=r['site_label']
        assert site not in by_id[sid],(family,sid,site)
        by_id[sid][site]=r
        match=sid in tips
        status=r['site_status']
        reason=('mapped_observed' if status=='observed' else 'mapped_missing_site') if match else ('accepted_observation_absent_tree' if status=='observed' else 'accepted_missing_site_absent_tree')
        coverage.append(dict(alignment_set='de_novo',family=family,source_id=sid,normalized_id=sid,tree_match=int(match),site=site,observed_state=r['observed_aa'],site_status=status,eligible=int(match and status=='observed'),reason=reason))
        summary[(family,reason)]+=1
    for r in rows:
        if r['family']!=family or r['alignment_set']=='de_novo':continue
        sid=r['sequence_id'];match=sid in tips;status=r['site_status']
        reason=('mapped_observed' if status=='observed' else 'mapped_missing_site') if match else ('accepted_observation_absent_tree' if status=='observed' else 'accepted_missing_site_absent_tree')
        coverage.append(dict(alignment_set=r['alignment_set'],family=family,source_id=sid,normalized_id=sid,tree_match=int(match),site=r['site_label'],observed_state=r['observed_aa'],site_status=status,eligible=int(match and status=='observed'),reason=reason))
        summary[(family,reason)]+=1
    for sid in sorted(tips.keys()-by_id.keys()):
        coverage.append(dict(alignment_set='tree_only',family=family,source_id='',normalized_id=sid,tree_match=1,site='',observed_state='',site_status='',eligible=0,reason='tree_tip_without_accepted_call'))
        summary[(family,'tree_tip_without_accepted_call')]+=1
    def score(states):
        state=[3]*len(children); changes=0
        for sid,v in states.items(): state[tips[sid]]=2 if v else 1
        for i,kids in enumerate(children):
            if kids:
                v=state[kids[0]]
                for child in kids[1:]:
                    both=v&state[child]
                    if both:v=both
                    else:v|=state[child]; changes+=1
                state[i]=v
        return changes
    fam_result={'tree_tips':len(tips),'mapped_accepted_proteins':len(tips&by_id.keys()),'accepted_proteins_absent_tree':sorted(by_id.keys()-tips.keys()),'sites':{}}
    for site,expected in accepted_sites.items():
        all_observed={sid:r[site] for sid,r in by_id.items() if site in r and r[site]['site_status']=='observed'}
        got=(sum(r['conserved']=='True' for r in all_observed.values()),len(all_observed))
        assert got==expected,(family,site,got,expected)
        observed={sid:r['conserved']=='False' for sid,r in all_observed.items() if sid in tips}
        ids=sorted(observed); k=sum(observed.values()); actual=score(observed)
        rng=random.Random(seed+sum(map(ord,family+site)))
        simulations=[]
        if k:
            for _ in range(N):
                selected=set(rng.sample(ids,k));simulations.append(score({sid:sid in selected for sid in ids}))
        fam_result['sites'][site]=dict(published_retained=got[0],published_total=got[1],mapped_retained=len(ids)-k,mapped_total=len(ids),absent_observed=sorted(all_observed.keys()-tips.keys()),exceptions=k,observed_changes=actual,permutations=len(simulations),seed=seed+sum(map(ord,family+site)),null_mean=sum(simulations)/N if simulations else None,p_clustered=(1+sum(v<=actual for v in simulations))/(N+1) if simulations else None,null_histogram=dict(sorted(collections.Counter(simulations).items())))
    # Preserve each observed site state and its unknown mask as one vector.
    vectors=collections.defaultdict(list)
    for sid in sorted(tips&by_id.keys()):
        vector=tuple(None if s not in by_id[sid] or by_id[sid][s]['site_status']!='observed' else by_id[sid][s]['conserved']=='False' for s in accepted_sites)
        if any(v is not None for v in vector): vectors[tuple(v is not None for v in vector)].append((sid,vector))
    actual_states={sid:any(v is True for v in vector) for items in vectors.values() for sid,vector in items}
    actual=score(actual_states);k=sum(actual_states.values())
    rng=random.Random(seed+sum(map(ord,family)))
    simulations=[]
    if k:
        for _ in range(N):
            permuted={}
            for items in vectors.values():
                states=[v for _,v in items];rng.shuffle(states)
                permuted.update({sid:any(x is True for x in vector) for (sid,_),vector in zip(items,states)})
            simulations.append(score(permuted))
    fam_result['any_exception']=dict(population=len(actual_states),exceptions=k,observed_changes=actual,mask_strata={''.join('1' if x else '0' for x in mask):len(items) for mask,items in vectors.items()},permutations=len(simulations),seed=seed+sum(map(ord,family)),null_mean=sum(simulations)/N if simulations else None,p_clustered=(1+sum(v<=actual for v in simulations))/(N+1) if simulations else None,null_histogram=dict(sorted(collections.Counter(simulations).items())))
    results[family]=fam_result
args.out.mkdir(parents=True,exist_ok=True)
with (args.out/'CALL_COVERAGE.tsv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=list(coverage[0]),delimiter='\t');w.writeheader();w.writerows(coverage)
(args.out/'ANALYSIS.json').write_text(json.dumps({'accepted_call_source':str(args.calls),'tree_source':str(ROOT/'figure3_s1_origin_current'),'permutation_count':N,'families':results,'coverage_reason_counts':{f'{f}/{r}':n for (f,r),n in sorted(summary.items())}},indent=2)+'\n')
print(json.dumps({f:{'sites':{s:(x['mapped_total'],x['exceptions'],x['observed_changes'],x['p_clustered']) for s,x in v['sites'].items()},'any_exception':v['any_exception']} for f,v in results.items()},indent=2))
