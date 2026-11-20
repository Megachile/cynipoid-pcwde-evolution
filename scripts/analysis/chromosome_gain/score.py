#!/usr/bin/env python3
"""Score retained BK-chromosome occupancy gains on the paper species tree."""
from __future__ import annotations
import argparse, ast, csv, hashlib, json, statistics
from collections import Counter, defaultdict
from pathlib import Path
from tree import parse_newick, prune, fitch

CHROMS=tuple(f'chr{i}' for i in range(1,11))
FOCAL=('GH5_2','PL1_8','PL4_4','B-lectin')

def rows(path):
    with open(path,newline='') as handle:yield from csv.DictReader(handle,delimiter='\t')

def write(path,columns,records):
    with open(path,'w',newline='') as handle:
        w=csv.DictWriter(handle,columns,delimiter='\t',lineterminator='\n',extrasaction='ignore')
        w.writeheader();w.writerows(records)

def descendants(node):
    if not node.children:return frozenset([node.name])
    return frozenset().union(*(descendants(x) for x in node.children))

def node_id(node):
    if not node.children:return node.name
    members=','.join(sorted(descendants(node))).encode()
    return 'node_'+hashlib.sha1(members).hexdigest()[:12]

def anchor_clade(root,names):
    eligible=[]
    def visit(node):
        members=descendants(node)
        if names<=members:eligible.append((len(members),members))
        for child in node.children:visit(child)
    visit(root)
    return min(eligible,key=lambda x:x[0])[1]

def score_character(tree,ones,mode,root_override=None):
    """Sankoff DP; equal-cost Fitch with fixed ACCTRAN/DELTRAN depth tie-break."""
    gain_cost,loss_cost=(2,1) if mode=='sankoff' else (1,1)
    memo={}
    def visit(node,depth):
        key=id(node)
        if not node.children:
            actual=1 if node.name in ones else 0
            memo[key]={actual:(0,0,[])}
            return
        for child in node.children:visit(child,depth+1)
        choices={}
        for state in (0,1):
            cost=0;weighted=0;assign=[]
            for child in node.children:
                options=[]
                for child_state,(child_cost,child_weight,child_assign) in memo[id(child)].items():
                    change=child_state!=state
                    extra=(gain_cost if state==0 else loss_cost) if change else 0
                    depth_term=(depth+1) if change else 0
                    if mode=='deltran':depth_term=-depth_term
                    if mode=='sankoff':depth_term=0
                    options.append((child_cost+extra,child_weight+depth_term,child_state,child_assign))
                best=min(options,key=lambda x:(x[0],x[1],x[2]))
                cost+=best[0];weighted+=best[1]
                assign.append((child,best[2],best[3]))
            choices[state]=(cost,weighted,assign)
        memo[key]=choices
    visit(tree,0)
    options=memo[id(tree)]
    if root_override is not None:
        state=root_override
    else:
        state=min(options,key=lambda s:(options[s][0],options[s][1],s))
    gains=[];losses=[]
    def walk(node,state,depth):
        for child,child_state,unused in memo[id(node)][state][2]:
            if child_state!=state:
                event=dict(parent=node_id(node),child=node_id(child),depth=depth+1,
                           descendant_tips=','.join(sorted(descendants(child))))
                (gains if state==0 else losses).append(event)
            walk(child,child_state,depth+1)
    walk(tree,state,0)
    return state,options[state][0],gains,losses

def taxon_groups(source):
    module=ast.parse(source.read_text())
    for node in module.body:
        if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='TAXON_LABEL_GROUPS' for t in node.targets):
            return {name:set(members) for name,members in ast.literal_eval(node.value)}
    raise ValueError('Figure 2 taxon groups not found')

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--states',type=Path,required=True)
    p.add_argument('--focal',type=Path,required=True)
    p.add_argument('--panel',type=Path,required=True)
    p.add_argument('--tree',type=Path,required=True)
    p.add_argument('--figure-code',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    panel={r['genome']:r['species'].replace(' ','_') for r in rows(a.panel)}
    groups=taxon_groups(a.figure_code)
    assert set(panel.values())<=set().union(*groups.values(),{'Telenomus_remus'})
    ingroup={g:s for g,s in panel.items() if s!='Telenomus_remus'}
    assert len(panel)==41 and len(ingroup)==40
    tree=parse_newick(a.tree.read_text())
    assert len(descendants(tree))==49
    tree=prune(tree,set(ingroup.values()))
    cynipini=anchor_clade(tree,{'Dryocosmus_kuriphilus','Acraspis_erinacei'})
    assert len(cynipini)==10,cynipini  # Eleven in the full tree; Cynips lacks CynOGv3 membership.
    membership=defaultdict(lambda:defaultdict(lambda:dict(placed=set(),copies=0)))
    all_copy_counts=defaultdict(Counter)
    all_group_ids=set()
    placed_proteins=0
    occupied_cells=set()
    focal_proteins={}
    for r in rows(a.focal):
        for protein in r['protein_ids'].split(','):
            focal_proteins[(r['genome'],protein)]=r['family']
    focal_groups=set(); parity=Counter(); mismatches=[]
    focal_placements={(g,p):(r['bk_chromosome'],r['placement_route'])
        for r in rows(a.focal) for g in [r['genome']] for p in r['protein_ids'].split(',')}
    for r in rows(a.states):
        g=r['genome'];og=r['cynog_id'];key=g,r['protein_id']
        all_group_ids.add(og)
        if r['state'] in CHROMS:
            placed_proteins+=1
            occupied_cells.add((og,g))
        if key in focal_proteins:
            focal_groups.add(og)
            expect=focal_placements[key]
            got=('unplaced' if r['state']=='?' else r['state'],r['route'])
            parity['shared']+=1
            if expect!=got:mismatches.append((key,expect,got))
        all_copy_counts[og][g]+=1
        if g not in ingroup:continue
        item=membership[og][g];item['copies']+=1
        if r['state'] in CHROMS:item['placed'].add(r['state'])
    assert not mismatches,mismatches[:5]
    assert parity['shared']==457,parity
    focal_tips=defaultdict(lambda:defaultdict(lambda:dict(placed=set(),copies=0)))
    for r in rows(a.focal):
        g=r['genome']
        if g not in ingroup:continue
        item=focal_tips[r['family']][g]
        item['copies']+=len(r['protein_ids'].split(','))
        if r['bk_chromosome'] in CHROMS:item['placed'].add(r['bk_chromosome'])
    assert len(all_group_ids)==22448 and placed_proteins==454288
    eligibility=[]
    reason_counts=Counter()
    cross_counts=Counter()
    for og in sorted(all_group_ids):
        tips=membership.get(og,{})
        placed={g:v for g,v in tips.items() if v['placed']}
        member_sources={next(k for k,v in groups.items() if s in v) for g in tips for s in [ingroup[g]]}
        placed_sources={next(k for k,v in groups.items() if s in v) for g in placed for s in [ingroup[g]]}
        member_crown='Cynipidae' in member_sources and bool({'Figitidae','Diplolepididae'} & member_sources)
        placed_crown='Cynipidae' in placed_sources and bool({'Figitidae','Diplolepididae'} & placed_sources)
        if og in focal_groups:reason='focal_group'
        elif not tips:reason='no_ingroup_member'
        elif len(placed)<5:reason='lt5_placed_tips'
        elif not placed_crown:reason='not_crown_present'
        elif prune(tree,{ingroup[g] for g in placed}) is None:reason='prune_none'
        else:reason='eligible'
        counterfactual='meets_both' if len(tips)>=5 and member_crown else ('absent_by_membership_crown_scope' if not member_crown else 'lt5_member_tips')
        if reason not in ('eligible','focal_group') and counterfactual=='meets_both':counterfactual='excluded_only_unplaced'
        eligibility.append(dict(orthogroup=og,reason=reason,member_genomes=len(tips),placed_tips=len(placed),
                                member_side_sets=','.join(sorted(member_sources)),placed_side_sets=','.join(sorted(placed_sources)),
                                member_ge5=int(len(tips)>=5),member_crown_present=int(member_crown),counterfactual=counterfactual))
        reason_counts[reason]+=1;cross_counts[(reason,counterfactual)]+=1
    assert reason_counts['eligible']==10942 and reason_counts['focal_group']==11
    write(a.out/'background_eligibility.tsv',['orthogroup','reason','member_genomes','placed_tips','member_side_sets','placed_side_sets',
          'member_ge5','member_crown_present','counterfactual'],eligibility)
    records=[]; changes=[]; tiprows=[];forced_rows=[]
    selected={}
    for og,tips in membership.items():
        if og in focal_groups:continue
        placed={g:v for g,v in tips.items() if v['placed']}
        if len(placed)<5:continue
        sources={next(k for k,v in groups.items() if s in v) for g in placed for s in [ingroup[g]]}
        if 'Cynipidae' not in sources or not ({'Figitidae','Diplolepididae'}&sources):continue
        selected[og]=(tips,placed)
    for fam,tips in focal_tips.items():
        selected['FOCAL:'+fam]=(tips,{g:v for g,v in tips.items() if v['placed']})
    for ident,(tips,placed) in sorted(selected.items()):
        is_focal=ident.startswith('FOCAL:')
        if len(placed)<5:continue
        names={ingroup[g] for g in placed};sub=prune(tree,names)
        if sub is None:continue
        for g,v in tips.items():
            tiprows.append(dict(orthogroup=ident,genome=g,species=ingroup[g],
                                chromosomes=','.join(sorted(v['placed'])),placed=int(bool(v['placed'])),copies=v['copies']))
        root_sets={};totals={mode:Counter() for mode in ('deltran','acctran','sankoff')}
        for ch in CHROMS:
            ones={ingroup[g] for g,v in placed.items() if ch in v['placed']}
            possible,fitch_cost=fitch(sub,ones)
            root_sets[ch]=''.join(map(str,sorted(possible[id(sub)])))
            for mode in totals:
                root_state,cost,gains,losses=score_character(sub,ones,mode)
                if mode!='sankoff' and cost!=fitch_cost:
                    raise AssertionError((ident,ch,mode,cost,fitch_cost))
                totals[mode]['gains']+=len(gains);totals[mode]['losses']+=len(losses);totals[mode]['cost']+=cost
                if mode=='deltran':
                    for event in gains:
                        members=set(event['descendant_tips'].split(','))
                        totals[mode]['cynipini_gains' if members<=cynipini else 'elsewhere_gains']+=1
                if is_focal:
                    for kind,events in (('gain',gains),('loss',losses)):
                        for event in events:
                            changes.append(dict(orthogroup=ident,mode=mode,chromosome=ch,change=kind,
                                                root_state=root_state,**event))
            if ident=='FOCAL:GH5_2':
                for label,root_chroms in [('chr1',{'chr1'}),('chr6_chr10',{'chr6','chr10'})]:
                    root_state,cost,gains,losses=score_character(sub,ones,'deltran',int(ch in root_chroms))
                    forced_rows.append(dict(root_scenario=label,chromosome=ch,cost=cost,
                                            gains=len(gains),losses=len(losses)))
                    for kind,events in (('gain',gains),('loss',losses)):
                        for event in events:
                            changes.append(dict(orthogroup=ident,mode='forced_'+label,chromosome=ch,
                                                change=kind,root_state=root_state,**event))
        copy_counts=[v['copies'] for v in tips.values()]
        core=(not is_focal and sum(n==1 for n in all_copy_counts[ident].values())>=37)
        rec=dict(orthogroup=ident,kind='focal' if is_focal else 'background',
                 placed_tips=len(placed),member_genomes=len(tips),median_copies=statistics.median(copy_counts),
                 single_copy_core=int(core),root_present=','.join(ch for ch,s in root_sets.items() if s=='1'),
                 root_ambiguous=','.join(ch for ch,s in root_sets.items() if s=='01'),
                 cynipini_placed_tips=sum(ingroup[g] in cynipini for g in placed),
                 elsewhere_placed_tips=sum(ingroup[g] not in cynipini for g in placed),
                 cynipini_gains=totals['deltran']['cynipini_gains'],
                 elsewhere_gains=totals['deltran']['elsewhere_gains'])
        for mode in totals:
            rec[mode+'_gains']=totals[mode]['gains'];rec[mode+'_losses']=totals[mode]['losses']
            rec[mode+'_rate']=round(totals[mode]['gains']/len(placed),8)
        records.append(rec)
    columns=['orthogroup','kind','placed_tips','member_genomes','median_copies','single_copy_core','root_present','root_ambiguous',
             'cynipini_placed_tips','elsewhere_placed_tips','cynipini_gains','elsewhere_gains']
    for mode in ('deltran','acctran','sankoff'):columns += [mode+'_gains',mode+'_losses',mode+'_rate']
    focal_check={r['orthogroup']:r for r in records if r['kind']=='focal'}
    for family,chrom in [('PL1_8','chr8'),('PL4_4','chr7')]:
        assert chrom in focal_check['FOCAL:'+family]['root_present'].split(','),f'{family} root control failed'
    write(a.out/'orthogroup_scores.tsv',columns,records)
    write(a.out/'tip_states.tsv',['orthogroup','genome','species','chromosomes','placed','copies'],tiprows)
    write(a.out/'focal_node_changes.tsv',['orthogroup','mode','chromosome','change','root_state','parent','child','depth','descendant_tips'],changes)
    write(a.out/'gh5_forced_roots.tsv',['root_scenario','chromosome','cost','gains','losses'],forced_rows)
    root_diag=[]
    for r in records:
        if r['kind']!='focal':continue
        for ch in CHROMS:
            state='1' if ch in r['root_present'].split(',') else ('01' if ch in r['root_ambiguous'].split(',') else '0')
            root_diag.append(dict(panel='41_member_genomes',family=r['orthogroup'].removeprefix('FOCAL:'),
                                  chromosome=ch,root_state_set=state,placed_tips=r['placed_tips']))
    all_blectin=defaultdict(set)
    rename={'Ceroptres_kovalevi':'Ceroptres_masudai','Eschatocerus_acaciae':'Eschatocerus_acacia',
            'Eubothrus_centaureae':'Isocolus_centaureae'}
    for r in rows(a.focal):
        if r['family']=='B-lectin' and r['bk_chromosome'] in CHROMS:
            species=r['species'].replace(' ','_');species=rename.get(species,species)
            all_blectin[species].add(r['bk_chromosome'])
    assert len(all_blectin)==20
    full_blectin_tree=prune(parse_newick(a.tree.read_text()),set(all_blectin))
    for ch in CHROMS:
        sets,_=fitch(full_blectin_tree,{s for s,labels in all_blectin.items() if ch in labels})
        root_diag.append(dict(panel='all_figure1_placed',family='B-lectin',chromosome=ch,
                              root_state_set=''.join(map(str,sorted(sets[id(full_blectin_tree)]))),placed_tips=20))
    write(a.out/'root_diagnostic.tsv',['panel','family','chromosome','root_state_set','placed_tips'],root_diag)
    focal_availability=[]
    for fam in FOCAL:
        tips=focal_tips[fam]
        absent=sum(g not in tips for g in ingroup)
        unplaced=sum(g in tips and not tips[g]['placed'] for g in ingroup)
        focal_availability.append(dict(family=fam,no_copy_genomes=absent,all_copies_unplaced_genomes=unplaced,
                                        placed_genomes=sum(bool(v['placed']) for v in tips.values())))
    write(a.out/'focal_availability.tsv',['family','no_copy_genomes','all_copies_unplaced_genomes','placed_genomes'],focal_availability)
    summary={'genomes_with_membership':len(panel),'cynipoid_genomes':len(ingroup),
             'eligible_background_groups':sum(r['kind']=='background' for r in records),
             'single_copy_core_groups':sum(r['single_copy_core'] for r in records),
             'focal_og_exclusions':sorted(focal_groups),'focal_shared_protein_route_parity':parity['shared'],
             'cynipini_species':sorted(cynipini),
             'total_member_proteins':623796,'placed_member_proteins':placed_proteins,
             'occupied_genome_group_cells':len(occupied_cells),
             'eligibility_reason_counts':dict(sorted(reason_counts.items())),
             'eligibility_reason_counterfactual':[{ 'reason':k[0], 'counterfactual':k[1], 'groups':v}
                 for k,v in sorted(cross_counts.items())],
             'focal_availability':focal_availability}
    (a.out/'score_qc.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,sort_keys=True))

if __name__=='__main__':main()
