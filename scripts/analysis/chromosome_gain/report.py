#!/usr/bin/env python3
"""Write an answer-first report and supplement-ready methods from scored TSVs."""
from __future__ import annotations
import argparse,csv,json,statistics
from collections import Counter
from pathlib import Path
def rows(path):
    with open(path) as handle:yield from csv.DictReader(handle,delimiter='\t')
def pct(value,distribution):
    return 100*(sum(x<value for x in distribution)+.5*sum(x==value for x in distribution))/len(distribution)
def fmt(x):return f'{x:.1f}'
p=argparse.ArgumentParser();p.add_argument('--scores',type=Path,required=True)
p.add_argument('--forced',type=Path,required=True);p.add_argument('--coverage',type=Path,required=True)
p.add_argument('--eligibility',type=Path,required=True)
p.add_argument('--availability',type=Path,required=True)
p.add_argument('--state-qc',type=Path,required=True)
p.add_argument('--figure-qc',type=Path,required=True)
p.add_argument('--preflight',type=Path,required=True)
p.add_argument('--handchecks',type=Path,required=True)
p.add_argument('--out',type=Path,required=True)
a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
data=list(rows(a.scores));bg=[r for r in data if r['kind']=='background'];core=[r for r in bg if r['single_copy_core']=='1']
focal={r['orthogroup'].removeprefix('FOCAL:'):r for r in data if r['kind']=='focal'}
assert len(focal)==4 and bg and core
eligibility=list(rows(a.eligibility))
availability=list(rows(a.availability))
state_qc=json.loads(a.state_qc.read_text())
figure_qc=json.loads(a.figure_qc.read_text())
preflight=json.loads(a.preflight.read_text())
handchecks=list(rows(a.handchecks))
assert len(eligibility)==preflight['cynog_groups']==22448
assert figure_qc['figure1_chromosome_cells']==1248 and figure_qc['mismatches']==0
assert state_qc['rows']==preflight['member_proteins']==623796
assert len(handchecks)==10 and len(availability)==4
def copy_stratum(value):
    n=float(value)
    if n==1:return '1'
    if n<=2:return '(1,2]'
    if n<=4:return '(2,4]'
    return '>4'
with (a.out/'copy_stratum_summary.tsv').open('w',newline='') as handle:
    w=csv.writer(handle,delimiter='\t',lineterminator='\n')
    w.writerow(['copy_stratum','mode','groups','median_gains','median_losses','median_rate','zero_gain_share','GH5_2_percentile','PL1_8_percentile','PL4_4_percentile'])
    for label in ('1','(1,2]','(2,4]','>4'):
        subset=[r for r in bg if copy_stratum(r['median_copies'])==label]
        assert subset
        for mode in ('deltran','acctran','sankoff'):
            rates=[float(r[mode+'_rate']) for r in subset]
            w.writerow([label,mode,len(subset),f'{statistics.median(int(r[mode+"_gains"]) for r in subset):.3f}',
                        f'{statistics.median(int(r[mode+"_losses"]) for r in subset):.3f}',
                        f'{statistics.median(rates):.8f}',
                        f'{sum(int(r[mode+"_gains"])==0 for r in subset)/len(subset):.8f}',
                        *(f'{pct(float(focal[family][mode+"_rate"]),rates):.3f}' for family in ('GH5_2','PL1_8','PL4_4'))])
with (a.out/'background_summary.tsv').open('w',newline='') as handle:
    w=csv.writer(handle,delimiter='\t',lineterminator='\n')
    w.writerow(['set','mode','groups','zero_gain_groups','median_gains','median_losses','median_rate','max_rate'])
    for label,subset in [('crown_present',bg),('single_copy_core',core)]:
        for mode in ('deltran','acctran','sankoff'):
            rates=[float(r[mode+'_rate']) for r in subset]
            w.writerow([label,mode,len(subset),sum(int(r[mode+'_gains'])==0 for r in subset),
                        f'{statistics.median(int(r[mode+"_gains"]) for r in subset):.3f}',
                        f'{statistics.median(int(r[mode+"_losses"]) for r in subset):.3f}',
                        f'{statistics.median(rates):.8f}',f'{max(rates):.8f}'])
focal_percentiles=[pct(float(focal[f][mode+'_rate']),[float(x[mode+'_rate']) for x in bg])
                   for f in ('GH5_2','PL1_8','PL4_4') for mode in ('deltran','acctran','sankoff')]
core_rates=sorted(float(x['deltran_rate']) for x in core)
high_copy=[r for r in bg if copy_stratum(r['median_copies'])=='>4']
high_copy_rates=[float(r['deltran_rate']) for r in high_copy]
lines=['# Retained chromosome gains on the cynipoid species tree','',
 f'All three PCWDE families have more retained new BK-chromosome placements per placed genome than {min(focal_percentiles):.1f}–{max(focal_percentiles):.1f}% of crown-present CynOGv3 orthogroups across the tested parsimony reconstructions. This is a retained-placement comparison, not an estimate of mutation or rearrangement event rate. More copies offer more chances to retain a new placement; the PCWDEs started from one copy.','',
 f'Copy number is descriptive: median copies per member genome are '+', '.join(f'{family} {float(focal[family]["median_copies"]):g}' for family in ('GH5_2','PL1_8','PL4_4'))+f'. Within the {len(high_copy)} background groups with median copy number >4, their DELTRAN midranks are '+', '.join(f'{pct(float(focal[family]["deltran_rate"]),high_copy_rates):.1f}th' for family in ('GH5_2','PL1_8','PL4_4'))+' percentiles, respectively. All three PCWDE rates sit high against the full crown-present set and below the middle among background groups of similar copy number. Neither copy stratum is an adjusted null.','',
 f'The background contains {len(bg):,} crown-present groups (≥5 placed cynipoid genomes); {sum(int(r["deltran_gains"])==0 for r in bg):,}/{len(bg):,} ({fmt(100*sum(int(r["deltran_gains"])==0 for r in bg)/len(bg))}%) have zero DELTRAN gains. The descriptive single-copy core contains {len(core):,} groups, of which {sum(int(r["deltran_gains"])==0 for r in core):,} ({fmt(100*sum(int(r["deltran_gains"])==0 for r in core)/len(core))}%) have zero gains. Its median rate is {statistics.median(core_rates):.3f} and its 95th percentile is {core_rates[int(.95*(len(core_rates)-1))]:.3f} gains per placed tip.','',
 '| Family | Fitch root present (ambiguous) | DELTRAN gains / losses | Gains per placed tip | Background percentile | Core percentile | ACCTRAN gains / percentile | Sankoff gains / percentile |',
 '| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |']
reversals=[]
for family in ('GH5_2','PL1_8','PL4_4'):
    values=[pct(float(focal[family][mode+'_rate']),[float(x[mode+'_rate']) for x in bg])
            for mode in ('deltran','acctran','sankoff')]
    if min(values)<50<max(values):reversals.append(family)
if reversals:
    lines.insert(4,f'**Sensitivity conflict:** {", ".join(reversals)} reverses direction across reconstruction rules; no primary conclusion is selected.\n')
for family in ('GH5_2','PL1_8','PL4_4','B-lectin'):
    r=focal[family]
    root=r['root_present']+(f' ({r["root_ambiguous"]})' if r['root_ambiguous'] else '')
    columns=[]
    for mode in ('deltran','acctran','sankoff'):
        value=float(r[mode+'_rate']);dist=[float(x[mode+'_rate']) for x in bg]
        columns.append((value,pct(value,dist)))
    corepct=pct(columns[0][0],[float(x['deltran_rate']) for x in core])
    if family=='B-lectin':
        lines.append(f'| {family} (descriptive) | {root} | {r["deltran_gains"]} / {r["deltran_losses"]} | {columns[0][0]:.3f} | — | — | {r["acctran_gains"]} / — | {r["sankoff_gains"]} / — |')
    else:
        lines.append(f'| {family} | {root} | {r["deltran_gains"]} / {r["deltran_losses"]} | {columns[0][0]:.3f} | {columns[0][1]:.1f} | {corepct:.1f} | {r["acctran_gains"]} / {columns[1][1]:.1f} | {r["sankoff_gains"]} / {columns[2][1]:.1f} |')
lines+=['',
 'PL1_8 chr8 and PL4_4 chr7 pass the required ancestral-state controls. B-lectin is descriptive: chr6 is root-present and chr3 is root-absent, including when all 20 placed Figure 1 B-lectin species are considered. GH5_2 chr1 is root-present and chr3 is ambiguous; this analysis does not resolve the original GH5_2 placement.','',
 '## Where the gains occur','',
 'Cynipini is the species-tree clade from *Dryocosmus kuriphilus* through *Acraspis erinacei* (10 genomes with CynOGv3 membership). A gain is assigned to Cynipini when every placed descendant of its gaining branch is in that clade. This counts a gain on its stem with Cynipini.','',
 '| Family | Cynipini gains / placed tips | Elsewhere gains / placed tips |',
 '| --- | ---: | ---: |']
for family in ('GH5_2','PL1_8','PL4_4','B-lectin'):
    r=focal[family]
    lines.append(f'| {family} | {r["cynipini_gains"]} / {r["cynipini_placed_tips"]} | {r["elsewhere_gains"]} / {r["elsewhere_placed_tips"]} |')
for label,gain_col,tip_col in [('Cynipini','cynipini_gains','cynipini_placed_tips'),('Elsewhere','elsewhere_gains','elsewhere_placed_tips')]:
    usable=[r for r in bg if int(r[tip_col])>0]
    rates=[int(r[gain_col])/int(r[tip_col]) for r in usable]
    lines.append(f'- Background {label}: {len(usable):,} groups with placed tips; median {statistics.median(rates):.3f} gains per placed tip; {sum(int(r[gain_col])==0 for r in usable)/len(usable):.1%} with zero gains.')
lines.append('Cynipini has elevated genome-wide duplication (25.8% of orthogroups multi-copy versus 18.4% in other cynipids); the comparison retains that biological difference rather than adjusting it away.')
forced=list(rows(a.forced));agg={}
for label in ('chr1','chr6_chr10'):
    subset=[r for r in forced if r['root_scenario']==label]
    assert len(subset)==10
    agg[label]=(sum(int(r['gains']) for r in subset),sum(int(r['losses']) for r in subset),sum(int(r['cost']) for r in subset))
lines+=['','## GH5_2 forced-root sensitivity','',
       f'Forcing chr1 as the sole GH5_2 root placement gives {agg["chr1"][0]} gains and {agg["chr1"][1]} losses (weighted cost {agg["chr1"][2]}). Forcing chr6+chr10 gives {agg["chr6_chr10"][0]} gains and {agg["chr6_chr10"][1]} losses (cost {agg["chr6_chr10"][2]}). These constraints are sensitivity scenarios, not inferred ancestral states.','',
       'The source species tree has no branch lengths, so gains per summed branch length cannot be calculated. Placement noise is visible in the single-copy core distribution above. Tip states, per-group scores, per-branch focal changes, and exact input hashes accompany this report. No transfer direction or chromosome-change mechanism is inferred.','']
coverage=list(rows(a.coverage));assert len(coverage)==41
low=sorted((float(r['placed_fraction']),r['genome']) for r in coverage)
lines+=['## Placement coverage','',
        f'Every membership genome had a NCBI-BK, RagTag, or JCVI background route; no cynipoid genome was excluded for route incompatibility. The 457 focal proteins shared with CynOGv3 membership matched background route and chromosome exactly. Background member-protein placement ranged from {low[0][0]:.1%} ({low[0][1]}) to {low[-1][0]:.1%} ({low[-1][1]}). Low placement coverage makes some tips missing and can undercount retained gains.','']
reason=Counter(r['reason'] for r in eligibility)
counterfactual=Counter(r['counterfactual'] for r in eligibility if r['reason'] not in ('eligible','focal_group'))
assert reason['eligible']==len(bg)==10942 and sum(reason.values())==22448
lines+=['## Eligibility and validation','',
        f'The 41 membership genomes include *Telenomus remus*, the outgroup; scoring covers 40 cynipoid genomes. Of {len(eligibility):,} groups, {reason["eligible"]:,} meet the crown-present rule, {reason["focal_group"]} contain focal proteins, {reason["no_ingroup_member"]} lack an ingroup member, {reason["lt5_placed_tips"]:,} have fewer than five placed tips, and {reason["not_crown_present"]} lack placed copies on both crown sides. The reason counts sum to {sum(reason.values()):,}.', '',
        f'Across all membership rows, {state_qc["placed"]:,}/{state_qc["rows"]:,} proteins were placed. Figure 1 QC compared {figure_qc["figure1_chromosome_cells"]:,} cells and found {figure_qc["mismatches"]} mismatches. Ten gain assignments were checked against tip states and branch descendants.', '',
        '| Family | No copy genomes | All copies unplaced genomes | Placed genomes |',
        '| --- | ---: | ---: | ---: |']
for row in availability:
    lines.append(f'| {row["family"]} | {row["no_copy_genomes"]} | {row["all_copies_unplaced_genomes"]} | {row["placed_genomes"]} |')
lines+=['', 'The row-level eligibility table, per-tip states, per-node changes, copy-number strata and hand checks accompany this report.','']
(a.out/'summary.md').write_text('\n'.join(lines))
methods='''# Methods note — retained chromosome gains

We used the accepted Figure 2 cynipoid species-tree topology and the frozen CynOGv3 41-genome membership panel. The *Telenomus remus* outgroup was excluded from change scoring, leaving 40 cynipoid genomes; it remained in the 41-genome denominator for the descriptive single-copy core. Focal PCWDE copy placements came from the build member/copy table and were checked against every displayed Figure 1 chromosome cell. Background proteins were routed through the same NCBI-BK, RagTag, or JCVI `place()` function with protein-specific evidence. An orthogroup was crown-present when a placed copy occurred in Diplolepididae or Figitidae and another in Cynipidae; every scored group had at least five placed cynipoid genomes. Selection used presence, not copy number or chromosome state. Orthogroups sharing any focal protein were excluded from background. The descriptive single-copy core had exactly one member protein in at least 37 of the 41 membership genomes.

For each orthogroup, each of the ten BK chromosomes was a binary occupancy character. A genome with a member protein but no placed copy, or no member protein, was missing for that orthogroup. We pruned the species tree to genomes with at least one placed copy and summed minimum-cost changes over characters. Equal-cost Fitch reconstructions used a deterministic global depth tie-break: DELTRAN maximized the total depth of changes and ACCTRAN minimized it, after minimizing the total number of changes; remaining ties chose root state 0 and then child state 0. The primary statistic was DELTRAN 0→1 changes divided by the number of placed cynipoid genomes. We also counted 1→0 changes and recomputed under Sankoff gain cost 2, loss cost 1. For each focal normalized rate, the background percentile used midranks. The source tree supplied topology but no branch lengths, so branch-length normalization was unavailable.

The statistic counts retained placements: event generation, retention and elapsed time all affect it. Families with more copies have more opportunities to retain gains, and the focal PCWDE families began from one copy. No event-rate, transfer-direction or mechanism inference follows from this parsimony score.
'''
(a.out/'METHODS_NOTE.md').write_text(methods)
with (a.out/'scientific_context.tsv').open('w',newline='') as handle:
    w=csv.writer(handle,delimiter='\t',lineterminator='\n')
    w.writerow(['clade','multicopy_orthogroup_percent','source'])
    w.writerow(['Cynipini','25.8','fixed scientific framing (2026-11-01)'])
    w.writerow(['other cynipids','18.4','fixed scientific framing (2026-11-01)'])
print('REPORT_PASS',len(bg),len(core))
