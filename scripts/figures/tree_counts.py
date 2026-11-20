#!/usr/bin/env python3
"""Exploratory combined Figure 1/2 view; source loaders remain authoritative."""
import csv
import os
import math
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from Bio import Phylo
import presence_counts as f1
import chromosome_counts as f2
from tree_counts_tree import graft, verify_graft

ROOT = Path(__file__).resolve().parents[2]
ECOLOGY = ROOT / 'config/tree_counts_outgroup_ecology.tsv'
REMNANT = ROOT / 'config/tree_counts_remnant.tsv'
FAMILIES = ('PL1_8', 'PL4_4', 'GH5_2', 'CBM13')
MAIN_LABELS = '''Drosophila melanogaster|Contarinia nasturtii|Mayetiola destructor|Athalia rosae|Microplitis demolitor|Apis mellifera|Habropoda laboriosa|Ceratina calcarata|Nasonia vitripennis|Telenomus remus|Pediaspis aceris|Diplolepis rosae|Diplolepis spinosa|Parnips nigripes|Euceroptres primus|Alloxysta arcuata|Ganaspis brasiliensis|Leptopilina boulardi|Eschatocerus acaciae|Timaspis phoenixopodos|Fumariphilus hypecoi|Eubothrus centaureae|Aulacidea tavakolii|Qwaqwaia scolopiae|Synophrus hispanicus|Synergus laeviventris|Synergus itoensis|Synergus japonicus|Periclistus caninae|Diastrophus rubi|Aylax minor|Ceroptres kovalevi|Protobalandricus spectabilis|Dryocosmus kuriphilus|Melikaiella ostensackeni|Amphibolips quercusostensackenii|Andricus quercusramuli|Cynips quercusfolii|Belonocnema kinseyi|Druon quercuslanigerum'''.split('|')
SCAN = f1.INPUT_DIR

def load_extended_gh1():
    """Use the GH1 counts produced from this build's full competitive scans."""
    path = Path(os.environ['PCWDE_GH1_DERIVED_COUNTS'])
    with path.open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    if len(rows) != 64 or len({r['species'] for r in rows}) != 64:
        raise ValueError('Expected 64 unique derived GH1 counts')
    counts = {r['species']: int(r['GH1_e72']) for r in rows}
    if any(value < 0 for value in counts.values()):
        raise ValueError('Negative derived GH1 count')
    sources = {sp: 'cazy_counts/derived_gh1e72_gene_counts.tsv' for sp in counts}
    return counts, sources

def bracket_spans(tips, inner):
    """Return exact contiguous clade spans, broadest lane first."""
    diptera = {'Drosophila_melanogaster','Bradysia_coprophila','Contarinia_nasturtii','Mayetiola_destructor'}
    sawflies = {'Athalia_rosae','Neodiprion_lecontei','Cephus_cinctus','Orussus_abietinus'}
    groups = [
        ('Diptera', diptera, 0),
        ('Hymenoptera', set(tips)-diptera, 0),
        ('Apocrita', set(tips)-diptera-sawflies, 1),
        ('Cynipoidea', set(inner), 2),
    ] + [(label, set(members), 3) for label,members in f2.TAXON_LABEL_GROUPS]
    spans=[]
    for label, members, lane in groups:
        members &= set(tips)
        if not members: continue
        rows=[i for i,sp in enumerate(tips) if sp in members]
        if set(tips[min(rows):max(rows)+1]) != members:
            raise ValueError(f'{label}: members are not a contiguous tip span')
        spans.append((label,min(rows),max(rows),lane))
    return spans

def cell_shade(value, cap, family):
    """Use the installed Figure 2 count gradient."""
    return f2.lerp_color(min(value / cap, 1.0) if cap else 0, family)

def numeral_color(value, cap):
    return 'white' if value / cap > .4 else '#333'

def percentile95(values):
    values=sorted(v for v in values if isinstance(v,int) and v>0)
    if not values: return 1
    rank=.95*(len(values)-1); lo=math.floor(rank); hi=math.ceil(rank)
    return max(1,values[lo]+(rank-lo)*(values[hi]-values[lo]))

def build(outdir, roster='full', layout='portrait'):
    outdir.mkdir(parents=True, exist_ok=True)
    t1 = Phylo.read(str(f1.TREE_FILE), 'newick')
    t2 = f2.build_display_tree(f2.DESIRED_ORDER)
    tree = graft(t1, t2, f2.PRUNE_SPECIES)
    verify_graft(t1, t2, tree, f2.PRUNE_SPECIES)
    f2.reorder_tree(tree, [t.name for t in t1.get_terminals() if t.name not in {x.name for x in t1.common_ancestor('Diplolepis_rosae','Belonocnema_kinseyi').get_terminals()}] + f2.DESIRED_ORDER[1:])
    # Reorder only the root children after the inherited within-clade ordering.
    if 'Drosophila_melanogaster' not in {t.name for t in tree.root.clades[0].get_terminals()}:
        raise ValueError('Expected Diptera in first root child')
    tree.root.clades.reverse()
    tips = [x.name for x in tree.get_terminals()]
    if roster == 'main':
        labels = {f2.short_label(sp, {}): sp for sp in tips}
        if any(label not in labels for label in MAIN_LABELS):
            raise ValueError(f'Main labels absent: {set(MAIN_LABELS)-set(labels)}')
        keep = {labels[label] for label in MAIN_LABELS}
        for terminal in list(tree.get_terminals()):
            if terminal.name not in keep: tree.prune(terminal)
        tips = [x.name for x in tree.get_terminals()]
        if len(tips) != 40: raise ValueError('Main roster must have 40 tips')
    inner = {x.name for x in t2.get_terminals()} - {'Telenomus_remus'}
    outer = set(tips) - inner
    counts = f1.load_counts()
    gh1, gh1_sources = load_extended_gh1()
    with (outdir/'gh1_cynipoid_sources.tsv').open('w') as handle:
        writer=csv.writer(handle,delimiter='\t'); writer.writerow(('species','GH1_e72','source'))
        for sp in sorted(inner): writer.writerow((sp,gh1[sp] if gh1.get(sp) is not None else 'NA',gh1_sources.get(sp,'not_scanned_in_qualifying_lineage')))
    corrected = f2.load_copy_unit_cells(f2.COPY_TABLE, f2.GENOME_ROSTER)
    famdata = {fam: corrected[fam] for fam in FAMILIES}
    chroms = {fam: f2.get_active_chromosomes(famdata[fam]) for fam in FAMILIES}
    if roster == 'main':
        chroms['CBM13'] = [ch for ch in chroms['CBM13'] if ch not in {'chr2','chr4','chr7'}]
    caps={fam:max([1] + [a.get(ch,0) for a in famdata[fam].values() for ch in chroms[fam] if isinstance(a.get(ch,0),int)] + [counts.get(sp,{}).get(fam,0) for sp in outer]) for fam in FAMILIES}
    gh1_cap=max([1] + [v for v in gh1.values() if isinstance(v,int)])
    eco = f2.load_ecology_states([x.name for x in t2.get_terminals()])
    with REMNANT.open() as handle:
        remnant = list(csv.DictReader(handle, delimiter='\t'))
    if len(remnant) != 1 or (remnant[0]['species'],remnant[0]['family'],remnant[0]['chromosome'],remnant[0]['mark']) != ('Parnips_nigripes','PL1_8','chr8','ψ'):
        raise ValueError('Unexpected remnant fixture')
    with ECOLOGY.open() as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            if row['species'] in eco:
                raise ValueError('Overlapping ecology fixtures')
            eco[row['species']] = {k: row[k] for k in ('gall_associated','herbivorous')}
    if not set(tips) <= set(eco):
        raise ValueError(f'Ecology missing: {set(tips)-set(eco)}')
    cells=[]
    display=[]
    for sp in tips:
        for fam in FAMILIES:
            if sp in outer:
                value=counts.get(sp,{}).get(fam,0)
                cells.append((sp,fam,'genome',str(value)))
                display.append((sp,fam,'genome',str(value),'count' if value else 'zero'))
            else:
                attrs=famdata[fam].get(sp,{})
                unplaced,_=f2.unresolved_display_state(attrs)
                pending=attrs.get('annotation')=='pending_no_hit'
                for chrom in chroms[fam]:
                    value=attrs.get(chrom,0)
                    cells.append((sp,fam,chrom,str(value)))
                    display.append((sp,fam,chrom,str(value),'count' if value else 'pending' if pending else 'unplaced' if unplaced else 'zero'))
    with (outdir/f'{roster}_combined_cells.tsv').open('w') as handle:
        writer=csv.writer(handle,delimiter='\t'); writer.writerow(('species','family','position','value')); writer.writerows(cells)
    with (outdir/f'{roster}_combined_display_states.tsv').open('w') as handle:
        writer=csv.writer(handle,delimiter='\t'); writer.writerow(('species','family','position','value','state')); writer.writerows(display)
    # Recheck every exported cell against the original source loaders.
    for sp,fam,pos,val in cells:
        expected = counts.get(sp,{}).get(fam,0) if pos=='genome' else famdata[fam].get(sp,{}).get(pos,0)
        assert val == str(expected)
    widths = {fam:max(1,len(chroms[fam])) for fam in FAMILIES}
    starts={}; x=-6.15
    for fam in FAMILIES:
        starts[fam]=x; x+=widths[fam]+0.35
    ghx=x+(1.2 if roster == 'main' else 0.8); gx=-10.1; hx=-7.8
    eco_width=1.4
    count_size=7.5 if roster == 'main' else 6.5
    size = (9,6.5) if layout == 'landscape' else ((7,7) if roster == 'main' else (7,9))
    fig,ax=plt.subplots(figsize=size)
    if roster == 'full':
        plt.rcParams['font.family'] = 'Nimbus Sans'
    # Match r11's x-span exactly: one chromosome unit has the same PDF-point width.
    r11_span = 59.65
    left = -34.0
    ax.set_xlim(left,left+r11_span); ax.set_ylim(len(tips)+(.45 if roster == 'full' else 1.0),-4.65); ax.axis('off')
    # The longest full-panel species name needs clearance before ecology.
    name_x = -25.0 if roster == 'full' else -23.7 if roster == 'main' and layout == 'portrait' else -21.3
    # Cladogram topology, aligned to the accepted display order.
    ymap={sp:i for i,sp in enumerate(tips)}
    def height(node):
        return 0 if node.is_terminal() else 1 + max(height(c) for c in node.clades)
    max_height=height(tree.root)
    tip_x=name_x-.3
    spans = bracket_spans(tips, inner)
    def bracket_x(label, lane):
        return -32.45+lane*1.3 + {1:.2,2:.55,3:.55}.get(lane,0)
    # Bracket labels sit left of their lines; reserve the existing .1-unit
    # clearance plus .22 for the bracket stroke and a visible .08-unit gap.
    root_x=max(bracket_x(label,lane) for label,_,_,lane in spans)+.4
    step=(tip_x-root_x)/max_height
    assert step>0
    assert all(height(c)<height(n) for n in tree.find_clades() for c in n.clades)
    def yy(node): return ymap[node.name] if node.is_terminal() else sum(yy(c) for c in node.clades)/len(node.clades)
    def xx(node): return tip_x-step*height(node)
    assert len({xx(n) for n in tree.get_terminals()}) == 1
    assert all(xx(c)-xx(n)>=step-1e-10 for n in tree.find_clades() for c in n.clades)
    point_step=step*size[0]*72/r11_span
    print(f'tree geometry {roster} {layout}: old=[{name_x-3.5:.4f},{name_x-1:.4f}] new=[{root_x:.4f},{tip_x:.4f}] height={max_height} step={step:.4f} units/{point_step:.4f} pt')
    if roster=='main' and layout=='landscape' and point_step<2.5:
        print(f'landscape per-level step below 2.5 pt: {point_step:.3f} pt')
    for node in tree.find_clades():
        if node.clades:
            ax.plot([xx(node),xx(node)],[min(yy(c) for c in node.clades),max(yy(c) for c in node.clades)],color='black',lw=.4)
            for child in node.clades: ax.plot([xx(node),xx(child)],[yy(child),yy(child)],color='black',lw=.4)
    # Every terminal node ends at tip_x; no separate leader is needed.
    for label, first, last, lane in spans:
        bx=bracket_x(label,lane)
        ancestor=tree.common_ancestor(tips[first], tips[last])
        # The root is the leftmost point of every tree segment, so this
        # excludes both horizontal and vertical segment intersections.
        assert bx+.32 < xx(tree.root), (label, 'bracket overlaps clade tree')
        ax.plot([bx,bx],[first-.35,last+.35],color='#888888',lw=1.2)
        label_y = (first+last)/2
        ax.text(bx-.22,label_y,{'Diplolepididae':'Diplolepidid.','Figitidae':'Figitid.'}.get(label,label),rotation=90,ha='right',va='center',fontsize=7.44,color='#444',weight='bold')
    colors={fam:'#%02x%02x%02x'%f2.FAMILY_COLORS[fam] for fam in FAMILIES}
    for fam in FAMILIES:
        st=starts[fam]; n=widths[fam]
        ax.text(st+n/2,-2.25,'B-lectin' if fam=='CBM13' else fam,ha='center',fontsize=6.0 if fam=='CBM13' and roster=='full' else 6.7,color=colors[fam],weight='bold')
        for j,ch in enumerate(chroms[fam]): ax.text(st+j+.5,-.85,ch.replace('chr',''),ha='center',fontsize=6.4,weight='bold')
    block_left=starts[FAMILIES[0]]; block_right=starts[FAMILIES[-1]]+widths[FAMILIES[-1]]
    ax.text((block_left+block_right)/2,-3.45,'Reference chromosome',ha='center',va='center',fontsize=6.5)
    ax.text(ghx+.82,-2.25,'GH1_e72',ha='center',fontsize=6.7)
    ax.text(gx+.02,-.65,'Gall-inhabiting',rotation=45,ha='left',va='bottom',fontsize=6.5); ax.text(hx+.02,-1.45,'Herbivorous',rotation=45,ha='left',va='bottom',fontsize=6.5)
    for i,sp in enumerate(tips):
        ax.text(name_x,i,f2.short_label(sp, {}),va='center',fontsize=6.5,style='italic',color='#555555' if sp in outer else 'black')
        for fam in FAMILIES:
            st=starts[fam]; n=widths[fam]; attrs=famdata[fam].get(sp,{})
            if sp in outer:
                value=counts.get(sp,{}).get(fam,0)
                ax.add_patch(Rectangle((st,i-.49),n-.005,.98,facecolor=cell_shade(value,caps[fam],fam),edgecolor='#ddd',lw=.2))
                if value: ax.text(st+n/2,i,str(value),ha='center',va='center',fontsize=count_size,fontfamily='Nimbus Sans',color=numeral_color(value,caps[fam]))
            else:
                unplaced,_=f2.unresolved_display_state(attrs)
                pending=attrs.get('annotation')=='pending_no_hit'
                vals=[attrs.get(ch,0) for ch in chroms[fam]]
                all_unplaced = attrs.get('total',0)>0 and attrs.get('placed_total',sum(vals))==0
                if all_unplaced:
                    value=attrs['total']
                    ax.add_patch(Rectangle((st,i-.49),n-.005,.98,facecolor=cell_shade(value,caps[fam],fam),edgecolor='#ddd',lw=.2))
                    ax.text(st+n/2,i,str(value),ha='center',va='center',fontsize=count_size,fontfamily='Nimbus Sans',color=numeral_color(value,caps[fam]))
                for j,val in enumerate(vals if not all_unplaced else []):
                    fill='#f3f3f3' if pending else f2.lerp_color(.20,fam) if unplaced and val==0 else cell_shade(val,caps[fam],fam)
                    ax.add_patch(Rectangle((st+j,i-.49),.995,.98,facecolor=fill,edgecolor='#222' if unplaced and val==0 else '#ddd',lw=.8 if unplaced and val==0 else .2,linestyle=(0,(2,1.5)) if unplaced and val==0 else 'solid'))
                    mark=str(val) if isinstance(val,int) and val>0 else ('…' if pending and j==n//2 else ('?' if unplaced else ''))
                    if mark: ax.text(st+j+.4975,i,mark,ha='center',va='center',fontsize=count_size,fontfamily='Nimbus Sans',color=numeral_color(val,caps[fam]) if val else '#999' if pending else '#222')
                is_galler=sp not in {'Telenomus_remus','Alloxysta_arcuata','Phaenoglyphis_villosa','Ganaspis_brasiliensis','Leptopilina_boulardi','Leptopilina_clavipes','Leptopilina_drosophilae','Leptopilina_heterotoma','Trichoplasta_sp','Trybliographa_sp','Parnips_nigripes','Periclistus_JH2016','Callirhytis_spRG'}
                asterisk=(fam,sp) in f2.ASTERISK_ABSENCES
                if attrs.get('total',0)==0 and not unplaced and not pending and is_galler and (fam!='CBM13' or asterisk):
                    ax.add_patch(Rectangle((st,i-.49),n-.005,.98,fill=False,edgecolor='black',lw=.7))
                    if asterisk: ax.text(st+n/2,i-.04,'*',ha='center',va='center',fontsize=6.2)
                if fam==remnant[0]['family'] and sp==remnant[0]['species']:
                    j=chroms[fam].index(remnant[0]['chromosome']); ax.text(st+j+.4975,i-.04,remnant[0]['mark'],ha='center',va='center',fontsize=6.2,color='#111')
        v=gh1.get(sp)
        ax.add_patch(Rectangle((ghx-.08,i-.49),1.80,.98,facecolor=f2.lerp_color(v/gh1_cap,'GH1_e72') if v is not None else '#f3f3f3',edgecolor='#ddd',lw=.2))
        ax.text(ghx+.82,i,'NA' if v is None else str(v),ha='center',va='center',fontsize=count_size,fontfamily='Nimbus Sans',color='#999' if v is None else numeral_color(v,gh1_cap))
        for key,cx,col in [('gall_associated',gx,'#8a7967'),('herbivorous',hx,'#4f4f4f')]:
            state=eco[sp][key]; ax.add_patch(Rectangle((cx,i-.49),eco_width,.98,facecolor=col if state=='yes' else '#ebebeb',edgecolor='#ddd',lw=.2))
            if state=='uncertain': ax.text(cx+eco_width/2,i,'?',ha='center',va='center',fontsize=6.0)
    ax.set_position([.005, .005, .99, .99])
    stem='tree_counts' if roster=='main' else 'tree_counts_full'
    fig.savefig(outdir/f'{stem}.pdf')
    fig.savefig(outdir/f'{stem}.png',dpi=300,facecolor='white')
    plt.close(fig)
    return len(tips),sum(sp in inner and gh1.get(sp) is not None for sp in tips),sum(sp in inner and gh1.get(sp) is None for sp in tips)

if __name__=='__main__':
    out=Path(os.environ['PCWDE_TREE_COUNTS_OUTPUT_DIR'])
    for roster,layout in [('main','portrait'),('full','portrait')]:
        print(roster,layout,build(out,roster,layout))
