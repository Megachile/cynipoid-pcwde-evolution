#!/usr/bin/env python3
"""Render review trees and the crown-present background rate distribution."""
from __future__ import annotations
import argparse,csv
from collections import defaultdict
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from tree import parse_newick,prune
from score import node_id

COLORS={'GH5_2':'#1b7837','PL1_8':'#2166ac','PL4_4':'#b2182b'}
def rows(path):
    with open(path) as handle:yield from csv.DictReader(handle,delimiter='\t')
def draw_tree(family,tree,tip_states,events,out,title_prefix=''):
    tip_font = 10 if title_prefix else 8
    gain_font = 9 if title_prefix else 7
    sub=prune(tree,set(tip_states));order=[]
    def tips(node):
        if not node.children:order.append(node.name)
        else:
            for child in node.children:tips(child)
    tips(sub)
    ypos={name:len(order)-1-i for i,name in enumerate(order)}
    xy={}
    def coords(node,depth=0):
        if not node.children:y=ypos[node.name]
        else:
            for child in node.children:coords(child,depth+1)
            y=sum(xy[id(child)][1] for child in node.children)/len(node.children)
        xy[id(node)]=(depth,y)
    coords(sub)
    fig,ax=plt.subplots(figsize=(14,max(7,len(order)*.36+1.8)))
    ax.set_facecolor('#fffdf9')
    gains=defaultdict(list)
    for e in events:gains[e['child']].append(e['chromosome'])
    for node in list_nodes(sub):
        x,y=xy[id(node)]
        if node.children:
            ys=[xy[id(c)][1] for c in node.children]
            ax.plot([x,x],[min(ys),max(ys)],color='#687174',lw=1)
            for child in node.children:
                cx,cy=xy[id(child)]
                ax.plot([x,cx],[cy,cy],color='#687174',lw=1)
                marks=gains.get(node_id(child),[])
                if marks:
                    ax.scatter([(x+cx)/2],[cy],s=115+32*(len(marks)-1),color=COLORS[family],zorder=5)
                    ax.text((x+cx)/2,cy+.16,','.join(marks),fontsize=gain_font,color=COLORS[family],ha='center')
        else:
            ax.text(x+.15,y,node.name.replace('_',' ')+'  '+','.join(sorted(tip_states[node.name])),
                    va='center',fontsize=tip_font)
    max_depth=max(x for x,y in xy.values())
    ax.set_xlim(-.3,max_depth+6.5);ax.set_ylim(-1,len(order));ax.axis('off')
    ax.set_title(f'{title_prefix}{family}: retained chromosome gains on the species tree',loc='left',fontsize=13,pad=16)
    ax.text(0,-.75,'Filled circles mark DELTRAN 0→1 gains; labels give gained BK chromosomes. Tips show placed chromosome sets.',fontsize=8,color='#555')
    fig.tight_layout()
    if out is None:
        return fig
    fig.savefig(out,dpi=160,facecolor='white');plt.close(fig)

def list_nodes(root):
    yield root
    for child in root.children:yield from list_nodes(child)

def main():
    p=argparse.ArgumentParser();p.add_argument('--scores',type=Path,required=True)
    p.add_argument('--tips',type=Path,required=True);p.add_argument('--changes',type=Path,required=True)
    p.add_argument('--tree',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    scores=list(rows(a.scores));background=[float(r['deltran_rate']) for r in scores if r['kind']=='background']
    focal={r['orthogroup'].removeprefix('FOCAL:'):float(r['deltran_rate']) for r in scores if r['kind']=='focal'}
    fig,ax=plt.subplots(figsize=(9,5));ax.hist(background,bins=40,color='#8699a0',edgecolor='white')
    for family,color in COLORS.items():ax.axvline(focal[family],color=color,lw=2.3,label=family)
    ax.set_xlabel('Retained chromosome gains per placed genome');ax.set_ylabel('Crown-present orthogroups')
    ax.set_title('Chromosome-gain distribution');ax.legend(frameon=False)
    fig.tight_layout();fig.savefig(a.out/'background_histogram.png',dpi=180,facecolor='white');plt.close(fig)
    tree=parse_newick(a.tree.read_text())
    states=defaultdict(dict)
    for r in rows(a.tips):
        if r['orthogroup'].startswith('FOCAL:') and r['placed']=='1':
            states[r['orthogroup'].removeprefix('FOCAL:')][r['species']]=set(r['chromosomes'].split(','))
    changes=defaultdict(list)
    for r in rows(a.changes):
        if r['mode']=='deltran' and r['change']=='gain':changes[r['orthogroup'].removeprefix('FOCAL:')].append(r)
    for family in COLORS:draw_tree(family,tree,states[family],changes[family],a.out/f'{family}_tree.png')
    print('FIGURES_PASS',len(COLORS)+1)
if __name__=='__main__':main()
