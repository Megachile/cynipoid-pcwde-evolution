#!/usr/bin/env python3
"""
08_synteny_svg.py - SVG synteny figures from spreadsheets (CynOG v2 port)

Creates one SVG per block spreadsheet showing gene arrangements across genomes.
Consolidates tandem repeats and trims to N flanking genes per side.
Includes phylogenetic tree on left side.

Uses direct SVG generation (no matplotlib) for speed and low memory.

Input:
    - spreadsheet_*.tsv files (in SCRIPT_DIR)
    - ASTRAL species tree (Newick)
    - outputs/region_viz_subset.tsv (for CynOG lookup)

Output:
    - figures/synteny_{CynOG}_{chr#letter}.svg
"""

import pandas as pd
import os
from pathlib import Path
import re
from collections import defaultdict
import hashlib
import colorsys
from xml.sax.saxutils import escape

# Gem-tone PCWDE colors by CynOG family
PCWDE_COLORS = {
    'CynOGv3_000234': '#1B7837',  # emerald (GH5)
    'CynOGv3_000148': '#2166AC',  # sapphire (PL1)
    'CynOGv3_000168': '#B2182B',  # garnet (PL4)
}

BASE = Path(os.environ['PCWDE_PROJECT_ROOT']).resolve()
SCRIPT_DIR = Path(__file__).parent
OUTPUT_DIR = SCRIPT_DIR / 'outputs'
FIGURE_DIR = SCRIPT_DIR / 'figures'

VIZ_FILE = OUTPUT_DIR / 'region_viz_subset.tsv'
TREE_FILE = BASE / 'analysis/pancynipoid_phylogenomics/phase3_astral/astral_species_tree.tre'

# Layout parameters (pixels)
N_FLANKING = 4
GENE_W = 80
GENE_H = 22
GENE_GAP = 4
ROW_H = 28
LOCUS_W = 100
TREE_W = 120
LABEL_W = 145
MARGIN_TOP = 15
MARGIN_BOTTOM = 10
MARGIN_LEFT = 10
MARGIN_RIGHT = 10
FONT_SIZE = 10
FONT_SIZE_SMALL = 8


# === Tree parsing ===

def parse_newick(newick_str):
    newick_str = newick_str.strip().rstrip(';')

    def parse_node(s, pos=0):
        node = {'name': None, 'children': [], 'branch_length': 0}
        if pos >= len(s):
            return node, pos
        if s[pos] == '(':
            pos += 1
            while True:
                child, pos = parse_node(s, pos)
                node['children'].append(child)
                if pos >= len(s):
                    break
                if s[pos] == ',':
                    pos += 1
                elif s[pos] == ')':
                    pos += 1
                    break
        name_chars = []
        while pos < len(s) and s[pos] not in ',:();':
            name_chars.append(s[pos])
            pos += 1
        name_str = ''.join(name_chars)
        if ':' in name_str:
            name_str = name_str.split(':')[0]
        name_str = re.sub(r'^\d+\.?\d*/\d+$', '', name_str)
        node['name'] = name_str.strip() if name_str.strip() else None
        if pos < len(s) and s[pos] == ':':
            pos += 1
            bl_chars = []
            while pos < len(s) and s[pos] not in ',:();':
                bl_chars.append(s[pos])
                pos += 1
            try:
                node['branch_length'] = float(''.join(bl_chars))
            except ValueError:
                pass
        return node, pos

    tree, _ = parse_node(newick_str)
    return tree


def get_all_tips(node):
    if not node['children']:
        return [node['name']] if node['name'] else []
    tips = []
    for child in node['children']:
        tips.extend(get_all_tips(child))
    return tips


OUTGROUP = 'Telenomus_remus'


def reroot_at_outgroup(node, outgroup_name):
    from copy import deepcopy
    node = deepcopy(node)

    def find_path(n, target, path=None):
        if path is None:
            path = []
        path = path + [n]
        if not n['children'] and n['name'] == target:
            return path
        for child in n['children']:
            result = find_path(child, target, path)
            if result:
                return result
        return None

    path = find_path(node, outgroup_name)
    if not path or len(path) < 2:
        return node

    for i in range(len(path) - 2, 0, -1):
        child_node = path[i]
        parent_node = path[i - 1]
        parent_node['children'] = [
            c for c in parent_node['children'] if c is not child_node
        ]
        child_node['children'].append(parent_node)

    def collapse_unary(n):
        new_children = []
        for child in n['children']:
            child = collapse_unary(child)
            if child['children'] and len(child['children']) == 1 and not child['name']:
                grandchild = child['children'][0]
                grandchild['branch_length'] += child.get('branch_length', 0)
                new_children.append(grandchild)
            else:
                new_children.append(child)
        n['children'] = new_children
        return n

    new_root = collapse_unary(path[-2])
    return new_root


def prune_tree(node, keep_tips):
    keep_set = set(keep_tips)

    def prune_node(n):
        if not n['children']:
            if n['name'] in keep_set:
                return {'name': n['name'], 'children': [], 'branch_length': n['branch_length']}
            return None
        pruned_children = []
        for child in n['children']:
            pruned = prune_node(child)
            if pruned is not None:
                pruned_children.append(pruned)
        if not pruned_children:
            return None
        elif len(pruned_children) == 1:
            return pruned_children[0]
        else:
            return {'name': n['name'], 'children': pruned_children, 'branch_length': n['branch_length']}

    return prune_node(node)


def get_tree_depth(node, depth=0):
    if not node['children']:
        return depth
    return max(get_tree_depth(c, depth + 1) for c in node['children'])


def calculate_tree_segments(node, tip_y_map, x_scale):
    segments = []

    def get_coords(n, depth=0):
        x = depth * x_scale
        if not n['children']:
            y = tip_y_map.get(n['name'], 0)
            return x, y
        child_coords = []
        for child in n['children']:
            cx, cy = get_coords(child, depth + 1)
            child_coords.append((cx, cy))
        y = sum(c[1] for c in child_coords) / len(child_coords)
        for cx, cy in child_coords:
            segments.append((x, cy, cx, cy))
        min_cy = min(c[1] for c in child_coords)
        max_cy = max(c[1] for c in child_coords)
        segments.append((x, min_cy, x, max_cy))
        return x, y

    get_coords(node)
    return segments


def load_species_tree():
    if not TREE_FILE.exists():
        print(f"  Warning: tree file not found: {TREE_FILE}")
        return None
    with open(TREE_FILE) as f:
        newick = f.read().strip()
    newick = re.sub(r"'\[[^\]]*\]'", '', newick)
    tree = parse_newick(newick)
    tree = reroot_at_outgroup(tree, OUTGROUP)
    return tree


_SPECIES_TREE = None

def get_species_tree():
    global _SPECIES_TREE
    if _SPECIES_TREE is None:
        _SPECIES_TREE = load_species_tree()
    return _SPECIES_TREE


# === Annotation to color ===

ANNOT_TO_CYNOG = None

def build_annotation_to_cynog_map():
    viz_df = pd.read_csv(VIZ_FILE, sep='\t')
    annot_cynog_counts = defaultdict(lambda: defaultdict(int))
    for annot, cynog in zip(viz_df['annotation'].fillna(''), viz_df['cynog_id'].fillna('')):
        if annot and cynog:
            if len(str(annot)) > 40:
                annot = str(annot)[:37] + '...'
            annot_cynog_counts[str(annot)][str(cynog)] += 1
    annot_to_cynog = {}
    for annot, cynog_counts in annot_cynog_counts.items():
        annot_to_cynog[annot] = max(cynog_counts, key=cynog_counts.get)
    return annot_to_cynog


def get_cynog_for_annotation(annot):
    global ANNOT_TO_CYNOG
    if ANNOT_TO_CYNOG is None:
        ANNOT_TO_CYNOG = build_annotation_to_cynog_map()
    return ANNOT_TO_CYNOG.get(annot, annot)


def annotation_to_color(annot):
    if not annot or annot in ('', '[empty]', '[n/d]', '[hyp]'):
        return '#E8E8E8'
    cynog = get_cynog_for_annotation(annot)
    h = hashlib.md5(cynog.encode()).hexdigest()
    hue = int(h[:2], 16) / 255
    sat = 0.30 + (int(h[2:4], 16) / 255) * 0.20
    r, g, b = colorsys.hls_to_rgb(hue, 0.82, sat)
    return f'#{int(r*255):02x}{int(g*255):02x}{int(b*255):02x}'


# === Gene processing ===

def consolidate_tandems(genes):
    result = []
    current_annot = None
    current_count = 0
    for annot in genes:
        if not annot or (isinstance(annot, float) and pd.isna(annot)) or str(annot).strip() == '':
            continue
        annot = str(annot)
        if annot == current_annot:
            current_count += 1
        else:
            if current_annot is not None:
                result.append((current_annot, current_count))
            current_annot = annot
            current_count = 1
    if current_annot is not None:
        result.append((current_annot, current_count))
    return result


def truncate(text, max_len=12):
    if not text:
        return ''
    if len(text) > max_len:
        return text[:max_len - 2] + '..'
    return text


# === SVG generation ===

def svg_rect(x, y, w, h, fill, stroke='#333', stroke_width=0.5, rx=3):
    return f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{rx}" fill="{fill}" stroke="{stroke}" stroke-width="{stroke_width}"/>'


def svg_text(x, y, text, size=FONT_SIZE_SMALL, anchor='middle', style='', weight='normal', fill='#333'):
    escaped = escape(str(text))
    style_attr = f' font-style="{style}"' if style else ''
    return f'<text x="{x}" y="{y}" font-size="{size}" text-anchor="{anchor}" font-weight="{weight}" fill="{fill}"{style_attr} font-family="Arial, sans-serif">{escaped}</text>'


def svg_line(x1, y1, x2, y2, stroke='#333', width=0.8):
    return f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{stroke}" stroke-width="{width}"/>'


def render_block_svg(spreadsheet_file, output_file, cynog_family=None):
    df = pd.read_csv(spreadsheet_file, sep='\t')
    if len(df) == 0:
        return False

    fname = spreadsheet_file.stem
    # Match CynOGv3-style filenames: spreadsheet_CynOGv3_000148_chr2a etc.
    match = re.match(r'spreadsheet_(CynOGv3_\d+)_(.+)', fname)
    if not cynog_family and match:
        cynog_family = match.group(1)

    pcwde_color = PCWDE_COLORS.get(cynog_family, '#CC0000')

    upstream_cols = sorted(
        [c for c in df.columns if c.startswith('U') and c[1:].isdigit()],
        key=lambda x: int(x[1:]), reverse=True)
    downstream_cols = sorted(
        [c for c in df.columns if c.startswith('D') and c[1:].isdigit()],
        key=lambda x: int(x[1:]))

    # Exclude species with uncertain phylogenetic placement (LBA)
    SVG_EXCLUDE = {'Eschatocerus acacia'}
    df = df[~df['species'].isin(SVG_EXCLUDE)].reset_index(drop=True)
    if len(df) == 0:
        return False

    # Sort rows by registry tree_order (phylo_order column in spreadsheet)
    if 'phylo_order' in df.columns:
        df = df.sort_values('phylo_order').reset_index(drop=True)
    n_genomes = len(df)

    # Tree
    species_tree = get_species_tree()
    pruned_tree = None
    tree_segments = []
    has_tree = False

    if species_tree:
        tree_tips = get_all_tips(species_tree)
        sp_to_tip = {}
        for sp in df['species']:
            tree_name = str(sp).replace(' ', '_')
            if tree_name in tree_tips:
                sp_to_tip[sp] = tree_name

        keep = list(sp_to_tip.values())
        if keep:
            pruned_tree = prune_tree(species_tree, keep)
        if pruned_tree:
            has_tree = True
            df = df.reset_index(drop=True)
            n_genomes = len(df)

            tip_y_map = {}
            for row_idx, (_, row) in enumerate(df.iterrows()):
                tip = sp_to_tip.get(row['species'])
                if tip:
                    tip_y_map[tip] = MARGIN_TOP + row_idx * ROW_H + ROW_H / 2
            depth = get_tree_depth(pruned_tree)
            x_scale = TREE_W / (depth + 1) if depth > 0 else TREE_W
            tree_segments = calculate_tree_segments(
                pruned_tree, tip_y_map, x_scale)

    # Process rows
    processed = []
    for _, row in df.iterrows():
        up = [row.get(c, '') for c in upstream_cols]
        up_cons = consolidate_tandems(up)
        up_trim = up_cons[-N_FLANKING:] if len(up_cons) > N_FLANKING else up_cons

        dn = [row.get(c, '') for c in downstream_cols]
        dn_cons = consolidate_tandems(dn)
        dn_trim = dn_cons[:N_FLANKING] if len(dn_cons) > N_FLANKING else dn_cons

        processed.append({
            'species': str(row.get('species', '')),
            'status': str(row.get('status', '')),
            'locus': str(row.get('LOCUS', '')),
            'upstream': up_trim,
            'downstream': dn_trim,
        })

    max_up = max((len(r['upstream']) for r in processed), default=0)
    max_dn = max((len(r['downstream']) for r in processed), default=0)

    # Calculate dimensions
    tree_area = (MARGIN_LEFT + TREE_W + 5) if has_tree else MARGIN_LEFT
    genes_start_x = tree_area + LABEL_W + 5
    n_gene_slots = max_up + max_dn
    content_w = n_gene_slots * (GENE_W + GENE_GAP) + LOCUS_W + GENE_GAP
    total_w = genes_start_x + content_w + MARGIN_RIGHT
    total_h = MARGIN_TOP + n_genomes * ROW_H + MARGIN_BOTTOM

    elements = []

    # Background
    elements.append(f'<rect width="{total_w}" height="{total_h}" fill="white"/>')

    # Tree lines
    if has_tree:
        for x1, y1, x2, y2 in tree_segments:
            elements.append(svg_line(
                MARGIN_LEFT + x1, y1, MARGIN_LEFT + x2, y2))

    # Rows
    for row_idx, row_data in enumerate(processed):
        y = MARGIN_TOP + row_idx * ROW_H
        cy = y + ROW_H / 2

        # Species label
        label_x = tree_area + LABEL_W
        sp = row_data['species']
        if len(sp) > 30:
            sp = sp[:27] + '...'
        elements.append(svg_text(
            label_x, cy + 3, sp, size=FONT_SIZE, anchor='end', style='italic'))

        # Locus x position (fixed for alignment)
        locus_x = genes_start_x + max_up * (GENE_W + GENE_GAP)

        # Upstream genes (right-aligned to locus)
        n_up = len(row_data['upstream'])
        x = locus_x - n_up * (GENE_W + GENE_GAP)
        for annot, count in row_data['upstream']:
            color = annotation_to_color(annot)
            gy = cy - GENE_H / 2
            elements.append(svg_rect(x, gy, GENE_W, GENE_H, color, '#ccc', 0.3))
            x += GENE_W + GENE_GAP

        # Locus — gem-tone color for PCWDE, muted for tier2, gray for absent/n.d.
        locus = row_data['locus']
        if locus == '[n/d]':
            locus_fill, locus_stroke = '#CCCCCC', '#666666'
            locus_label = 'n.d.'
        elif locus == '[empty]':
            locus_fill, locus_stroke = 'none', '#999999'
            locus_label = 'absent'
        elif locus.startswith('[t2'):
            # Tier 2 — muted version of PCWDE color
            locus_fill, locus_stroke = '#AAAAAA', '#666666'
            locus_label = 'tier2'
        else:
            locus_fill, locus_stroke = pcwde_color, pcwde_color
            locus_label = truncate(locus, 14)

        gy = cy - GENE_H / 2
        elements.append(svg_rect(
            locus_x, gy, LOCUS_W, GENE_H, locus_fill, locus_stroke, 1.5))
        text_color = 'white' if locus_fill not in ('none', '#CCCCCC', '#E8E8E8', '#AAAAAA') else '#666'
        elements.append(svg_text(
            locus_x + LOCUS_W / 2, cy + 3, locus_label,
            size=FONT_SIZE, weight='bold', fill=text_color))

        # Downstream genes
        x = locus_x + LOCUS_W + GENE_GAP
        for annot, count in row_data['downstream']:
            color = annotation_to_color(annot)
            gy = cy - GENE_H / 2
            elements.append(svg_rect(x, gy, GENE_W, GENE_H, color, '#ccc', 0.3))
            x += GENE_W + GENE_GAP

    # Assemble SVG
    svg = (
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'width="{total_w}" height="{total_h}" '
        f'viewBox="0 0 {total_w} {total_h}">\n'
    )
    svg += '\n'.join(elements)
    svg += '\n</svg>\n'

    with open(output_file, 'w') as f:
        f.write(svg)

    return True


def main():
    print("=" * 60)
    print("Step 08: Synteny SVG Figures (CynOG v2)")
    print("=" * 60)
    print(f"Flanking genes per side: {N_FLANKING}")
    print(f"Output directory: {FIGURE_DIR}")

    FIGURE_DIR.mkdir(exist_ok=True)

    spreadsheets = sorted(SCRIPT_DIR.glob('spreadsheet_*.tsv'))
    print(f"\nFound {len(spreadsheets)} spreadsheet files")

    n_created = 0
    n_skipped = 0
    for ss_file in spreadsheets:
        fname = ss_file.stem
        match = re.match(r'spreadsheet_(CynOGv3_\d+)_(.+)', fname)
        if not match:
            print(f"  Skipping {fname} (unrecognized format)")
            n_skipped += 1
            continue

        family, block_id = match.groups()
        # Sanitize block_id for filename (? -> Q)
        safe_block = block_id.replace('?', 'Q')
        output_file = FIGURE_DIR / f"synteny_{family}_{safe_block}.svg"

        if render_block_svg(ss_file, output_file, cynog_family=family):
            n_created += 1
            print(f"  Created {output_file.name}")
        else:
            print(f"  Empty/failed: {fname}")

    print(f"\nCreated {n_created} SVG figures, skipped {n_skipped}")
    print(f"Output: {FIGURE_DIR}/")


if __name__ == '__main__':
    main()
