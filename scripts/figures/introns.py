#!/usr/bin/env python3
"""Figure 3 combined intron-position renderer.

Top row (A-D): Intron position frequency stem plots along BK reference protein.
Bottom row (E-H+): Schematic panels for dominant intron sites with consensus
    flanking context and locus/chromosome distribution.

Order: PL1_8, PL4_4, GH5_2.
Data source: built intron groups, context and whole-genome background.

All positions mapped through pairwise alignment to a BK reference protein
per family (USE_TRIMMED=false — full-length coordinates).
"""

import sys, os
import pandas as pd
import numpy as np
import json
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyBboxPatch
from pathlib import Path
from collections import Counter

INPUT_VALUE = os.environ.get("PCWDE_INTRONS_INPUT_DIR", "")
if not INPUT_VALUE:
    raise SystemExit("PCWDE_INTRONS_INPUT_DIR is required")
DATA_DIR = Path(INPUT_VALUE).resolve()
if not DATA_DIR.is_dir():
    raise SystemExit(f"PCWDE_INTRONS_INPUT_DIR is not a directory: {DATA_DIR}")
REF_TSV = DATA_DIR / "ref_proteins.tsv"
BASELINE_INTRONS = DATA_DIR / "baseline_intron_details.tsv"

# Family display order and colors
FAMILIES = [
    ('PL1_8', '#2166AC'),
    ('PL4_4', '#B2182B'),
    ('GH5_2', '#1B7837'),
    
]
FAM_COLOR = {f: c for f, c in FAMILIES}

# How many dominant groups to show as schematics per family
N_SCHEMATIC = {'PL1_8': 2, 'PL4_4': 2, 'GH5_2': 1}

LOCUS_COLORS = {
    'chr6': '#264653', 'chr10': '#EF4444', 'chr3': '#2A9D8F',
    'unknown': '#D4D4D4',
    'chr1': '#DB2777', 'chr2': '#457B9D', 'chr7': '#A855F7',
    'chr8': '#3B82F6', 'chr9': '#F97316', 'chrunplaced': '#D4D4D4',
}

# Map the pipeline's "single_locus" sentinel to the family's actual chromosome
# (families with copies at one locus only — PL1_8 chr8, PL4_4 chr7 from synteny)
SINGLE_LOCUS_CHR = {'PL1_8': 'chr8', 'PL4_4': 'chr7'}


def lighten(hex_color, amount=0.75):
    rgb = mcolors.to_rgb(hex_color)
    return tuple(c + (1 - c) * amount for c in rgb)


def load_data():
    """Load all pipeline outputs.

    Reads:
      {fam}_groups.tsv           — clustered position groups (from 02b)
      {fam}_context_summary.tsv  — per-group consensus_pretty + cross-locus stats (from 04b)
      {fam}_context_rich.tsv     — per-column residue frequency (from 04b)
      {fam}_intron_positions.tsv — for n_intron_genes total
    """
    ref_df = pd.read_csv(REF_TSV, sep='\t')
    ref_len = dict(zip(ref_df['family'], ref_df['ref_len_full']))

    groups, ctx_sum, ctx_rich, n_intron_genes = {}, {}, {}, {}
    for fam, _ in FAMILIES:
        groups[fam] = pd.read_csv(DATA_DIR / f'{fam}_groups.tsv', sep='\t')

        sum_path = DATA_DIR / f'{fam}_context_summary.tsv'
        ctx_sum[fam] = pd.read_csv(sum_path, sep='\t') if sum_path.exists() and sum_path.stat().st_size > 1 else pd.DataFrame()

        rich_path = DATA_DIR / f'{fam}_context_rich.tsv'
        ctx_rich[fam] = pd.read_csv(rich_path, sep='\t') if rich_path.exists() and rich_path.stat().st_size > 1 else pd.DataFrame()

        pos_path = DATA_DIR / f'{fam}_intron_positions.tsv'
        pos_df = pd.read_csv(pos_path, sep='\t')
        n_intron_genes[fam] = pos_df[pos_df['ref_aa_pos'] > 0]['protein_id'].nunique()

    return ref_len, groups, ctx_sum, ctx_rich, n_intron_genes


def draw_stem_panel(ax, fam, groups_df, ref_len, n_genes_total, fcolor,
                    baseline_positions=None, show_left_label=True,
                    show_right_label=True, show_x_label=True):
    """Draw intron position frequency stem plot."""
    agg = groups_df.groupby('ref_aa_pos')['n_genes'].sum().reset_index()

    # Genome-wide baseline histogram
    if baseline_positions is not None:
        aa_pos = baseline_positions * ref_len
        n_bins = max(ref_len // 6, 30)
        counts_bl, edges = np.histogram(aa_pos, bins=n_bins, range=(0, ref_len))
        centers = (edges[:-1] + edges[1:]) / 2
        ax._baseline_centers = centers
        ax._baseline_counts = counts_bl

    stem_positions = agg['ref_aa_pos'].tolist()
    stem_counts = agg['n_genes'].tolist()

    if not stem_positions:
        ax.text(0.5, 0.5, 'No intron\nposition groups', transform=ax.transAxes,
                ha='center', va='center', fontsize=11, color='#999')
        ax.set_xlim(-5, ref_len + 10)
        ax.set_ylim(0, 10)
        return

    markerline, stemlines, _ = ax.stem(
        stem_positions, stem_counts, linefmt='-', markerfmt='o', basefmt=' ')
    markerline.set_color(fcolor)
    markerline.set_markersize(3.5)
    markerline.set_zorder(3)
    stemlines.set_color(fcolor)
    stemlines.set_alpha(0.7)
    stemlines.set_linewidth(1.2)

    ax.axhline(y=n_genes_total, color='#ccc', linestyle='--', linewidth=0.8, alpha=0.8)
    ax.text(0.97, 0.95, f'n={n_genes_total} mapped genes', fontsize=9, color='black',
            va='top', ha='right', transform=ax.transAxes)

    ymax = max(stem_counts) * 1.2 if stem_counts else 10
    ax.set_xlim(-5, ref_len + 10)
    ax.set_ylim(0, ymax)
    if show_x_label:
        ax.set_xlabel('Reference AA position', fontsize=11)
    if show_left_label:
        ax.set_ylabel('PCWDE genes', fontsize=11)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(labelsize=9)

    # Secondary y-axis: genome-wide baseline
    if hasattr(ax, '_baseline_counts') and ax._baseline_counts.max() > 0:
        ax2 = ax.twinx()
        bl_scaled = ax._baseline_counts / 1000
        ax2.fill_between(ax._baseline_centers, 0, bl_scaled,
                         color='#cccccc', alpha=0.45, step='mid')
        ax2.set_ylim(0, bl_scaled.max() * 1.2)
        if show_right_label:
            ax2.set_ylabel('All genes (thousands)', fontsize=10, color='black')
        ax2.tick_params(labelsize=8, colors='black')
        ax2.spines['right'].set_color('#ccc')
        ax2.spines['top'].set_visible(False)
        ax.set_zorder(ax2.get_zorder() + 1)
        ax.patch.set_visible(False)
        del ax._baseline_centers, ax._baseline_counts


def draw_schematic_panel(ax, fam, group_row, ctx_sum_df, ctx_rich_df):
    """Per-group intron schematic: ±12 AA consensus with conservation shading,
    numeric locus bar, cross-locus sharing stat."""
    fcolor = FAM_COLOR[fam]
    group_id = str(group_row['group_id'])
    ref_aa = int(group_row['ref_aa_pos'])
    n_genes = int(group_row['n_genes'])
    n_species = int(group_row['n_species'])
    phase_str = str(int(float(group_row['phase_mode']))) if pd.notna(group_row['phase_mode']) else '?'
    gt_ag_pct = group_row.get('gt_ag_pct', float('nan'))

    ax.set_title(f'{fam}  \u2014  AA {ref_aa}  (n={n_genes})', fontsize=11,
                 fontweight='bold', color=fcolor, pad=4)

    # --- ±6 AA consensus from context_summary.consensus_pretty (already has | marker) ---
    N_FLANK = 6
    sum_sub = ctx_sum_df[ctx_sum_df['group_id'] == group_id] if not ctx_sum_df.empty else pd.DataFrame()
    if not sum_sub.empty:
        pretty = str(sum_sub.iloc[0].get('consensus_pretty', ''))
        if '|' in pretty:
            cons_left, cons_right = pretty.split('|', 1)
        else:
            mid = len(pretty) // 2
            cons_left, cons_right = pretty[:mid], pretty[mid:]
        n_cross = int(sum_sub.iloc[0].get('n_species_cross_locus', 0))
    else:
        cons_left = '.' * N_FLANK
        cons_right = '.' * N_FLANK
        n_cross = 0
    cons_left = cons_left[-N_FLANK:].rjust(N_FLANK, '.')
    cons_right = cons_right[:N_FLANK].ljust(N_FLANK, '.')

    # Per-column conservation (top_pct) from context_rich
    pct_by_offset = {}
    if not ctx_rich_df.empty:
        rich_sub = ctx_rich_df[ctx_rich_df['group_id'] == group_id]
        for _, r in rich_sub.iterrows():
            pct_by_offset[int(r['col_offset'])] = float(r['top_pct'])

    # Draw cells (sized to fill the now-taller schematic row — chr bar removed,
    # so the cell row gets all the vertical real estate it can use).
    total_cells = N_FLANK * 2 + 1
    cell_w = 0.94 / total_cells
    start_x = 0.03
    cons_y = 0.50
    cell_h = 0.40

    def cell(x, aa, offset, is_splice):
        pct = pct_by_offset.get(offset, 0.0)
        if is_splice:
            bg = lighten(fcolor, 0.80); tc = fcolor; edge = fcolor; lw = 1.2
        elif aa in '.X_':
            bg = '#efefef'; tc = '#999'; edge = 'none'; lw = 0
        else:
            alpha = max(0.25, min(1.0, pct / 100.0))
            bg = lighten(fcolor, 1.0 - alpha * 0.70)
            tc = 'white' if alpha > 0.55 else '#1a1a1a'
            edge = 'none'; lw = 0
        rect = FancyBboxPatch((x, cons_y), cell_w * 0.94, cell_h,
                               boxstyle='square,pad=0',
                               facecolor=bg, edgecolor=edge,
                               linewidth=lw, transform=ax.transAxes)
        ax.add_patch(rect)
        display = '|' if is_splice else (aa if aa != '.' else '\u00B7')
        ax.text(x + cell_w * 0.47, cons_y + cell_h * 0.52,
                display, transform=ax.transAxes, ha='center', va='center',
                fontsize=14, fontweight='bold', color=tc, family='monospace')

    for i, aa in enumerate(cons_left):
        cell(start_x + i * cell_w, aa, i - N_FLANK, False)
    cell(start_x + N_FLANK * cell_w, '|', 0, True)
    for i, aa in enumerate(cons_right):
        cell(start_x + (N_FLANK + 1 + i) * cell_w, aa, i + 1, False)

    # --- Locus bar with numeric counts per segment ---
    raw_breakdown = json.loads(group_row['locus_class_breakdown'])
    # Remap single_locus to the family's real chromosome for display
    locus_breakdown = {}
    for k, v in raw_breakdown.items():
        if k == 'single_locus':
            k = SINGLE_LOCUS_CHR.get(fam, 'chr?')
        locus_breakdown[k] = locus_breakdown.get(k, 0) + v
    locus_counts = {k: v for k, v in locus_breakdown.items() if k not in ('unknown', 'other')}
    other_count = locus_breakdown.get('other', 0) + locus_breakdown.get('unknown', 0)
    total_locus = sum(locus_counts.values()) + other_count

    # Categorical chr list (replaces the proportional chr bar). Fixed grid
    # placement — each segment gets its own (col, row) cell with guaranteed
    # spacing, so labels never collide regardless of font width estimates.
    if total_locus > 0:
        locus_order = sorted(locus_counts.keys())
        if other_count > 0:
            locus_order.append('other')
            locus_counts['other'] = other_count
        segments = []
        for locus in locus_order:
            if locus in ('unplaced', 'other'):
                continue
            cnt = locus_counts[locus]
            color = LOCUS_COLORS.get(locus, '#888')
            segments.append((f"{locus}:{cnt}", color))
        if segments:
            n = len(segments)
            # Compact list layout. The old fixed grid still let wide
            # chromosome labels collide in panels with 4-5 loci.
            ncols = 1 if n <= 2 else 2
            nrows = (n + ncols - 1) // ncols
            chr_fontsize = 11 if nrows <= 2 else 9.5
            col_x = [0.31] if ncols == 1 else [0.16, 0.58]
            row_y0 = 0.42
            row_dy = 0.085
            for idx, (txt, color) in enumerate(segments):
                col = idx % ncols
                row = idx // ncols
                ax.text(col_x[col], row_y0 - row * row_dy, txt,
                        transform=ax.transAxes,
                        ha='left', va='center', fontsize=chr_fontsize,
                        color=color, fontweight='bold')

    # Stats line: total species + cross-locus sharing only
    stats = f'{n_species} spp'
    if n_cross > 0:
        stats += f'  \u00B7  {n_cross} spp share across loci'
    ax.text(0.5, 0.10, stats, transform=ax.transAxes,
            ha='center', va='bottom', fontsize=11, color='#222')

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')


def main():
    np.random.seed(42)

    output_value = os.environ.get("PCWDE_INTRONS_OUTPUT_DIR", "")
    if not output_value:
        raise SystemExit("PCWDE_INTRONS_OUTPUT_DIR is required")
    output_dir = Path(output_value)
    if output_dir.exists() and any(output_dir.iterdir()):
        raise SystemExit(f"refusing non-empty output directory: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)

    ref_len, groups, ctx_sum, ctx_rich, n_intron_genes = load_data()

    # Load genome-wide baseline
    if not BASELINE_INTRONS.is_file():
        raise SystemExit(f"missing built whole-genome intron baseline: {BASELINE_INTRONS}")
    bl = pd.read_csv(BASELINE_INTRONS, sep='\t', usecols=['relative_cds_position'])
    baseline_pos = bl['relative_cds_position'].dropna().values
    if not len(baseline_pos):
        raise SystemExit("built whole-genome intron baseline is empty")
    print(f"Loaded {len(baseline_pos):,} baseline introns for background")

    # Identify dominant groups for schematic panels
    targets = []
    for fam, _ in FAMILIES:
        n_show = N_SCHEMATIC[fam]
        if n_show > 0 and not groups[fam].empty:
            top = groups[fam].nlargest(n_show, 'n_genes').sort_values('ref_aa_pos')
            for _, row in top.iterrows():
                targets.append((fam, row))

    print("\nSchematic targets:")
    for fam, row in targets:
        print(f"  {fam} group {row['group_id']}: "
              f"ref_aa={int(row['ref_aa_pos'])}, n_genes={int(row['n_genes'])}, "
              f"n_species={int(row['n_species'])}")

    n_schematics = len(targets)

    # --- Figure layout (original landscape): 2 rows × 15 cols.  Stems on top,
    # schematic strips on bottom, family-aligned 5-col blocks. The categorical
    # chr text (replacing the proportional chr bar) is the only behavior change
    # vs the pre-headless version. ---
    fig = plt.figure(figsize=(18, 8))
    fig.patch.set_facecolor('white')

    gs = gridspec.GridSpec(2, 15, figure=fig, hspace=0.40, wspace=1.5,
                           height_ratios=[1, 1.45])
    fig.subplots_adjust(left=0.05, right=0.96, top=0.92, bottom=0.08)

    # Top row: one stem plot per family, each exactly 5 cols wide
    ax_stems = []
    for i in range(len(FAMILIES)):
        ax_stems.append(fig.add_subplot(gs[0, i * 5:(i + 1) * 5]))

    # Bottom row: schematic panels on the 15-col grid
    ax_schematics = []
    if n_schematics > 0:
        cols_per = 15 // max(n_schematics, 1)
        for i in range(n_schematics):
            c0 = i * cols_per
            c1 = min(c0 + cols_per, 15) if i < n_schematics - 1 else 15
            ax_schematics.append(fig.add_subplot(gs[1, c0:c1]))

    panel_label_idx = 0
    panel_labels = "ABCDEFGHIJKL"

    for i, (fam, fcolor) in enumerate(FAMILIES):
        ax = ax_stems[i]
        ax.set_facecolor('white')
        ax.set_title(fam, fontsize=15, fontweight='bold', color=fcolor)
        draw_stem_panel(ax, fam, groups[fam], ref_len[fam], n_intron_genes[fam],
                        fcolor, baseline_positions=baseline_pos,
                        show_left_label=(i == 0),
                        show_right_label=(i == len(FAMILIES) - 1),
                        show_x_label=True)
        ax.text(-0.02, 1.10, panel_labels[panel_label_idx],
                transform=ax.transAxes, fontsize=14, fontweight='bold', va='top')
        panel_label_idx += 1

    for j, (fam, group_row) in enumerate(targets):
        ax = ax_schematics[j]
        ax.set_facecolor('white')
        draw_schematic_panel(ax, fam, group_row, ctx_sum[fam], ctx_rich[fam])
        ax.text(0.0, 1.08, panel_labels[panel_label_idx],
                transform=ax.transAxes, fontsize=14, fontweight='bold', va='top')
        panel_label_idx += 1

    out_svg = output_dir / 'introns.svg'
    out_pdf = output_dir / 'introns.pdf'
    out_png = output_dir / 'introns.png'
    plt.savefig(out_svg, format='svg', bbox_inches='tight', facecolor='white')
    plt.savefig(out_pdf, bbox_inches='tight', facecolor='white')
    plt.savefig(out_png, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"\nSaved: {out_svg.name}, {out_pdf.name}, {out_png.name}")


if __name__ == '__main__':
    main()
