#!/usr/bin/env python3
"""Figure 8: expression by tissue.

Rows: DC (D. cinerosa transcriptome), BK (B. kinseyi genome)
Columns: PL1_8, PL4_4, GH5_2
GH5_2 shows grouped boxplots per locus (chr10_venom, chr10_larval, chr6).

Data sources are the supplied expression tables, family membership, chromosome
assignments and tissue crosswalks. The expression.py entry point binds
these inputs from PCWDE_FIGURE_DATA_ROOT/figure8 and configures this module.
"""

import sys, os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from pathlib import Path

INPUT_VALUE = os.environ.get('PCWDE_EXPRESSION_INPUT_DIR', '')
if not INPUT_VALUE:
    raise SystemExit('PCWDE_EXPRESSION_INPUT_DIR is required')
INPUT_DIR = Path(INPUT_VALUE).resolve()
if not INPUT_DIR.is_dir():
    raise SystemExit(f'PCWDE_EXPRESSION_INPUT_DIR is not a directory: {INPUT_DIR}')
LOCUS_CLASS = INPUT_DIR / 'gh5_locus_classification.tsv'
SHARED_VALUE = os.environ.get('PCWDE_FIGURE_SHARED_INPUT_DIR', '')
if not SHARED_VALUE:
    raise SystemExit('PCWDE_FIGURE_SHARED_INPUT_DIR is required')
SHARED_INPUT_DIR = Path(SHARED_VALUE).resolve()

FAMILIES = [
    ('PL1_8', '#2166AC'),   # sapphire
    ('PL4_4', '#B2182B'),   # garnet
    ('GH5_2', '#1B7837'),   # emerald
]

# GH5_2 locus display order and colors (emerald shades)
GH5_LOCUS_ORDER = ['chr10_larval', 'chr10_venom', 'chr6']
GH5_LOCUS_STYLE = {
    'chr10_larval': ('#0d5c23', 'chr10 (larval)'),
    'chr10_venom':  ('#1B7837', 'chr10 (venom)'),
    'chr6':         ('#4ca96b', 'chr6 (salivary)'),
}

# PL1_8 lineage breakdown (homology-inferred DC tissue lineage; sapphire shades).
# PL1_8 is a single-chromosome (chr8) tandem array, so unlike GH5_2 there is no
# chromosome split to group on — paralogs are grouped by the nearest
# DC-tissue-labelled neighbour in the PL1_8 locus tree (crosswalk below).
# Head lineage is intentionally excluded (negligible expression in the shown
# tissues; the DC Head tissue column is not part of the figure axis).
PL18_XWALK = INPUT_DIR / 'pl1_8_homology_tissue_crosswalk.tsv'
PL18_LINEAGE_ORDER = ['Salivary', 'Gut', 'Malpighian', 'Venom']
PL18_LINEAGE_STYLE = {
    'Salivary':   ('#103a66', 'Salivary'),
    'Gut':        ('#2166AC', 'Gut'),
    'Malpighian': ('#4a90c4', 'Malpighian'),
    'Venom':      ('#86b6da', 'Venom'),
}
# A lineage is shown in a species' row only if at least one of its paralogs
# reaches this peak tissue-mean TPM — drops paralogs with no meaningful
# expression in any tissue (e.g. the silent BK venom-lineage ortholog).
PL18_DROP_TPM = 10.0

# Vertically inherited beta-glucosidase candidate for GH5_2 product processing
NONHGT_ORDER = ['GH1_e72']
NONHGT_STYLE = {
    'GH1_e72': ('#555555', 'GH1_e72'),
}
NONHGT_LABEL = 'GH1_e72'

# CBM13 / B-lectin domain HGT family (previewed 4th HGT family; full analysis shelved)
CBM13_COLOR = '#D98C00'  # amber
BLECTIN_MEMBERSHIP = INPUT_DIR / 'b_lectin_membership.tsv'

# DC tissues
DC_TISSUES = {
    'Gut':     ['tpm_DCL_G3', 'tpm_DCL_G5', 'tpm_DCL_G1L', 'tpm_DCL_G2L'],
    'Malp.':   ['tpm_DCL_M1E', 'tpm_DCL_M1L', 'tpm_DCL_M2L', 'tpm_DCL_M4'],
    'Saliv.':  ['tpm_DCL_S8', 'tpm_DCL_S9'],
    'Venom':   ['tpm_DC_VGX_01', 'tpm_DC_VGX_02', 'tpm_DC_VGX_03', 'tpm_DC_VGX_04'],
    'Carc.':   ['tpm_DC_CRX_01', 'tpm_DC_CRX_02', 'tpm_DC_CRX_03', 'tpm_DC_CRX_04'],
    'Acc.':    ['tpm_DC_AGX_01', 'tpm_DC_AGX_02', 'tpm_DC_AGX_03', 'tpm_DC_AGX_04'],
}
DC_ORDER = ['Gut', 'Malp.', 'Saliv.', 'Venom', 'Carc.', 'Acc.']
DC_LARVAL_N = 3  # first 3 are retained larval tissues

# BK tissues
BK_TISSUES = {
    'Larvae': ['tpm_ALA_BK_1', 'tpm_ALA_BK_3', 'tpm_ALA_BK_4', 'tpm_ALA_BK_5', 'tpm_ALA_BK_8'],
    'Venom':  ['tpm_BKS_VGX_02', 'tpm_BKS_VGX_03', 'tpm_BKS_VGX_04', 'tpm_BKS_VGX_06', 'tpm_BKS_VGX_08'],
    'Acc.':   ['tpm_BKS_AGX_02', 'tpm_BKS_AGX_03', 'tpm_BKS_AGX_04', 'tpm_BKS_AGX_06', 'tpm_BKS_AGX_08'],
    'Carc.':  ['tpm_BKS_CRX_02', 'tpm_BKS_CRX_03', 'tpm_BKS_CRX_04', 'tpm_BKS_CRX_06', 'tpm_BKS_CRX_08'],
}
BK_ORDER = ['Larvae', 'Venom', 'Acc.', 'Carc.']
BK_LARVAL_N = 1  # first 1 is larval


def lighten(hex_color, amount=0.75):
    rgb = mcolors.to_rgb(hex_color)
    return tuple(c + (1 - c) * amount for c in rgb)


def load_data():
    """Load family membership, TPMs, and GH5_2 locus classifications."""
    meta = pd.read_csv(Path(os.environ['PCWDE_EXPRESSION_MEMBER_METADATA']), sep='\t')
    dc_ann = pd.read_csv(INPUT_DIR / 'dc_annotation.tsv',
                         sep='\t', low_memory=False)
    bk_ann = pd.read_csv(INPUT_DIR / 'bk_annotation.tsv',
                         sep='\t', low_memory=False)
    p2c = pd.read_csv(SHARED_INPUT_DIR / 'protein_to_chromosome.tsv',
                      sep='\t')

    # GH5_2 DC locus classification (DIAMOND BLAST + diagnostic amino acid positions)
    cls_df = pd.read_csv(LOCUS_CLASS, sep='\t', comment='#',
                         names=['protein_id', 'chr', 'classification',
                                'v_score', 'l_score', 'covered', 'profile'])
    dc_cls_map = cls_df.set_index('protein_id')['classification'].to_dict()

    # BK chromosome map (strip version suffix to match p2c)
    bk_chr_map = p2c.set_index('protein_id')['bk_chr'].to_dict()

    # PL1_8 homology-inferred lineage maps (DC = own peak tissue / lineage anchor;
    # BK = nearest DC tissue lineage in the locus tree)
    xw = pd.read_csv(PL18_XWALK, sep='\t')
    xw = xw[xw['family'] == 'PL1_8']
    pl18_dc_lineage = dict(zip(xw['nearest_dc_protein_id'], xw['homology_inferred_tissue']))
    pl18_bk_lineage = dict(zip(xw['bk_protein_id'], xw['homology_inferred_tissue']))

    family_data = {}
    for family, _color in FAMILIES:
        # DC members from metadata
        dc_pids = set(meta[(meta['family'] == family) &
                           (meta['source_file'] == 'DC')]['protein_id'])
        dc_rows = dc_ann[dc_ann['protein_id'].isin(dc_pids)].copy()
        dc_rows['_species'] = 'DC'

        # BK members from metadata (metadata uses XP_NNN, annotation uses XP_NNN.1)
        bk_pids = set(meta[(meta['family'] == family) &
                           (meta['source_file'] == 'BK_NCBI_genome')]['protein_id'])
        bk_ann['_pid_stripped'] = bk_ann['protein_id'].str.replace(r'\.\d+$', '', regex=True)
        bk_rows = bk_ann[bk_ann['_pid_stripped'].isin(bk_pids)].copy()
        bk_rows['_species'] = 'BK'

        if family == 'GH5_2':
            # DC: sequence-based classification from diagnostic positions
            dc_rows['classification'] = dc_rows['protein_id'].map(dc_cls_map).fillna('uncertain')

            # BK: chromosome-based classification
            bk_rows['_chr'] = bk_rows['_pid_stripped'].map(bk_chr_map).fillna('unknown').astype(str)
            bk_rows['classification'] = 'chr6'
            is_chr10 = bk_rows['_chr'] == '10'
            # Within chr10, split venom vs larval by expression
            vgx = [c for c in bk_rows.columns if c.startswith('tpm_BKS_VGX_')]
            ala = [c for c in bk_rows.columns if c.startswith('tpm_ALA_BK_')]
            bk_rows.loc[is_chr10, 'classification'] = 'chr10_larval'
            bk_rows.loc[is_chr10 & (bk_rows[vgx].max(axis=1) > bk_rows[ala].max(axis=1)),
                        'classification'] = 'chr10_venom'

        elif family == 'PL1_8':
            # Homology-inferred DC tissue lineage (single-chr8 tandem array)
            dc_rows['classification'] = dc_rows['protein_id'].map(pl18_dc_lineage).fillna('uncertain')
            bk_rows['classification'] = bk_rows['_pid_stripped'].map(pl18_bk_lineage).fillna('uncertain')

        combined = pd.concat([dc_rows, bk_rows], ignore_index=True)
        family_data[family] = combined
        n_dc = len(dc_rows)
        n_bk = len(bk_rows)
        print(f"{family}: {n_dc} DC + {n_bk} BK = {n_dc + n_bk} proteins")
        if family == 'GH5_2':
            for sp in ['DC', 'BK']:
                sub = combined[combined['_species'] == sp]
                counts = sub['classification'].value_counts().to_dict()
                print(f"  {sp} loci: {counts}")

    # ── Non-HGT reference controls: match dbcan_families for each subfamily ──
    def classify_nonhgt(dbcan_str):
        if not isinstance(dbcan_str, str):
            return None
        tokens = set(dbcan_str.split('|'))
        for sub in NONHGT_ORDER:
            if sub in tokens:
                return sub
        return None

    dc_nh = dc_ann.copy()
    dc_nh['classification'] = dc_nh['dbcan_families'].map(classify_nonhgt)
    dc_nh = dc_nh[dc_nh['classification'].notna()].copy()
    dc_nh['_species'] = 'DC'

    bk_nh = bk_ann.copy()
    bk_nh['classification'] = bk_nh['dbcan_families'].map(classify_nonhgt)
    bk_nh = bk_nh[bk_nh['classification'].notna()].copy()
    bk_nh['_species'] = 'BK'

    nonhgt = pd.concat([dc_nh, bk_nh], ignore_index=True)
    family_data['nonHGT'] = nonhgt
    dc_counts = dc_nh['classification'].value_counts().to_dict()
    bk_counts = bk_nh['classification'].value_counts().to_dict()
    print(f"non-HGT: DC {dc_counts}, BK {bk_counts}")

    # ── CBM13 / B-lectin domain HGT family: join curated membership to annotation TSVs ──
    blectin_members = pd.read_csv(BLECTIN_MEMBERSHIP, sep='\t')
    dc_blectin_ids = set(
        blectin_members.loc[
            blectin_members['species'] == 'Disholcaspis cinerosa', 'protein_id'
        ]
    )
    bk_blectin_ids = set(
        blectin_members.loc[
            blectin_members['species'] == 'Belonocnema kinseyi', 'protein_id'
        ]
    )

    dc_cbm = dc_ann[dc_ann['protein_id'].isin(dc_blectin_ids)].copy()
    dc_cbm['_species'] = 'DC'

    # BK membership uses short IDs (XP_033211071); BK.tsv uses XP_033211071.1
    if '_pid_stripped' not in bk_ann.columns:
        bk_ann['_pid_stripped'] = bk_ann['protein_id'].str.replace(r'\.\d+$', '', regex=True)
    bk_cbm = bk_ann[bk_ann['_pid_stripped'].isin(bk_blectin_ids)].copy()
    bk_cbm['_species'] = 'BK'

    cbm13 = pd.concat([dc_cbm, bk_cbm], ignore_index=True)
    family_data['CBM13'] = cbm13
    print(f"CBM13: DC {len(dc_cbm)} + BK {len(bk_cbm)} = {len(cbm13)} proteins")

    return family_data


def get_replicate_totals(df, species, tissue_samples, locus=None):
    """Sum TPM across copies for each replicate. Optionally filter by locus."""
    sp_df = df[df['_species'] == species]
    if locus is not None:
        sp_df = sp_df[sp_df['classification'] == locus]
    available = [s for s in tissue_samples if s in df.columns]
    if not available or len(sp_df) == 0:
        return np.array([])
    return sp_df[available].sum(axis=0).values


def draw_simple_panel(ax, df, species, tissue_map, tissue_order, fcolor):
    """Single-locus boxplot panel (PL1_8, PL4_4)."""
    light = lighten(fcolor, 0.75)
    data = []
    for tissue in tissue_order:
        vals = get_replicate_totals(df, species, tissue_map[tissue])
        data.append(vals if len(vals) > 0 else np.array([0]))

    bp = ax.boxplot(data, tick_labels=tissue_order, patch_artist=True,
                    widths=0.55,
                    medianprops=dict(color=fcolor, linewidth=1.5),
                    whiskerprops=dict(color='#aaa'),
                    capprops=dict(color='#aaa'),
                    flierprops=dict(marker=''))
    for patch in bp['boxes']:
        patch.set_facecolor(light)
        patch.set_edgecolor(fcolor)
        patch.set_linewidth(0.8)

    for i, vals in enumerate(data):
        if len(vals) > 0:
            x = np.random.normal(i + 1, 0.05, size=len(vals))
            ax.scatter(x, vals, color=fcolor, s=22, alpha=0.8,
                      edgecolor='white', linewidth=0.3, zorder=3)


# Short chr labels for inline box annotation
CHR_SHORT = {
    'chr10_larval': 'c10L',
    'chr10_venom':  'c10V',
    'chr6':         'c6',
}
TPM_THRESHOLD = 5000  # only label boxes where max TPM exceeds this


def draw_grouped_panel(ax, df, species, tissue_map, tissue_order, loci,
                       inferred=False, label_ymax=None, show_labels=True,
                       style_map=None, short_map=None):
    """Multi-locus grouped boxplot panel (GH5_2 loci / PL1_8 lineages).

    label_ymax: if set, cap chr label y-position at this value (for broken axis).
    show_labels: if False, skip chr labels entirely.
    style_map/short_map: override the GH5_2 locus styling (e.g. PL1_8 lineages).
    """
    if style_map is None:
        style_map = GH5_LOCUS_STYLE
    if short_map is None:
        short_map = CHR_SHORT
    n_loci = len(loci)
    if n_loci == 1:
        offsets = [0]
        bw = 0.5
    elif n_loci == 2:
        offsets = [-0.18, 0.18]
        bw = 0.28
    elif n_loci == 3:
        offsets = [-0.25, 0, 0.25]
        bw = 0.2
    else:
        span = 0.74
        offsets = list(np.linspace(-span / 2, span / 2, n_loci))
        bw = span / n_loci * 0.85

    for li, locus in enumerate(loci):
        style = style_map.get(locus, ('#1B7837', locus))
        fcolor = style[0]
        light = lighten(fcolor, 0.65)

        positions = [t + 1 + offsets[li] for t in range(len(tissue_order))]
        data = []
        for tissue in tissue_order:
            vals = get_replicate_totals(df, species, tissue_map[tissue], locus)
            data.append(vals if len(vals) > 0 else np.array([0]))

        bp = ax.boxplot(data, positions=positions, widths=bw, patch_artist=True,
                        medianprops=dict(color=fcolor, linewidth=1.2),
                        whiskerprops=dict(color='#bbb'),
                        capprops=dict(color='#bbb'),
                        flierprops=dict(marker=''))
        for patch in bp['boxes']:
            patch.set_facecolor(light)
            patch.set_edgecolor(fcolor)
            patch.set_linewidth(0.7)

        for i, vals in enumerate(data):
            if len(vals) > 0:
                x = np.random.normal(positions[i], 0.03, size=len(vals))
                ax.scatter(x, vals, color=fcolor, s=18, alpha=0.8,
                          edgecolor='white', linewidth=0.2, zorder=3)

        # Add chr labels on boxes with high expression
        if show_labels:
            short = short_map.get(locus, locus)
            if inferred:
                short = short + '*'
            fsize = 10 if len(short.rstrip('*')) <= 2 else 9
            for i, vals in enumerate(data):
                if len(vals) > 0 and vals.max() > TPM_THRESHOLD:
                    top = vals.max()
                    if label_ymax is not None:
                        top = min(top, label_ymax)
                    ax.text(positions[i], top * 1.08, short,
                            ha='center', va='bottom',
                            fontsize=fsize, color=fcolor, fontweight='bold')

    ax.set_xticks(range(1, len(tissue_order) + 1))
    ax.set_xticklabels(tissue_order)


def add_stage_brackets(ax, n_larval, n_total, y_line=-0.60, y_text=-0.65):
    """Add 'larval' / 'adult' brackets below x-axis labels."""
    trans = ax.get_xaxis_transform()
    if n_larval > 0:
        x0, x1 = 0.6, n_larval + 0.4
        ax.plot([x0, x1], [y_line, y_line],
                transform=trans, clip_on=False, color='black', linewidth=1)
        ax.text((x0 + x1) / 2, y_text, 'larval',
                transform=trans, ha='center', va='top', fontsize=9, color='black')
    if n_larval < n_total:
        x0, x1 = n_larval + 0.6, n_total + 0.4
        ax.plot([x0, x1], [y_line, y_line],
                transform=trans, clip_on=False, color='black', linewidth=1)
        ax.text((x0 + x1) / 2, y_text, 'adult',
                transform=trans, ha='center', va='top', fontsize=9, color='black')




def style_ax(ax):
    ax.set_facecolor('white')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_ylim(bottom=0)
    ax.tick_params(axis='x', rotation=45, labelsize=9)
    ax.tick_params(axis='y', labelsize=9)


BREAK_LO = 25000   # bottom axis goes 0 → BREAK_LO
BREAK_HI_PAD = 10000  # top axis goes (data_max - pad) → (data_max + pad)


def add_break_marks(ax_top, ax_bot):
    """Draw diagonal break marks on the y axis between two stacked axes.

    Only the left (y-axis) side: the right spine is hidden, so right-edge
    marks would float beside the panel.
    """
    d = 0.015
    kwargs = dict(transform=ax_top.transAxes, color='black',
                  clip_on=False, linewidth=0.8)
    ax_top.plot((-d, +d), (-d, +d), **kwargs)
    kwargs['transform'] = ax_bot.transAxes
    ax_bot.plot((-d, +d), (1 - d, 1 + d), **kwargs)


def draw_nonhgt_panel(ax, df, species, tissue_map, tissue_order):
    """Grouped boxplot panel for non-HGT reference enzymes.

    GH1_e72 paralogs are pooled (summed) per replicate and shown as a
    single vertically inherited beta-glucosidase comparison.
    """
    sub_order = [s for s in NONHGT_ORDER
                 if s in df[df['_species'] == species]['classification'].unique()]
    n_sub = len(sub_order)
    if n_sub == 0:
        return
    if n_sub == 1:
        offsets = [0]; bw = 0.5
    elif n_sub == 2:
        offsets = [-0.18, 0.18]; bw = 0.28
    else:
        offsets = [-0.25, 0, 0.25]; bw = 0.2

    for li, sub in enumerate(sub_order):
        fcolor, _ = NONHGT_STYLE[sub]
        light = lighten(fcolor, 0.6)
        positions = [t + 1 + offsets[li] for t in range(len(tissue_order))]
        data = []
        for tissue in tissue_order:
            vals = get_replicate_totals(df, species, tissue_map[tissue], locus=sub)
            data.append(vals if len(vals) > 0 else np.array([0]))
        bp = ax.boxplot(data, positions=positions, widths=bw, patch_artist=True,
                        medianprops=dict(color=fcolor, linewidth=1.2),
                        whiskerprops=dict(color='#bbb'),
                        capprops=dict(color='#bbb'),
                        flierprops=dict(marker=''))
        for patch in bp['boxes']:
            patch.set_facecolor(light)
            patch.set_edgecolor(fcolor)
            patch.set_linewidth(0.7)
        for i, vals in enumerate(data):
            if len(vals) > 0:
                x = np.random.normal(positions[i], 0.03, size=len(vals))
                ax.scatter(x, vals, color=fcolor, s=18, alpha=0.8,
                           edgecolor='white', linewidth=0.2, zorder=3)
    ax.set_xticks(range(1, len(tissue_order) + 1))
    ax.set_xticklabels(tissue_order)


def draw_nonhgt_legend(ax):
    """Small inline legend identifying GH1_e72."""
    handles = []
    labels = []
    from matplotlib.patches import Patch
    for sub in NONHGT_ORDER:
        fcolor, label = NONHGT_STYLE[sub]
        handles.append(Patch(facecolor=lighten(fcolor, 0.6), edgecolor=fcolor, linewidth=0.8))
        labels.append(label)
    ax.legend(handles, labels, loc='upper right', fontsize=8,
              frameon=False, handlelength=1.0, handletextpad=0.4,
              labelspacing=0.2, borderaxespad=0.3)


def pl18_lineages_present(df, species, tissue_map):
    """PL1_8 lineages shown in a species' row: those with >=1 paralog reaching
    PL18_DROP_TPM peak tissue-mean TPM. Drops paralogs with no meaningful
    expression in any tissue (e.g. the silent BK venom-lineage ortholog)."""
    present = []
    for lin in PL18_LINEAGE_ORDER:
        sub = df[(df['_species'] == species) & (df['classification'] == lin)]
        if len(sub) == 0:
            continue
        peak = 0.0
        for cols in tissue_map.values():
            avail = [c for c in cols if c in df.columns]
            if avail:
                peak = max(peak, sub[avail].sum(axis=0).mean())
        if peak >= PL18_DROP_TPM:
            present.append(lin)
    return present


def draw_pl18_legend(ax, lineages):
    """Inline legend identifying PL1_8 DC tissue lineages by sapphire shade."""
    from matplotlib.patches import Patch
    handles, labels = [], []
    for lin in lineages:
        fcolor, label = PL18_LINEAGE_STYLE[lin]
        handles.append(Patch(facecolor=lighten(fcolor, 0.65), edgecolor=fcolor, linewidth=0.8))
        labels.append(label)
    ax.legend(handles, labels, title='DC lineage', loc='upper right', fontsize=7.5,
              title_fontsize=7.5, frameon=False, handlelength=1.0, handletextpad=0.4,
              labelspacing=0.2, borderaxespad=0.3)


def main():
    np.random.seed(42)
    output_value = os.environ.get('PCWDE_EXPRESSION_OUTPUT_DIR', '')
    if not output_value:
        raise SystemExit('PCWDE_EXPRESSION_OUTPUT_DIR is required')
    output_dir = Path(output_value)
    if output_dir.exists() and any(output_dir.iterdir()):
        raise SystemExit(f'refusing non-empty output directory: {output_dir}')
    output_dir.mkdir(parents=True, exist_ok=True)

    family_data = load_data()

    # Determine GH5_2 loci per species (respecting display order)
    gh5 = family_data['GH5_2']
    skip = {'uncertain', 'chr10_uncertain'}
    dc_loci_raw = set(gh5[gh5['_species'] == 'DC']['classification'].dropna().unique())
    bk_loci_raw = set(gh5[gh5['_species'] == 'BK']['classification'].dropna().unique())
    dc_loci = [l for l in GH5_LOCUS_ORDER if l in dc_loci_raw - skip]
    bk_loci = [l for l in GH5_LOCUS_ORDER if l in bk_loci_raw - skip]
    print(f"GH5_2 DC loci: {dc_loci}")
    print(f"GH5_2 BK loci: {bk_loci}")

    # Five gene rows, with the two species side by side in each row.
    import matplotlib.gridspec as gridspec
    fig = plt.figure(figsize=(6.8, 10.0))
    fig.patch.set_facecolor('white')
    outer = gridspec.GridSpec(5, 2, figure=fig,
                              hspace=1.2, wspace=0.20, bottom=0.055, top=0.94,
                              left=0.12, right=0.98)

    axes = {}  # (grid_row, grid_col) -> ax or (ax_top, ax_bot)

    # DC GH5_2 retains its broken axis.
    for grow in range(5):
        for gcol in range(2):
            if grow == 2 and gcol == 0:
                inner = gridspec.GridSpecFromSubplotSpec(
                    2, 1, subplot_spec=outer[grow, gcol],
                    height_ratios=[1, 3], hspace=0.08)
                ax_top = fig.add_subplot(inner[0])
                ax_bot = fig.add_subplot(inner[1])
                axes[(grow, gcol)] = (ax_top, ax_bot)
            else:
                axes[(grow, gcol)] = fig.add_subplot(outer[grow, gcol])

    ROWS = [
        ('DC',  DC_TISSUES, DC_ORDER, DC_LARVAL_N, 'D. cinerosa'),
        ('BK',  BK_TISSUES, BK_ORDER, BK_LARVAL_N, 'B. kinseyi'),
    ]

    # Top section: 3 PCWDE families (cols 0-2 in grid rows 0-1)
    TOP_COLUMNS = [
        ('PL1_8',  'PL1_8',  '#2166AC', 'simple'),
        ('PL4_4',  'PL4_4',  '#B2182B', 'simple'),
        ('GH5_2',  'GH5_2',  '#1B7837', 'gh5_grouped'),
    ]

    # PL1_8 lineages present per species (drops non-expressed paralogs)
    pl18 = family_data['PL1_8']
    pl18_dc_lineages = pl18_lineages_present(pl18, 'DC', DC_TISSUES)
    pl18_bk_lineages = pl18_lineages_present(pl18, 'BK', BK_TISSUES)
    print(f"PL1_8 DC lineages: {pl18_dc_lineages}")
    print(f"PL1_8 BK lineages: {pl18_bk_lineages}")
    # Final columns: GH1_e72 comparison + B-lectin domain
    BOT_COLUMNS = [
        ('nonHGT', NONHGT_LABEL, '#333333', 'nonhgt'),
        ('CBM13',  'B-lectin',  CBM13_COLOR, 'cbm13'),
    ]

    font_title = 11
    font_tick = 9
    font_ylabel = 10

    ALL_COLUMNS = TOP_COLUMNS + BOT_COLUMNS  # 5 cols
    for gcol, (fname, display_name, fcolor, draw_mode) in enumerate(ALL_COLUMNS):
        df = family_data[fname]
        for sp_i, (species, tissue_map, tissue_order, n_larval, sp_label) in enumerate(ROWS):
            grow = gcol
            cell = axes[(grow, sp_i)]

            if grow == 2 and sp_i == 0:
                ax_top, ax_bot = cell
                loci = dc_loci
                draw_grouped_panel(ax_bot, df, species, tissue_map, tissue_order,
                                   loci, inferred=True, label_ymax=BREAK_LO * 0.85)
                draw_grouped_panel(ax_top, df, species, tissue_map, tissue_order,
                                   loci, inferred=True, show_labels=False)
                style_ax(ax_bot)
                ax_bot.set_ylim(0, BREAK_LO)
                ax_bot.tick_params(axis='x', rotation=45, labelsize=font_tick)
                ax_bot.tick_params(axis='y', labelsize=font_tick)
                add_stage_brackets(ax_bot, n_larval, len(tissue_order),
                                   y_line=-0.83, y_text=-0.91)
                outlier_vals = []
                for tissue in tissue_order:
                    for locus in loci:
                        v = get_replicate_totals(df, species, tissue_map[tissue], locus)
                        for val in v:
                            if val > BREAK_LO:
                                outlier_vals.append(val)
                hi_lo = min(outlier_vals) * 0.85 if outlier_vals else BREAK_LO
                hi_hi = max(outlier_vals) * 1.25 if outlier_vals else BREAK_LO * 2
                ax_top.set_ylim(hi_lo, hi_hi)
                ax_top.set_facecolor('white')
                ax_top.spines['top'].set_visible(False)
                ax_top.spines['right'].set_visible(False)
                ax_top.spines['bottom'].set_visible(False)
                ax_top.tick_params(axis='x', bottom=False, labelbottom=False)
                ax_top.tick_params(axis='y', labelsize=font_tick)
                ax_bot.spines['top'].set_visible(False)
                ax_top.spines['bottom'].set_visible(False)
                add_break_marks(ax_top, ax_bot)
                if sp_i == 0:
                    ax_top.set_title(display_name, fontsize=font_title, fontweight='bold',
                                     color=fcolor, loc='left', pad=28)
            else:
                ax = cell
                if draw_mode == 'gh5_grouped':
                    loci = dc_loci if species == 'DC' else bk_loci
                    draw_grouped_panel(ax, df, species, tissue_map, tissue_order,
                                       loci, inferred=(species == 'DC'))
                elif draw_mode == 'pl18_lineage':
                    lins = pl18_dc_lineages if species == 'DC' else pl18_bk_lineages
                    draw_grouped_panel(ax, df, species, tissue_map, tissue_order,
                                       lins, inferred=(species == 'BK'), show_labels=False,
                                       style_map=PL18_LINEAGE_STYLE)
                    if sp_i == 0:
                        draw_pl18_legend(ax, pl18_dc_lineages)
                elif draw_mode == 'nonhgt':
                    draw_nonhgt_panel(ax, df, species, tissue_map, tissue_order)
                    if sp_i == 0:
                        draw_nonhgt_legend(ax)
                else:
                    draw_simple_panel(ax, df, species, tissue_map, tissue_order, fcolor)
                style_ax(ax)
                ax.tick_params(axis='x', rotation=45, labelsize=font_tick)
                ax.tick_params(axis='y', labelsize=font_tick)
                add_stage_brackets(ax, n_larval, len(tissue_order))
                if sp_i == 0:
                    ax.set_title(display_name, fontsize=font_title, fontweight='bold', color=fcolor, loc='left')

            ax_label = cell[1] if isinstance(cell, tuple) else cell
            if sp_i == 0:
                ax_label.set_ylabel('Total TPM', fontsize=font_ylabel)

    fig.text(0.32, 0.995, 'D. cinerosa', ha='center', va='top', fontsize=11,
             fontweight='bold', fontstyle='italic')
    fig.text(0.75, 0.995, 'B. kinseyi', ha='center', va='top', fontsize=11,
             fontweight='bold', fontstyle='italic')

    out_svg = output_dir / 'expression.svg'
    out_pdf = output_dir / 'expression.pdf'
    out_png = output_dir / 'expression.png'
    plt.savefig(out_svg, format='svg', bbox_inches='tight', facecolor='white')
    plt.savefig(out_pdf, bbox_inches='tight', facecolor='white')
    plt.savefig(out_png, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"Saved {out_svg}, {out_pdf.name}, {out_png.name}")


if __name__ == '__main__':
    main()
