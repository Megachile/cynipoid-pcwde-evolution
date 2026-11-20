#!/usr/bin/env python3
"""
06_region_spreadsheets.py - Per-block spreadsheets for visualization

Creates ONE spreadsheet per block (all genomes as rows), showing flanking
genes on each side of the target locus. Genomes sorted by phylogenetic order.

Uses gene-count windows for flanking context, bp coordinates for display columns.

Inputs:
    outputs/region_viz_subset.tsv    Filtered genes
    outputs/region_status.tsv        Region status + block_id

Output:
    spreadsheet_{family}_{chr#letter}.tsv    Per-block visual layout
"""

import pandas as pd
import re
import os
from pathlib import Path
from collections import defaultdict

BASE = Path(os.environ['PCWDE_PROJECT_ROOT']).resolve()
SCRIPT_DIR = Path(__file__).parent
OUTPUT_DIR = SCRIPT_DIR / 'outputs'
PHYLO_FILE = BASE / 'genomes/pancynipoid_helixer/gca_to_species.tsv'

VIZ_FILE = OUTPUT_DIR / 'region_viz_subset.tsv'
STATUS_FILE = OUTPUT_DIR / 'region_status.tsv'

N_FLANKING = 12
EVALUE_THRESHOLD = 1e-10

CYNOG_TO_CAZY = {
    'CynOGv3_000234': 'GH5',
    'CynOGv3_000148': 'PL1',
    'CynOGv3_000168': 'PL4',
    'CynOGv3_010310': 'CBM',
}

CHR_TO_SCAFFOLD = {
    'chr1': 'CM021338.1', 'chr2': 'CM021339.1', 'chr3': 'CM021340.1',
    'chr4': 'CM021341.1', 'chr5': 'CM021342.1', 'chr6': 'CM021343.1',
    'chr7': 'CM021344.1', 'chr8': 'CM021345.1', 'chr9': 'CM021346.1',
    'chr10': 'CM021347.1',
}
SCAFFOLD_TO_CHR = {}
for chrom, base_acc in CHR_TO_SCAFFOLD.items():
    SCAFFOLD_TO_CHR[base_acc] = chrom
    SCAFFOLD_TO_CHR[base_acc + '_RagTag'] = chrom


PROTEIN_TO_CHR_FILE = BASE / 'good_data/catalogs/cynog/catalog/protein_to_chromosome.tsv'
_PTC_SCAFFOLD_TO_CHR = None

def _load_ptc_lookup():
    global _PTC_SCAFFOLD_TO_CHR
    print(f"  Loading scaffold->bk_chr from {PROTEIN_TO_CHR_FILE.name}...")
    df = pd.read_csv(PROTEIN_TO_CHR_FILE, sep='	', dtype={'bk_chr': str})
    df = df[df['bk_chr'].notna() & (df['bk_chr'] != 'unplaced')]
    modes = df.groupby('scaffold')['bk_chr'].agg(
        lambda s: s.mode().iloc[0] if len(s.mode()) else None
    )
    table = {}
    for scaffold, bk in modes.items():
        if bk:
            chrom = f'chr{bk}'
            table[scaffold] = chrom
            if scaffold.endswith('_RagTag'):
                table[scaffold[:-len('_RagTag')]] = chrom
    _PTC_SCAFFOLD_TO_CHR = table

def get_chromosome(scaffold):
    global _PTC_SCAFFOLD_TO_CHR
    if _PTC_SCAFFOLD_TO_CHR is None:
        _load_ptc_lookup()
    if not scaffold or pd.isna(scaffold):
        return '?'
    s = str(scaffold)
    if s in _PTC_SCAFFOLD_TO_CHR:
        return _PTC_SCAFFOLD_TO_CHR[s]
    base = s.replace('_RagTag', '')
    if base in _PTC_SCAFFOLD_TO_CHR:
        return _PTC_SCAFFOLD_TO_CHR[base]
    for acc, chrom in SCAFFOLD_TO_CHR.items():
        if acc in s:
            return chrom
    return 'frag'


def chr_sort_key(c):
    if c.startswith('chr'):
        try:
            return (0, int(c[3:]))
        except ValueError:
            return (0, 99)
    elif c == 'frag':
        return (1, 0)
    return (2, 0)


def clean_annotation(annot, max_len=40):
    if not annot or pd.isna(annot):
        return ''
    annot = str(annot)
    if 'hypothetical' in annot.lower() or 'uncharacterized' in annot.lower():
        return '[hyp]'
    if len(annot) > max_len:
        return annot[:max_len - 3] + '...'
    return annot


def _is_cynog_fallback(annot):
    return bool(annot and re.match(r'^CynOGv3_\d+$', str(annot)))


def _select_flanking(genes_df, n_slots):
    """Select best flanking genes, preferring real annotations over CynOG-only."""
    real = []
    fallback = []
    for _, gene in genes_df.iterrows():
        annot = gene['annotation']
        if _is_cynog_fallback(annot):
            fallback.append((gene['gene_rank'], annot))
        else:
            real.append((gene['gene_rank'], annot))

    selected = real[:n_slots]
    remaining = n_slots - len(selected)
    if remaining > 0:
        selected.extend(fallback[:remaining])

    selected.sort(key=lambda x: x[0], reverse=(
        genes_df['gene_rank'].iloc[0] > genes_df['gene_rank'].iloc[-1]
        if len(genes_df) > 1 else False
    ))
    return [annot for _, annot in selected]


def build_block_chrletter_map(status_df):
    """Build block_id -> chr#letter mapping for consistent naming."""
    block_to_chrletter = {}

    for family in status_df['family'].unique():
        fam_df = status_df[status_df['family'] == family]
        blocks_fam = fam_df[fam_df['region_type'] == 'block']

        # Get unique block IDs and their chromosomes
        block_info = {}
        for _, row in blocks_fam.iterrows():
            bid = row['block_id']
            if bid not in block_info and bid != 'SINGLETON':
                present_rows = blocks_fam[(blocks_fam['block_id'] == bid) &
                                          (blocks_fam['status'] == 'PRESENT')]
                if len(present_rows) > 0:
                    chrom = get_chromosome(present_rows.iloc[0]['scaffold'])
                else:
                    chrom = '?'
                block_info[bid] = chrom

        # Group by chromosome, assign letters
        chr_blocks = defaultdict(list)
        for bid, chrom in block_info.items():
            chr_blocks[chrom].append(bid)

        for chrom in sorted(chr_blocks.keys(), key=chr_sort_key):
            bids = sorted(chr_blocks[chrom])
            for i, bid in enumerate(bids):
                if i < 26:
                    letter = chr(ord('a') + i)
                else:
                    letter = chr(ord('a') + i // 26 - 1) + chr(ord('a') + i % 26)
                block_to_chrletter[bid] = f"{chrom}{letter}"

    return block_to_chrletter


def main():
    print("=" * 60)
    print("Step 06: Region Spreadsheets (CynOG)")
    print("=" * 60)

    viz_df = pd.read_csv(VIZ_FILE, sep='\t')
    status_df = pd.read_csv(STATUS_FILE, sep='\t')
    print(f"Loaded {len(viz_df)} genes, {len(status_df)} regions")

    # Load phylogenetic order
    phylo_df = pd.read_csv(PHYLO_FILE, sep='\t')
    phylo_map = {}
    species_map = {}
    for _, row in phylo_df.iterrows():
        if row.get('status') == 'active':
            acc = row['accession']
            phylo_map[acc] = row.get('tree_order', 999)
            species_map[acc] = row.get('species', acc)

    block_to_chrletter = build_block_chrletter_map(status_df)
    print(f"Built chr#letter mapping for {len(block_to_chrletter)} blocks")

    # Build present_chrs lookup (include UNCONFIRMED — they have proteins at the locus)
    present_chrs_lookup = {}
    for _, row in status_df.iterrows():
        if row['status'] in ('PRESENT', 'UNCONFIRMED'):
            key = (row['family'], row['genome'])
            chrom = get_chromosome(row['scaffold'])
            if key not in present_chrs_lookup:
                present_chrs_lookup[key] = set()
            present_chrs_lookup[key].add(chrom)

    blocks_status = status_df[
        (status_df['region_type'] == 'block') &
        (status_df['block_id'] != 'SINGLETON')
    ].copy()

    n_spreadsheets = 0

    for block_id in sorted(blocks_status['block_id'].unique()):
        block_regions = blocks_status[blocks_status['block_id'] == block_id]
        family = block_regions['family'].iloc[0]

        rows = []
        for _, status_row in block_regions.iterrows():
            region_id = status_row['region_id']
            genome = status_row['genome']
            status = status_row['status']
            scaffold = status_row.get('scaffold', '')

            phylo_order = phylo_map.get(genome, 999)
            species = species_map.get(genome, genome)

            region_genes = viz_df[viz_df['region_id'] == region_id].sort_values('gene_rank')
            targets_all = region_genes[region_genes['gene_type'] == 'TARGET']
            flanking = region_genes[region_genes['gene_type'] == 'FLANKING']

            # Filter weak hits
            targets = targets_all[
                targets_all['hit_evalue'].apply(
                    lambda x: pd.isna(x) or float(x) <= EVALUE_THRESHOLD
                )
            ]
            actual_n_tandem = len(targets)

            # Downgrade EMPTY -> N_D when we do not have enough shared
            # flanking context to claim we recognized the region. The block
            # detection in 03 requires MIN_FLANKING_CYNOGS=3 conserved
            # matches to define a region; if 05 leaves us with fewer than
            # that after the conserved filter, we did not locate it well
            # enough and the honest semantic is N_D, not EMPTY.
            MIN_FLANKING_FOR_EMPTY = 3
            if status == 'EMPTY' and len(targets) == 0 and len(flanking) < MIN_FLANKING_FOR_EMPTY:
                status = 'N_D' 

            # bp display columns
            target_start_bp = ''
            target_end_bp = ''
            if len(targets) > 0 and 'gene_start' in targets.columns:
                valid_starts = targets['gene_start'].dropna()
                valid_ends = targets['gene_end'].dropna()
                if len(valid_starts) > 0:
                    target_start_bp = int(valid_starts.min())
                    target_end_bp = int(valid_ends.max())

            row_data = {
                'genome': genome,
                'species': species,
                'status': status,
                'n_tandem': actual_n_tandem,
                'phylo_order': phylo_order,
                'scaffold': scaffold,
                'start_bp': target_start_bp,
                'end_bp': target_end_bp,
            }

            if len(targets) > 0:
                center_rank = targets['gene_rank'].mean()
            else:
                # For empty regions use midpoint of flanking range
                if len(flanking) > 0:
                    center_rank = flanking['gene_rank'].mean()
                else:
                    center_rank = 0

            # Upstream flanking (lower rank)
            upstream = flanking[flanking['gene_rank'] < center_rank].sort_values(
                'gene_rank', ascending=False)
            upstream_selected = _select_flanking(upstream, N_FLANKING) if len(upstream) > 0 else []
            for i in range(N_FLANKING):
                col = f'U{N_FLANKING - i}'
                row_data[col] = clean_annotation(upstream_selected[i]) if i < len(upstream_selected) else ''

            # LOCUS column
            cazy_code = CYNOG_TO_CAZY.get(family, family)
            if status == 'UNCONFIRMED':
                locus_text = f"[unc:{len(targets_all)}]" if len(targets_all) > 0 else '[unc]'
            elif len(targets) > 0:
                locus_text = f"[{len(targets)}x] {cazy_code}" if len(targets) > 1 else cazy_code
            elif status == 'EMPTY':
                locus_text = '[empty]'
            elif len(targets_all) > 0:
                locus_text = 'Pseudogene'
            else:
                locus_text = '[n/d]'
            row_data['LOCUS'] = locus_text

            # Downstream flanking (higher rank)
            downstream = flanking[flanking['gene_rank'] > center_rank].sort_values('gene_rank')
            downstream_selected = _select_flanking(downstream, N_FLANKING) if len(downstream) > 0 else []
            for i in range(N_FLANKING):
                col = f'D{i + 1}'
                row_data[col] = clean_annotation(downstream_selected[i]) if i < len(downstream_selected) else ''

            # present_chrs
            key = (family, genome)
            if key in present_chrs_lookup:
                sorted_chrs = sorted(present_chrs_lookup[key], key=chr_sort_key)
                row_data['present_chrs'] = ','.join(sorted_chrs)
            else:
                row_data['present_chrs'] = ''

            rows.append(row_data)

        if not rows:
            continue

        df = pd.DataFrame(rows).sort_values('phylo_order')

        flanking_cols = ([f'U{i}' for i in range(N_FLANKING, 0, -1)]
                         + ['LOCUS']
                         + [f'D{i}' for i in range(1, N_FLANKING + 1)])
        col_order = ['genome', 'species', 'status', 'n_tandem',
                     'start_bp', 'end_bp',
                     'present_chrs', 'phylo_order', 'scaffold'] + flanking_cols
        df = df[[c for c in col_order if c in df.columns]]

        chrletter = block_to_chrletter.get(block_id, block_id)
        output_file = SCRIPT_DIR / f"spreadsheet_{family}_{chrletter}.tsv"
        df.to_csv(output_file, sep='\t', index=False)
        n_spreadsheets += 1
        print(f"  {output_file.name}: {len(df)} genomes")

    print(f"\nCreated {n_spreadsheets} block spreadsheets")


if __name__ == '__main__':
    main()
