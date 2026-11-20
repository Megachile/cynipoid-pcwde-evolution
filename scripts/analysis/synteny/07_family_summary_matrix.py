#!/usr/bin/env python3
"""
07_family_summary_matrix.py - Pan-family overview showing all loci (CynOG version)

Creates per-family summary matrix with columns named by chromosome + letter
(chr6a, chr6b, chr8a, ...). Rows are genomes in phylogenetic order.

Inputs:
    outputs/region_status.tsv    Region status for all regions

Output:
    summary_{CynOG_ID}.tsv       Per-family genome x region matrix
"""

import pandas as pd
import os
from pathlib import Path
from collections import defaultdict

BASE = Path(os.environ['PCWDE_PROJECT_ROOT']).resolve()
SCRIPT_DIR = Path(__file__).parent
OUTPUT_DIR = SCRIPT_DIR / 'outputs'
PHYLO_FILE = BASE / 'genomes/pancynipoid_helixer/gca_to_species.tsv'

STATUS_FILE = OUTPUT_DIR / 'region_status.tsv'

CHR_TO_SCAFFOLD = {
    'chr1': 'CM021338.1', 'chr2': 'CM021339.1', 'chr3': 'CM021340.1',
    'chr4': 'CM021341.1', 'chr5': 'CM021342.1', 'chr6': 'CM021343.1',
    'chr7': 'CM021344.1', 'chr8': 'CM021345.1', 'chr9': 'CM021346.1',
    'chr10': 'CM021347.1',
}
SCAFFOLD_TO_CHR = {}
for chrom, base in CHR_TO_SCAFFOLD.items():
    SCAFFOLD_TO_CHR[base] = chrom
    SCAFFOLD_TO_CHR[base + '_RagTag'] = chrom


def get_chromosome(scaffold):
    if not scaffold or pd.isna(scaffold):
        return '?'
    for acc, chrom in SCAFFOLD_TO_CHR.items():
        if acc in str(scaffold):
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


def main():
    print("=" * 60)
    print("Step 07: Family Summary Matrix (CynOG)")
    print("=" * 60)

    status_df = pd.read_csv(STATUS_FILE, sep='\t')
    print(f"Loaded {len(status_df)} region statuses")

    phylo_df = pd.read_csv(PHYLO_FILE, sep='\t')
    phylo_map = {}
    species_map = {}
    for _, row in phylo_df.iterrows():
        if row.get('status') == 'active':
            acc = row['accession']
            phylo_map[acc] = row.get('tree_order', 999)
            species_map[acc] = row.get('species', acc)

    families = sorted(status_df['family'].unique())

    for family in families:
        print(f"\nProcessing {family}...")
        fam_df = status_df[status_df['family'] == family].copy()

        blocks_df = fam_df[(fam_df['region_type'] == 'block') &
                           (fam_df['block_id'] != 'SINGLETON')]
        singletons_df = fam_df[fam_df['region_type'] == 'singleton']

        block_info = {}
        for bid in blocks_df['block_id'].unique():
            bdf = blocks_df[blocks_df['block_id'] == bid]
            present = bdf[bdf['status'] == 'PRESENT']
            if len(present) > 0:
                chrom = get_chromosome(present.iloc[0]['scaffold'])
            elif len(bdf) > 0:
                chrom = get_chromosome(bdf.iloc[0]['scaffold'])
            else:
                chrom = '?'
            block_info[bid] = chrom

        singleton_info = {}
        for bid in singletons_df['block_id'].unique():
            bdf = singletons_df[singletons_df['block_id'] == bid]
            owner_row = bdf[bdf['status'] == 'PRESENT'].iloc[0] if len(bdf[bdf['status'] == 'PRESENT']) > 0 else bdf.iloc[0]
            has_cross_genome = len(bdf) > 1
            singleton_info[bid] = {
                'genome': owner_row['genome'],
                'chrom': get_chromosome(owner_row['scaffold']),
                'n_tandem': owner_row.get('n_tandem', 1),
                'has_cross_genome': has_cross_genome,
            }

        frag_singletons = {}
        chr_loci = defaultdict(list)
        for bid, chrom in block_info.items():
            chr_loci[chrom].append(('block', bid, None))
        for sid, info in singleton_info.items():
            if info['chrom'] == 'frag':
                frag_singletons[sid] = info
            else:
                chr_loci[info['chrom']].append(('singleton', sid, info))

        locus_columns = []
        for chrom in sorted(chr_loci.keys(), key=chr_sort_key):
            loci = chr_loci[chrom]
            loci.sort(key=lambda x: (0 if x[0] == 'block' else 1, x[1]))
            for i, (ltype, lid, info) in enumerate(loci):
                if i < 26:
                    letter = chr(ord('a') + i)
                else:
                    letter = chr(ord('a') + i // 26 - 1) + chr(ord('a') + i % 26)
                col_name = f"{chrom}{letter}"
                locus_columns.append((col_name, ltype, lid, info))

        has_frag = len(frag_singletons) > 0

        all_genomes = sorted(fam_df['genome'].unique())
        rows = []

        for genome in all_genomes:
            phylo_order = phylo_map.get(genome, 999)
            species = species_map.get(genome, genome)
            present_chrs = []
            total_proteins = 0

            row = {
                'genome': genome,
                'species': species,
                'phylo_order': phylo_order,
            }

            for col_name, ltype, lid, info in locus_columns:
                if ltype == 'block':
                    region_id = f"{lid}__{genome}"
                    entry = blocks_df[blocks_df['region_id'] == region_id]

                    if len(entry) > 0:
                        e = entry.iloc[0]
                        if e['status'] == 'PRESENT':
                            count = int(e['n_tandem']) if pd.notna(e['n_tandem']) else 1
                            cell = str(count)
                            present_chrs.append(get_chromosome(e['scaffold']))
                            total_proteins += count
                        elif e['status'] == 'UNCONFIRMED':
                            count = int(e['n_tandem']) if pd.notna(e['n_tandem']) else 1
                            cell = f'unc({count})'
                        elif e['status'] == 'EMPTY':
                            cell = 'E'
                        else:
                            cell = 'n/d'
                    else:
                        cell = '-'
                else:
                    if info and info.get('has_cross_genome'):
                        entry = singletons_df[(singletons_df['block_id'] == lid) &
                                              (singletons_df['genome'] == genome)]
                        if len(entry) > 0:
                            e = entry.iloc[0]
                            if e['status'] == 'PRESENT':
                                count = int(e['n_tandem']) if pd.notna(e['n_tandem']) else 1
                                cell = str(count)
                                present_chrs.append(get_chromosome(e['scaffold']))
                                total_proteins += count
                            elif e['status'] == 'TIER2_ONLY':
                                count = int(e['n_tandem']) if pd.notna(e['n_tandem']) else 1
                                cell = f't2({count})'
                            elif e['status'] == 'EMPTY':
                                cell = 'E'
                            else:
                                cell = 'n/d'
                        else:
                            cell = '-'
                    elif info and info['genome'] == genome:
                        count = int(info.get('n_tandem', 1)) if info.get('n_tandem') else 1
                        cell = str(count)
                        present_chrs.append(info['chrom'])
                        total_proteins += count
                    else:
                        cell = '-'

                row[col_name] = cell

            if has_frag:
                frag_count = 0
                for sid, info in frag_singletons.items():
                    if info['genome'] == genome:
                        count = int(info.get('n_tandem', 1)) if info.get('n_tandem') else 1
                        frag_count += count
                if frag_count > 0:
                    row['frag'] = str(frag_count)
                    present_chrs.append('frag')
                    total_proteins += frag_count
                else:
                    row['frag'] = '-'

            row['n_loci'] = len(present_chrs)
            row['n_proteins'] = total_proteins
            row['present_chrs'] = ','.join(sorted(set(present_chrs))) if present_chrs else ''
            rows.append(row)

        result_df = pd.DataFrame(rows).sort_values('phylo_order')

        locus_col_names = [c[0] for c in locus_columns]
        if has_frag:
            locus_col_names.append('frag')
        col_order = (['genome', 'species', 'phylo_order', 'n_loci',
                       'n_proteins', 'present_chrs'] + locus_col_names)
        for col in col_order:
            if col not in result_df.columns:
                result_df[col] = ''
        result_df = result_df[[c for c in col_order if c in result_df.columns]]

        output_file = SCRIPT_DIR / f"summary_{family}.tsv"
        result_df.to_csv(output_file, sep='\t', index=False)
        print(f"  Saved {len(result_df)} genomes to {output_file.name}")

        n_blocks = len(block_info)
        n_chr_sing = len(singleton_info) - len(frag_singletons)
        n_frag_sing = len(frag_singletons)
        print(f"  Loci: {n_blocks} blocks + {n_chr_sing} chr singletons + {n_frag_sing} frag singletons (1 col) = {len(locus_col_names)} columns")

    print("\nDone!")


if __name__ == '__main__':
    main()
