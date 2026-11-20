#!/usr/bin/env python3
"""
05_filter_for_visualization.py - Filter to synteny-informative genes

Keeps TARGET genes and FLANKING genes in CynOGs shared across >= 3 genomes
within the same block. Singletons keep all flanking genes (no filter).

Inputs:
    outputs/region_gene_inventory.tsv    Full gene inventory
    outputs/region_status.tsv            For block_id lookup

Output:
    outputs/region_viz_subset.tsv        Filtered for visualization
"""

import pandas as pd
from pathlib import Path
from collections import defaultdict

OUTPUT_DIR = Path(__file__).parent / 'outputs'

INVENTORY_FILE = OUTPUT_DIR / 'region_gene_inventory.tsv'
STATUS_FILE = OUTPUT_DIR / 'region_status.tsv'

MIN_GENOMES_FOR_CYNOG = 3


def main():
    print("=" * 60)
    print("Step 05: Filter for Visualization (CynOG)")
    print("=" * 60)

    inventory_df = pd.read_csv(INVENTORY_FILE, sep='\t')
    status_df = pd.read_csv(STATUS_FILE, sep='\t')
    print(f"Loaded {len(inventory_df)} genes from inventory")

    # Build region_id -> block_id mapping
    region_block = dict(zip(status_df['region_id'], status_df['block_id']))

    targets = inventory_df[inventory_df['gene_type'] == 'TARGET'].copy()
    flanking = inventory_df[inventory_df['gene_type'] == 'FLANKING'].copy()
    print(f"  TARGET: {len(targets)}, FLANKING: {len(flanking)}")

    # Add block_id to flanking
    flanking = flanking.copy()
    flanking['block_id'] = flanking['region_id'].map(region_block)

    # For blocks: count genomes per CynOG within each block
    block_cynog_genomes = defaultdict(lambda: defaultdict(set))
    for _, gene in flanking[flanking['block_id'] != 'SINGLETON'].iterrows():
        bid = gene['block_id']
        cynog = gene['cynog_id']
        genome = gene['genome']
        if cynog and pd.notna(cynog) and cynog != '':
            block_cynog_genomes[bid][cynog].add(genome)

    # Shared CynOGs per block
    shared_cynogs_by_block = {}
    for bid, cynog_genomes in block_cynog_genomes.items():
        shared_cynogs_by_block[bid] = {
            cynog for cynog, genomes in cynog_genomes.items()
            if len(genomes) >= MIN_GENOMES_FOR_CYNOG
        }

    # Filter flanking genes
    filtered_flanking = []
    for _, gene in flanking.iterrows():
        bid = gene.get('block_id')
        cynog = gene['cynog_id']
        annotation = str(gene['annotation']).lower()

        if bid == 'SINGLETON' or pd.isna(bid):
            # Singletons: include all flanking
            row = gene.drop('block_id').to_dict()
            filtered_flanking.append(row)
        else:
            # Blocks: only include if CynOG is shared across >= N genomes
            shared_cynogs = shared_cynogs_by_block.get(bid, set())
            if cynog and pd.notna(cynog) and cynog in shared_cynogs:
                row = gene.drop('block_id').to_dict()
                # Replace hypothetical annotations with CynOG ID
                if 'hypothetical' in annotation or 'uncharacterized' in annotation:
                    row['annotation'] = cynog
                filtered_flanking.append(row)

    filtered_df = pd.DataFrame(filtered_flanking)
    if len(filtered_df) > 0:
        result_df = pd.concat([targets, filtered_df], ignore_index=True)
    else:
        result_df = targets.copy()

    result_df = result_df.sort_values(['region_id', 'gene_rank'])

    output_file = OUTPUT_DIR / 'region_viz_subset.tsv'
    result_df.to_csv(output_file, sep='\t', index=False)
    print(f"\nSaved {len(result_df)} genes to {output_file}")
    print(f"  TARGET: {(result_df['gene_type'] == 'TARGET').sum()}")
    print(f"  FLANKING: {(result_df['gene_type'] == 'FLANKING').sum()}")

    # Stats per block
    for bid in sorted(shared_cynogs_by_block.keys()):
        n_shared = len(shared_cynogs_by_block[bid])
        print(f"  {bid}: {n_shared} shared CynOGs (>={MIN_GENOMES_FOR_CYNOG} genomes)")


if __name__ == '__main__':
    main()
