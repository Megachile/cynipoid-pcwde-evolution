#!/usr/bin/env python3
"""
04_region_gene_inventory.py - Build full gene inventory for each region

For each region (block x genome), lists ALL genes with annotations, CynOG
assignments, and HMM hits. Uses gene-count windows (not physical distance)
for flanking context, with bp coordinates for display.

Inputs:
    outputs/region_status.tsv            Region definitions + status
    chromosome_cynog_dict.tsv            Gene positions + CynOG + bp coords
    mc_assignments.domtblout             HMM hits for annotation
    Helixer protein FASTA files          Annotations from headers

Output:
    outputs/region_gene_inventory.tsv    Comprehensive gene listing per region
"""

import pandas as pd
import os
from pathlib import Path
from collections import defaultdict

BASE = Path(os.environ['PCWDE_PROJECT_ROOT']).resolve()
SCRIPT_DIR = Path(__file__).parent
OUTPUT_DIR = SCRIPT_DIR / 'outputs'
GFF_DIR = BASE / 'genomes/pancynipoid_helixer'

# Input files
STATUS_FILE = OUTPUT_DIR / 'region_status.tsv'
CHR_CYNOG_FILE = SCRIPT_DIR.parent / 'chromosome_cynog_dict.tsv'
DOMTBL_FILE = BASE / 'analysis/locus_catalog/phase2_locus_catalog/static_inputs/mc_assignments.domtblout'

# Parameters
N_FLANKING_GENES = 20  # genes per side (gene-count window)
EVALUE_THRESHOLD = 1e-5
TARGET_PREFIXES = ('MC_0006', 'MC_0030', 'MC_0047')

# Target CynOG families (for labeling)
CYNOG_TO_FAMILY = {
    'CynOGv3_000234': 'GH5',
    'CynOGv3_000148': 'PL1',
    'CynOGv3_000168': 'PL4',
    'CynOGv3_010310': 'CBM',
}


def load_chr_cynog_for_genomes(genomes):
    """Load chromosome_cynog_dict.tsv for specified genomes.

    Uses gene-count rank (position column) for ordering.
    Includes bp coordinates (start_bp, end_bp, strand) for display.

    Returns:
        genes_by_scaffold: {(genome, scaffold): [(rank, start_bp, end_bp, strand, protein_id, cynog_id), ...]}
                           sorted by rank
        pid_to_loc: {protein_id: ((genome, scaffold), list_index)}
    """
    print("Loading chromosome_cynog_dict.tsv...")
    result = defaultdict(list)

    with open(CHR_CYNOG_FILE) as f:
        header = f.readline().rstrip('\n').split('\t')
        col_idx = {name: i for i, name in enumerate(header)}

        for line in f:
            parts = line.rstrip('\n').split('\t')
            genome = parts[col_idx['genome']]
            if genome not in genomes:
                continue
            scaffold = parts[col_idx['chrom']]
            rank = int(parts[col_idx['position']])
            protein_id = parts[col_idx['protein_id']]
            cynog_id = parts[col_idx['cynog_id']] if parts[col_idx['cynog_id']] else ''

            # bp coords (may be empty)
            start_bp = parts[col_idx['start_bp']] if 'start_bp' in col_idx and parts[col_idx['start_bp']] else ''
            end_bp = parts[col_idx['end_bp']] if 'end_bp' in col_idx and parts[col_idx['end_bp']] else ''
            strand = parts[col_idx['strand']] if 'strand' in col_idx and parts[col_idx['strand']] else ''

            result[(genome, scaffold)].append((rank, start_bp, end_bp, strand, protein_id, cynog_id))

    # Sort each scaffold by rank
    for key in result:
        result[key].sort(key=lambda x: x[0])

    # Build protein_id -> (key, index) lookup
    pid_to_loc = {}
    for key, genes in result.items():
        for idx, (rank, sbp, ebp, strd, pid, cid) in enumerate(genes):
            pid_to_loc[pid] = (key, idx)

    n = sum(len(v) for v in result.values())
    print(f"  {n:,} genes for {len(genomes)} genomes")
    return result, pid_to_loc


def load_domtblout():
    """Load HMM hits for PCWDE families. Returns protein_id -> best hit."""
    print("Loading HMM hits...")
    hits = {}

    if not DOMTBL_FILE.exists():
        print(f"  WARNING: {DOMTBL_FILE} not found, skipping HMM hits")
        return hits

    with open(DOMTBL_FILE) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 13:
                continue

            mc = parts[0]
            protein_id = parts[3]
            evalue = float(parts[6])

            if not mc.startswith(TARGET_PREFIXES):
                continue

            # Normalize protein ID
            if protein_id.endswith('.1'):
                protein_id = protein_id[:-2]

            if evalue <= EVALUE_THRESHOLD:
                if protein_id not in hits or evalue < hits[protein_id][1]:
                    hits[protein_id] = (mc, evalue)

    print(f"  {len(hits)} proteins with PCWDE hits")
    return hits


CYNOG_NR_ANNOT_FILE = BASE / 'good_data/catalogs/cynog/catalog/cynog_nr_annotations.tsv'

def load_cynog_annotations():
    """Load CynOGv3 -> consensus nr annotation from the canonical family
    catalog lookup. Strip trailing '[species]' brackets to get clean text.
    """
    print("Loading CynOG consensus annotations...")
    if not CYNOG_NR_ANNOT_FILE.exists():
        print(f"  WARNING: {CYNOG_NR_ANNOT_FILE} not found")
        return {}
    df = pd.read_csv(CYNOG_NR_ANNOT_FILE, sep='\t')
    out = {}
    import re as _re
    for _, r in df.iterrows():
        text = r.get('consensus_nr_annotation', '')
        if pd.isna(text):
            continue
        # Strip trailing "[Species name]" bracket if present
        text = _re.sub(r'\s*\[[^\]]+\]\s*$', '', str(text)).strip()
        out[r['cynog_id']] = text
    print(f"  {len(out)} CynOG annotations")
    return out

def load_annotations(genomes):
    """Back-compat wrapper: builds a protein_id -> annotation dict by joining
    chromosome_cynog_dict's protein_id -> cynog_id against the CynOG-level
    consensus annotation. Replaces the old FASTA-header scraper, which
    returned empty strings on the v3 proteomes (bare IDs, no descriptions).
    """
    cynog_annot = load_cynog_annotations()
    print("Joining CynOG annotations to proteins via chromosome_cynog_dict...")
    chd = pd.read_csv(CHR_CYNOG_FILE, sep='\t')
    chd_sub = chd[chd['genome'].isin(genomes)]
    annotations = {}
    for pid, cynog in zip(chd_sub['protein_id'], chd_sub['cynog_id']):
        if pd.notna(cynog) and cynog in cynog_annot:
            annotations[pid] = cynog_annot[cynog]
    print(f"  {len(annotations)} annotated genes")
    return annotations


def main():
    print("=" * 60)
    print("Step 04: Region Gene Inventory (CynOG)")
    print(f"  Gene-count window: {N_FLANKING_GENES} genes per side")
    print("=" * 60)

    status_df = pd.read_csv(STATUS_FILE, sep='\t')
    print(f"Loaded {len(status_df)} regions")

    # Filter to regions that are PRESENT or EMPTY (not N_D)
    valid = status_df[status_df['status'] != 'N_D'].copy()
    print(f"  {len(valid)} regions with status PRESENT or EMPTY")

    genomes = set(valid['genome'])
    genes_by_scaffold, pid_to_loc = load_chr_cynog_for_genomes(genomes)
    hmm_hits = load_domtblout()
    annotations = load_annotations(genomes)

    all_genes = []
    n_present_no_targets = 0
    n_regions_processed = 0

    for idx, row in valid.iterrows():
        region_id = row['region_id']
        family = row['family']
        genome = row['genome']
        scaffold = row['scaffold']
        status = row['status']

        genes = genes_by_scaffold.get((genome, scaffold), [])
        if not genes:
            continue

        # Get target gene IDs from region_status
        target_ids = set()
        if pd.notna(row.get('target_gene_ids')) and row['target_gene_ids']:
            target_ids = set(str(row['target_gene_ids']).split(','))

        # Determine anchor indices using gene-count rank
        if status in ('PRESENT', 'UNCONFIRMED') and target_ids:
            target_indices = []
            for tid in target_ids:
                if tid in pid_to_loc:
                    loc_key, loc_idx = pid_to_loc[tid]
                    if loc_key == (genome, scaffold):
                        target_indices.append(loc_idx)

            if not target_indices:
                n_present_no_targets += 1
                continue

            min_target_idx = min(target_indices)
            max_target_idx = max(target_indices)
        elif status == 'EMPTY':
            # Use the cluster bounds that 03 wrote via flanking CynOG search.
            # These mark where the conserved anchors actually cluster on this
            # scaffold, which is the region 03 used to call this row EMPTY.
            region_start = row.get('region_start_rank', '')
            region_end   = row.get('region_end_rank', '')
            if pd.notna(region_start) and str(region_start).strip() != '':
                cluster_start = int(float(region_start))
                cluster_end   = int(float(region_end))
                cluster_mid   = (cluster_start + cluster_end) // 2
                # Find index of gene nearest cluster_mid by rank
                best_i, best_d = 0, abs(genes[0][0] - cluster_mid)
                for gi in range(1, len(genes)):
                    d = abs(genes[gi][0] - cluster_mid)
                    if d < best_d:
                        best_i, best_d = gi, d
                min_target_idx = best_i
                max_target_idx = best_i
            else:
                # Fallback: scaffold midpoint (should only happen for very
                # old region_status files that lack the columns).
                min_target_idx = len(genes) // 2
                max_target_idx = len(genes) // 2
        else:
            continue

        # Gene-count window: walk N_FLANKING_GENES in each direction
        start_idx = max(0, min_target_idx - N_FLANKING_GENES)
        end_idx = min(len(genes) - 1, max_target_idx + N_FLANKING_GENES)

        # Compute region center in bp (for relative positioning display)
        target_start_bp = genes[min_target_idx][1]  # start_bp
        target_end_bp = genes[max_target_idx][2]    # end_bp
        has_bp = target_start_bp and target_end_bp and target_start_bp != ''

        if has_bp:
            try:
                region_center = (int(target_start_bp) + int(target_end_bp)) / 2
            except (ValueError, TypeError):
                has_bp = False
                region_center = 0
        else:
            region_center = 0

        for i in range(start_idx, end_idx + 1):
            rank, start_bp, end_bp, strand, protein_id, cynog_id = genes[i]

            # Relative position in gene-counts from target center
            target_center_rank = (genes[min_target_idx][0] + genes[max_target_idx][0]) / 2
            relative_rank = rank - target_center_rank

            # Relative position in kb (for display)
            relative_pos_kb = None
            if has_bp and start_bp and end_bp and start_bp != '':
                try:
                    gene_center = (int(start_bp) + int(end_bp)) / 2
                    relative_pos_kb = round((gene_center - region_center) / 1000, 2)
                except (ValueError, TypeError):
                    pass

            gene_type = 'TARGET' if protein_id in target_ids else 'FLANKING'
            annotation = annotations.get(protein_id, '')

            hit_mc, hit_evalue = '', None
            if protein_id in hmm_hits:
                hit_mc, hit_evalue = hmm_hits[protein_id]

            all_genes.append({
                'region_id': region_id,
                'genome': genome,
                'family': family,
                'scaffold': scaffold,
                'gene_id': protein_id,
                'gene_rank': rank,
                'gene_start': int(start_bp) if start_bp and start_bp != '' else None,
                'gene_end': int(end_bp) if end_bp and end_bp != '' else None,
                'strand': strand if strand else None,
                'relative_rank': relative_rank,
                'relative_pos_kb': relative_pos_kb,
                'gene_type': gene_type,
                'annotation': annotation,
                'cynog_id': cynog_id,
                'hit_mc': hit_mc,
                'hit_evalue': hit_evalue,
            })

        n_regions_processed += 1

    result_df = pd.DataFrame(all_genes)
    result_df = result_df.sort_values(['region_id', 'gene_rank'])

    output_file = OUTPUT_DIR / 'region_gene_inventory.tsv'
    result_df.to_csv(output_file, sep='\t', index=False)
    print(f"\nSaved {len(result_df)} genes to {output_file}")

    print(f"\n  Regions processed: {n_regions_processed}")
    print(f"  TARGET genes: {(result_df['gene_type'] == 'TARGET').sum()}")
    print(f"  FLANKING genes: {(result_df['gene_type'] == 'FLANKING').sum()}")
    print(f"  Genes with CynOG: {(result_df['cynog_id'] != '').sum()}")
    print(f"  Genes with bp coords: {result_df['gene_start'].notna().sum()}")
    print(f"  Genes with PCWDE hits: {(result_df['hit_mc'] != '').sum()}")
    if n_present_no_targets > 0:
        print(f"  WARNING: {n_present_no_targets} PRESENT regions had no target genes in chr_cynog_dict")


if __name__ == '__main__':
    main()
