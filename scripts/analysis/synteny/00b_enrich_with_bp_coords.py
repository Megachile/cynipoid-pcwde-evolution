#!/usr/bin/env python3
"""
00b_enrich_with_bp_coords.py — Add bp coordinates from GFF files to chromosome_cynog_dict.tsv

Uses pipeline_utils.parse_gff_gene_positions() to extract (chrom, start, end, strand)
for each protein from the Helixer/NCBI/evidence GFF files. Joins on protein_id to
add start_bp, end_bp, strand columns.

GFF selection follows the same priority as phase2_synteny.py:
  NCBI > evidence > Helixer ab initio

Input:
    chromosome_cynog_dict.tsv: genome, chrom, position, protein_id, cynog_id

Output:
    chromosome_cynog_dict.tsv (overwritten): + start_bp, end_bp, strand columns
"""

import sys
import os
from pathlib import Path

BASE = Path(os.environ['PCWDE_PROJECT_ROOT']).resolve()
SCRIPT_DIR = Path(__file__).parent
GFF_DIR = BASE / 'genomes/pancynipoid_helixer'

# Add pipeline_utils to path
sys.path.insert(0, str(BASE / 'analysis/cynog_catalog'))
from pipeline_utils import parse_gff_gene_positions

CHR_CYNOG_FILE = SCRIPT_DIR.parent / 'chromosome_cynog_dict.tsv'


def find_gff(genome_dir: Path, accession: str) -> Path:
    """Locate the canonical GFF file for a genome (same as phase2_synteny.py)."""
    candidates = [
        genome_dir / f"{accession}_ncbi.gff",
        genome_dir / f"{accession}_evidence.gff3",
        genome_dir / f"{accession}_helixer.gff3",
        genome_dir / f"{accession}.gff3",
        genome_dir / f"{accession}_helixer.gff",
    ]
    for c in candidates:
        if c.exists():
            return c
    raise FileNotFoundError(
        f"No GFF found for {accession} in {genome_dir}. "
        f"Searched: {[c.name for c in candidates]}"
    )


def main():
    import pandas as pd

    print("=" * 60)
    print("Step 00b: Enrich chromosome_cynog_dict with bp coordinates")
    print("=" * 60)

    # Load current dict
    df = pd.read_csv(CHR_CYNOG_FILE, sep='\t')
    print(f"Loaded {len(df):,} rows, {df['genome'].nunique()} genomes")

    # Build protein_id -> (start_bp, end_bp, strand) from GFFs
    bp_map = {}  # protein_id -> (start, end, strand)
    genomes = sorted(df['genome'].unique())
    n_found = 0
    n_missing_gff = 0

    for genome in genomes:
        genome_dir = GFF_DIR / genome
        if not genome_dir.exists():
            print(f"  WARNING: no dir for {genome}")
            n_missing_gff += 1
            continue

        try:
            gff_path = find_gff(genome_dir, genome)
        except FileNotFoundError as e:
            print(f"  WARNING: {e}")
            n_missing_gff += 1
            continue

        # Get all protein IDs for this genome
        genome_pids = set(df.loc[df['genome'] == genome, 'protein_id'])

        positions = parse_gff_gene_positions(gff_path, protein_ids=genome_pids)
        for entry in positions:
            bp_map[entry['protein_id']] = (entry['start'], entry['end'], entry['strand'])

        n_found += len(positions)
        print(f"  {genome}: {len(positions)} genes from {gff_path.name}")

    print(f"\nTotal: {n_found:,} genes with bp coords from {len(genomes) - n_missing_gff} genomes")

    # Join
    df['start_bp'] = df['protein_id'].map(lambda pid: bp_map.get(pid, (None, None, None))[0])
    df['end_bp'] = df['protein_id'].map(lambda pid: bp_map.get(pid, (None, None, None))[1])
    df['strand'] = df['protein_id'].map(lambda pid: bp_map.get(pid, (None, None, None))[2])

    n_with_bp = df['start_bp'].notna().sum()
    n_total = len(df)
    print(f"  {n_with_bp:,}/{n_total:,} genes have bp coords ({100*n_with_bp/n_total:.1f}%)")

    # Convert to int where available
    for col in ['start_bp', 'end_bp']:
        df[col] = df[col].astype('Int64')  # nullable int

    # Save
    df.to_csv(CHR_CYNOG_FILE, sep='\t', index=False)
    print(f"\nWrote enriched {CHR_CYNOG_FILE}")
    print(f"  Columns: {list(df.columns)}")


if __name__ == '__main__':
    main()
