#!/usr/bin/env python3
"""
03_detect_empty_loci.py - Detect syntenic loci where the target gene is absent

For each synteny block, checks ALL genomes that have the block's chromosome:
- Member genomes: confirmed PRESENT with coordinates
- Non-member genomes: scanned for family-specific CynOG members in the synteny region
  -> PRESENT if hits found, EMPTY if absent, N_D if chromosome/region not detectable

Adapted for CynOG system: uses gene-rank positions and membership_table
instead of bp coordinates and HMM domtblout.

Inputs:
    outputs/synteny_blocks.tsv     Block/singleton membership
    outputs/pcwde_loci.tsv         Flanking CynOG signatures per locus
    chromosome_cynog_dict.tsv      All gene positions + CynOG assignments
    membership_table.tsv           Protein -> CynOG mapping
    gca_to_species.tsv             Genome registry

Output:
    outputs/region_status.tsv      Status (PRESENT/EMPTY/N_D) per block x genome
"""

import pandas as pd
import os
from pathlib import Path
from collections import defaultdict, Counter
import sys

BASE = Path(os.environ['PCWDE_PROJECT_ROOT']).resolve()
SCRIPT_DIR = Path(__file__).parent
OUTPUT_DIR = SCRIPT_DIR / 'outputs'

# Input files
BLOCKS_FILE = OUTPUT_DIR / 'synteny_blocks.tsv'
LOCI_FILE = OUTPUT_DIR / 'pcwde_loci.tsv'
CHR_CYNOG_FILE = SCRIPT_DIR.parent / 'chromosome_cynog_dict.tsv'
MEMBERSHIP_FILE = BASE / 'good_data/catalogs/cynog/catalog/membership_table.tsv'
REGISTRY_FILE = BASE / 'genomes/pancynipoid_helixer/gca_to_species.tsv'

# Parameters
REGION_BUFFER_GENES = 20      # +/- 20 genes around flanking CynOG cluster
MIN_FLANKING_CYNOGS = 3       # Minimum conserved flanking CynOGs to define a search region
CONSERVED_FREQ = 0.5          # CynOG must appear in >= this fraction of block members
CLUSTER_MAX_GAP = 50          # Max gap (in gene ranks) between flanking CynOGs in a cluster

# Target CynOG families
TARGET_FAMILIES = {'GH5_2', 'PL1_8', 'PL4_4', 'CBM13'}   # family names, not CynOGv3 IDs

# dbCAN domtblout for Tier 2 confirmation (legacy; now superseded by per-family hits TSVs)
DBCAN_FILE = SCRIPT_DIR / 'pcwde_dbcan.domtblout'
DBCAN_EVALUE = 1e-5

# dbCAN per-family hits TSVs — canonical PCWDE target definition
DBCAN_HITS_DIR = BASE / 'analysis/pcwde_effectors/pcwde_paper/data/dbcan_scan'
TARGET_FAMILY_NAMES = ['GH5_2', 'PL1_8', 'PL4_4']

# Protein -> BK chromosome lookup (has bk_chr resolved via JCVI for non-BK genomes)
PROTEIN_TO_CHR_FILE = BASE / 'good_data/catalogs/cynog/catalog/protein_to_chromosome.tsv'

# CynOG -> expected CAZy family prefix
CYNOG_TO_EXPECTED_CAZY = {
    'CynOGv3_000234': 'GH5',
    'CynOGv3_000148': 'PL1',
    'CynOGv3_000168': 'PL4',
    'CynOGv3_010310': 'CBM',
}

# Chromosome mapping
CHR_TO_SCAFFOLD = {
    'chr1': 'CM021338.1', 'chr2': 'CM021339.1', 'chr3': 'CM021340.1',
    'chr4': 'CM021341.1', 'chr5': 'CM021342.1', 'chr6': 'CM021343.1',
    'chr7': 'CM021344.1', 'chr8': 'CM021345.1', 'chr9': 'CM021346.1',
    'chr10': 'CM021347.1',
}


def load_active_genomes():
    """Load active genome accessions from registry."""
    reg = pd.read_csv(REGISTRY_FILE, sep='\t')
    active = reg[reg['status'] == 'active']
    return set(active['accession'])


def load_chr_cynog_dict():
    """Load chromosome_cynog_dict.tsv into indexed structures.

    Returns:
        genes_by_scaffold: {(genome, scaffold): [(rank, protein_id, cynog_id), ...]}
        cynog_positions: {(genome, cynog_id): [(scaffold, rank), ...]}
        genome_scaffolds: {genome: set of scaffolds}
    """
    print("Loading chromosome_cynog_dict.tsv...")
    genes_by_scaffold = defaultdict(list)
    cynog_positions = defaultdict(list)
    genome_scaffolds = defaultdict(set)

    chd = pd.read_csv(CHR_CYNOG_FILE, sep='\t')
    for _, row in chd.iterrows():
        genome = row['genome']
        scaffold = row['chrom']
        rank = int(row['position'])
        protein_id = row['protein_id']
        cynog_id = row['cynog_id'] if pd.notna(row['cynog_id']) else ''

        genes_by_scaffold[(genome, scaffold)].append(
            (rank, protein_id, cynog_id)
        )
        if cynog_id:
            cynog_positions[(genome, cynog_id)].append((scaffold, rank))
        genome_scaffolds[genome].add(scaffold)

    # Sort each scaffold by rank
    for key in genes_by_scaffold:
        genes_by_scaffold[key].sort(key=lambda x: x[0])

    n_genes = sum(len(v) for v in genes_by_scaffold.values())
    print(f"  {n_genes} genes across {len(genome_scaffolds)} genomes")
    return genes_by_scaffold, cynog_positions, genome_scaffolds


import re as _re

def load_dbcan_target_proteins():
    """Same loader as in 01_build_synteny.py — see that file's docstring."""
    print("Loading dbCAN-defined PCWDE target protein sets...")
    targets_per_family = {}
    pid_to_family = {}
    for fam in TARGET_FAMILY_NAMES:
        path = DBCAN_HITS_DIR / f"{fam}_hits.tsv"
        df = pd.read_csv(path, sep='\t')
        pids = set()
        for raw in df['protein_id']:
            if pd.isna(raw):
                continue
            if '|' in raw:
                _, pid_with_ver = raw.split('|', 1)
            else:
                pid_with_ver = raw
            pid = _re.sub(r'(?<=\d)\.\d+$', '', pid_with_ver)
            pids.add(pid)
            pid_to_family[pid] = fam
        targets_per_family[fam] = pids
        print(f"  {fam}: {len(pids)} dbCAN-defined target proteins")
    return targets_per_family, pid_to_family


def parse_dbcan_confirmed(membership_df):
    """Legacy stub — no longer used. Targets now come from dbCAN hits TSVs
    directly. Kept as a placeholder so any vestigial caller doesn't break."""
    print("(parse_dbcan_confirmed is now a no-op; targets come from dbCAN hits TSVs)")
    return set()


def _legacy_parse_dbcan_confirmed(membership_df):
    print("Parsing dbCAN domtblout for Tier 2 confirmation...")

    tier2 = membership_df[
        (membership_df['cynog_id'].isin(TARGET_FAMILIES)) &
        (membership_df['synteny_tier'] == 2)
    ]
    tier2_cynog = dict(zip(tier2['protein_id'], tier2['cynog_id']))

    confirmed = set()
    if not DBCAN_FILE.exists():
        print(f"  WARNING: {DBCAN_FILE} not found")
        return confirmed

    with open(DBCAN_FILE) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 13:
                continue
            cazy_hmm = parts[0].replace('.hmm', '')
            protein_id = parts[3]
            i_evalue = float(parts[12])

            if i_evalue > DBCAN_EVALUE:
                continue
            if protein_id not in tier2_cynog:
                continue

            cynog_id = tier2_cynog[protein_id]
            expected_prefix = CYNOG_TO_EXPECTED_CAZY.get(cynog_id, '')
            if expected_prefix and cazy_hmm.startswith(expected_prefix):
                confirmed.add(protein_id)

    print(f"  dbCAN-confirmed Tier 2: {len(confirmed)}")
    return confirmed


def load_target_members():
    """Load membership table and return target family members with tier info.

    Returns:
        family_map: protein_id -> cynog_id for target families
        tier_map: protein_id -> synteny_tier (1 or 2)
        dbcan_confirmed: set of protein_ids that are dbCAN-confirmed Tier 2
    """
    print("Loading membership_table.tsv for target families...")
    mem = pd.read_csv(MEMBERSHIP_FILE, sep='\t', usecols=['protein_id', 'cynog_id', 'synteny_tier'])
    # NEW: targets defined by dbCAN, not CynOGv3 cluster membership.
    # Build a synthetic 'targets' DataFrame with cynog_id column = family name
    # so the rest of this script (which keys off cynog_id == family) keeps
    # working. The chd / signatures still use real CynOGv3 IDs for flanking.
    targets_per_family, pid_to_family = load_dbcan_target_proteins()
    all_pids = set().union(*targets_per_family.values())
    mem_indexed = mem.set_index('protein_id', drop=False)
    target_rows = []
    for pid in all_pids:
        if pid in mem_indexed.index:
            row = mem_indexed.loc[pid]
            if isinstance(row, pd.DataFrame):
                row = row.iloc[0]
            target_rows.append({
                'protein_id': pid,
                'cynog_id': pid_to_family[pid],   # family name, not real CynOG ID
                'real_cynog_id': row.get('cynog_id', None),
                'genome': row.get('genome', None),
                'synteny_tier': row.get('synteny_tier', 1),
            })
    targets = pd.DataFrame(target_rows)
    print(f"  load_target_members: {len(targets)} dbCAN-defined targets ({len(all_pids) - len(targets)} not in membership table)")
    family_map = dict(zip(targets['protein_id'], targets['cynog_id']))
    tier_map = dict(zip(targets['protein_id'], targets['synteny_tier']))
    dbcan_confirmed = parse_dbcan_confirmed(mem)
    n_t1 = sum(1 for v in tier_map.values() if v == 1)
    n_t2 = sum(1 for v in tier_map.values() if v == 2)
    n_t2c = len(dbcan_confirmed)
    print(f"  {len(family_map)} target family proteins ({n_t1} Tier 1, {n_t2} Tier 2, {n_t2c} Tier 2 dbCAN-confirmed)")
    return family_map, tier_map, dbcan_confirmed


_GENOME_CHR_TO_SCAFFOLD = None
_SCAFFOLD_TO_CHR_NAME = None

def _load_protein_chr_lookups():
    global _GENOME_CHR_TO_SCAFFOLD, _SCAFFOLD_TO_CHR_NAME
    print(f"  Loading protein_to_chromosome.tsv for chromosome resolution...")
    df = pd.read_csv(PROTEIN_TO_CHR_FILE, sep='	', dtype={'bk_chr': str})
    df = df[df['bk_chr'].notna() & (df['bk_chr'] != 'unplaced')]
    g2cs = {}
    for (genome, bk_chr), grp in df.groupby(['genome', 'bk_chr']):
        chr_name = f'chr{bk_chr}'
        modal_scaffold = grp['scaffold'].mode().iloc[0]
        g2cs[(genome, chr_name)] = modal_scaffold
    s2c = {}
    scaff_modes = df.groupby('scaffold')['bk_chr'].agg(
        lambda s: s.mode().iloc[0] if len(s.mode()) else None
    )
    for scaffold, bk in scaff_modes.items():
        if bk and bk != 'unplaced':
            s2c[scaffold] = f'chr{bk}'
            if scaffold.endswith('_RagTag'):
                s2c[scaffold[:-len('_RagTag')]] = f'chr{bk}'
    _GENOME_CHR_TO_SCAFFOLD = g2cs
    _SCAFFOLD_TO_CHR_NAME = s2c
    print(f"  {len(g2cs)} (genome, chr) -> scaffold mappings; {len(s2c)} scaffold -> chr mappings")

_GENOMES_WITH_BK_MAPPING = None

def genome_has_bk_mapping(genome):
    """True if genome has any (genome, chr_name) entry in protein_to_chromosome.tsv.

    False for distant outgroups like Telenomus remus whose BK synteny is too
    eroded for JCVI to assign BK chromosomes. These genomes need all-scaffold
    fallback searching instead of giving up with N_D.
    """
    global _GENOMES_WITH_BK_MAPPING, _GENOME_CHR_TO_SCAFFOLD
    if _GENOME_CHR_TO_SCAFFOLD is None:
        _load_protein_chr_lookups()
    if _GENOMES_WITH_BK_MAPPING is None:
        _GENOMES_WITH_BK_MAPPING = {g for (g, _) in _GENOME_CHR_TO_SCAFFOLD.keys()}
    return genome in _GENOMES_WITH_BK_MAPPING


def find_scaffold_for_genome(genome, chr_name, genome_scaffolds):
    """Find the actual scaffold name for a genome on a given chromosome.

    Tries protein_to_chromosome.tsv first (has bk_chr resolved by JCVI for
    non-BK genomes like DK), falls back to legacy BK CHR_TO_SCAFFOLD pattern.
    """
    global _GENOME_CHR_TO_SCAFFOLD
    if _GENOME_CHR_TO_SCAFFOLD is None:
        _load_protein_chr_lookups()
    scaffold = _GENOME_CHR_TO_SCAFFOLD.get((genome, chr_name))
    if scaffold and scaffold in genome_scaffolds.get(genome, set()):
        return scaffold

    base = CHR_TO_SCAFFOLD.get(chr_name)
    if not base:
        return None

    scaffolds = genome_scaffolds.get(genome, set())
    if base in scaffolds:
        return base
    ragtag = base + '_RagTag'
    if ragtag in scaffolds:
        return ragtag
    for s in scaffolds:
        if base in s:
            return s
    return None


def get_conserved_flanking_cynogs(block_locus_ids, loci_df, col='cynogs_narrow'):
    """Get flanking CynOGs conserved across block members."""
    cynog_counts = Counter()
    n_with_data = 0

    for lid in block_locus_ids:
        rows = loci_df[loci_df['locus_id'] == lid]
        if len(rows) == 0:
            continue
        cynogs_str = rows.iloc[0].get(col, '')
        if pd.isna(cynogs_str) or not cynogs_str:
            continue
        cynogs = set(cynogs_str.split(','))
        cynogs.discard('')
        for c in cynogs:
            cynog_counts[c] += 1
        n_with_data += 1

    if n_with_data == 0:
        return set()

    conserved = {c for c, count in cynog_counts.items() if count >= n_with_data * CONSERVED_FREQ}

    # Fall back to mid if not enough from narrow
    if len(conserved) < MIN_FLANKING_CYNOGS and col == 'cynogs_narrow':
        return get_conserved_flanking_cynogs(block_locus_ids, loci_df, col='cynogs_mid')

    return conserved


def find_region_via_flanking_cynogs(genome, scaffold, conserved_cynogs, cynog_positions):
    """Find the expected locus region using conserved flanking CynOGs.

    Uses gene-rank positions and clustering to find where flanking CynOGs
    are concentrated.

    Returns (start_rank, end_rank) or None if not enough flanking CynOGs found.
    """
    positions_with_cynog = []
    for cynog_id in conserved_cynogs:
        for s, rank in cynog_positions.get((genome, cynog_id), []):
            if s == scaffold:
                positions_with_cynog.append((rank, cynog_id))

    if len(positions_with_cynog) < MIN_FLANKING_CYNOGS:
        return None

    # Cluster by proximity
    positions_with_cynog.sort()
    clusters = [[positions_with_cynog[0]]]
    for i in range(1, len(positions_with_cynog)):
        if positions_with_cynog[i][0] - positions_with_cynog[i - 1][0] > CLUSTER_MAX_GAP:
            clusters.append([])
        clusters[-1].append(positions_with_cynog[i])

    # Pick the cluster with the most distinct CynOGs
    best_cluster = max(clusters, key=lambda c: len({cynog for _, cynog in c}))
    n_distinct = len({cynog for _, cynog in best_cluster})

    if n_distinct < MIN_FLANKING_CYNOGS:
        return None

    cluster_start = best_cluster[0][0]
    cluster_end = best_cluster[-1][0]

    region_start = max(0, cluster_start - REGION_BUFFER_GENES)
    region_end = cluster_end + REGION_BUFFER_GENES

    return (region_start, region_end)


def find_region_any_scaffold(genome, conserved_cynogs, cynog_positions, genome_scaffolds):
    """Find the best conserved-flanking-CynOG cluster across ALL scaffolds of a genome.

    For genomes without BK chromosome mapping (Telenomus remus and similar
    distant outgroups), we can't pin search to a specific 'chr6' scaffold
    because JCVI couldn't synteny-map their contigs to BK. Instead we look
    for a flanking-CynOG cluster on any scaffold and return the best one.

    Returns (scaffold, start_rank, end_rank) or None.
    """
    all_hits = []
    for cynog_id in conserved_cynogs:
        for s, rank in cynog_positions.get((genome, cynog_id), []):
            all_hits.append((s, rank, cynog_id))

    if len(all_hits) < MIN_FLANKING_CYNOGS:
        return None

    by_scaffold = defaultdict(list)
    for s, rank, cynog_id in all_hits:
        by_scaffold[s].append((rank, cynog_id))

    best = None  # (n_distinct, scaffold, start, end)
    for scaffold, hits in by_scaffold.items():
        hits.sort()
        clusters = [[hits[0]]]
        for i in range(1, len(hits)):
            if hits[i][0] - hits[i - 1][0] > CLUSTER_MAX_GAP:
                clusters.append([])
            clusters[-1].append(hits[i])
        for cluster in clusters:
            n_distinct = len({c for _, c in cluster})
            if n_distinct < MIN_FLANKING_CYNOGS:
                continue
            cluster_start = cluster[0][0]
            cluster_end = cluster[-1][0]
            start = max(0, cluster_start - REGION_BUFFER_GENES)
            end = cluster_end + REGION_BUFFER_GENES
            if best is None or n_distinct > best[0]:
                best = (n_distinct, scaffold, start, end)

    if best is None:
        return None
    return (best[1], best[2], best[3])


def scan_region_for_targets(genome, scaffold, start_rank, end_rank,
                            family, genes_by_scaffold, target_members,
                            exclude_pids=frozenset()):
    """Scan a genomic region for target family members.

    Returns list of protein_ids that belong to the target family.
    """
    genes = genes_by_scaffold.get((genome, scaffold), [])
    found = []

    for rank, protein_id, cynog_id in genes:
        if rank < start_rank or rank > end_rank:
            continue
        if protein_id in exclude_pids:
            continue
        if protein_id in target_members and target_members[protein_id] == family:
            found.append(protein_id)

    return found


def main():
    print("=" * 60)
    print("Step 03: Detect Empty Loci (CynOG version)")
    print("=" * 60)

    all_genomes = load_active_genomes()
    print(f"Registry: {len(all_genomes)} active genomes")

    genes_by_scaffold, cynog_positions, genome_scaffolds = load_chr_cynog_dict()
    target_members, tier_map, dbcan_confirmed = load_target_members()

    blocks_df = pd.read_csv(BLOCKS_FILE, sep='\t')
    loci_df = pd.read_csv(LOCI_FILE, sep='\t')

    locus_genome = dict(zip(loci_df['locus_id'], loci_df['genome']))

    block_entries = blocks_df[blocks_df['status'] == 'in_block'].copy()
    singleton_entries = blocks_df[blocks_df['status'] == 'true_singleton'].copy()

    unique_blocks = sorted(block_entries['block_id'].unique())
    print(f"\n{len(unique_blocks)} synteny blocks to analyze")
    print(f"{len(singleton_entries)} singletons")

    # Build set of ALL protein IDs that are block members
    all_member_pids = set()
    for _, row in block_entries.iterrows():
        pids = [p.strip() for p in str(row['protein_ids']).split(',') if p.strip()]
        all_member_pids.update(pids)
    print(f"{len(all_member_pids)} total member proteins across all blocks")

    results = []

    # === PROCESS BLOCKS ===
    for block_id in unique_blocks:
        block = block_entries[block_entries['block_id'] == block_id]
        family = block['family'].iloc[0]
        modal_chr = block['modal_chr'].iloc[0]

        member_genomes = {}
        block_locus_ids = []

        for _, row in block.iterrows():
            lid = row['locus_id']
            genome = locus_genome.get(lid)
            if genome is None:
                print(f"  WARNING: no genome found for locus {lid}", file=sys.stderr)
                continue
            pids = [p.strip() for p in str(row['protein_ids']).split(',') if p.strip()]
            if genome not in member_genomes:
                member_genomes[genome] = []
            member_genomes[genome].extend(pids)
            block_locus_ids.append(lid)

        conserved_cynogs = get_conserved_flanking_cynogs(block_locus_ids, loci_df)

        print(f"\n{block_id} ({family}, {modal_chr}): "
              f"{len(member_genomes)} member genomes, "
              f"{len(conserved_cynogs)} conserved flanking CynOGs")

        # Member genomes: PRESENT
        for genome, pids in member_genomes.items():
            scaffold = find_scaffold_for_genome(genome, modal_chr, genome_scaffolds)
            results.append({
                'region_id': f"{block_id}__{genome}",
                'region_type': 'block',
                'block_id': block_id,
                'family': family,
                'genome': genome,
                'scaffold': scaffold or '',
                'status': 'PRESENT',
                'n_tandem': len(pids),
                'target_gene_ids': ','.join(pids),
                'region_start_rank': '',
                'region_end_rank': '',
            })

        # Non-member genomes
        for genome in sorted(all_genomes):
            if genome in member_genomes:
                continue

            scaffold = find_scaffold_for_genome(genome, modal_chr, genome_scaffolds)
            region = None

            # Fallback: if no per-chromosome scaffold mapping exists for this
            # (genome, modal_chr), search the genome's scaffolds for a cluster
            # of conserved flanking CynOGs anywhere. Handles distant outgroups
            # like Telenomus remus where JCVI synteny is too eroded to place
            # whole chromosomes, as well as fragmented assemblies.
            if not scaffold:
                any_result = find_region_any_scaffold(
                    genome, conserved_cynogs, cynog_positions, genome_scaffolds
                )
                if any_result is not None:
                    scaffold, rstart, rend = any_result
                    region = (rstart, rend)

            if not scaffold:
                results.append({
                    'region_id': f"{block_id}__{genome}",
                    'region_type': 'block',
                    'block_id': block_id,
                    'family': family,
                    'genome': genome,
                    'scaffold': '',
                    'status': 'N_D',
                    'n_tandem': 0,
                    'target_gene_ids': '',
                    'region_start_rank': '',
                    'region_end_rank': '',
                })
                continue

            if region is None:
                region = find_region_via_flanking_cynogs(
                    genome, scaffold, conserved_cynogs, cynog_positions
                )

            if not region:
                results.append({
                    'region_id': f"{block_id}__{genome}",
                    'region_type': 'block',
                    'block_id': block_id,
                    'family': family,
                    'genome': genome,
                    'scaffold': scaffold,
                    'status': 'N_D',
                    'n_tandem': 0,
                    'target_gene_ids': '',
                    'region_start_rank': '',
                    'region_end_rank': '',
                })
                continue

            search_start, search_end = region

            hits = scan_region_for_targets(
                genome, scaffold, search_start, search_end, family,
                genes_by_scaffold, target_members, exclude_pids=all_member_pids
            )

            if hits:
                has_confirmed = any(
                    tier_map.get(pid, 2) == 1 or pid in dbcan_confirmed
                    for pid in hits
                )
                if has_confirmed:
                    status = 'PRESENT'
                else:
                    status = 'UNCONFIRMED'
                n_tandem = len(hits)
                target_ids = ','.join(hits)
            else:
                status = 'EMPTY'
                n_tandem = 0
                target_ids = ''

            results.append({
                'region_id': f"{block_id}__{genome}",
                'region_type': 'block',
                'block_id': block_id,
                'family': family,
                'genome': genome,
                'scaffold': scaffold,
                'status': status,
                'n_tandem': n_tandem,
                'target_gene_ids': target_ids,
                'region_start_rank': search_start,
                'region_end_rank': search_end,
            })

    # === PROCESS SINGLETONS ===
    def _scaffold_to_chr(scaffold):
        global _SCAFFOLD_TO_CHR_NAME
        if _SCAFFOLD_TO_CHR_NAME is None:
            _load_protein_chr_lookups()
        if scaffold in _SCAFFOLD_TO_CHR_NAME:
            return _SCAFFOLD_TO_CHR_NAME[scaffold]
        base = scaffold.replace('_RagTag', '')
        if base in _SCAFFOLD_TO_CHR_NAME:
            return _SCAFFOLD_TO_CHR_NAME[base]
        for chr_name, scaf in CHR_TO_SCAFFOLD.items():
            if scaf == base:
                return chr_name
        return None

    n_chr_singletons = 0
    n_frag_singletons = 0

    for _, row in singleton_entries.iterrows():
        lid = row['locus_id']
        genome = locus_genome.get(lid)
        if genome is None:
            continue

        loci_row = loci_df[loci_df['locus_id'] == lid]
        scaffold = loci_row.iloc[0]['scaffold'] if len(loci_row) > 0 else ''
        family = row['family']
        pids = [p.strip() for p in str(row['protein_ids']).split(',') if p.strip()]
        chr_name = _scaffold_to_chr(scaffold) if scaffold else None

        sing_block_id = f"SINGLE_{lid}"
        results.append({
            'region_id': f"SINGLE__{genome}__{lid}",
            'region_type': 'singleton',
            'block_id': sing_block_id,
            'family': family,
            'genome': genome,
            'scaffold': scaffold,
            'status': 'PRESENT',
            'n_tandem': len(pids),
            'target_gene_ids': ','.join(pids),
            'region_start_rank': '',
            'region_end_rank': '',
        })

        if chr_name is None:
            n_frag_singletons += 1
            continue

        n_chr_singletons += 1

        # Use this singleton's own flanking CynOGs
        cynogs_str = loci_row.iloc[0].get('cynogs_narrow', '') if len(loci_row) > 0 else ''
        if pd.isna(cynogs_str) or not cynogs_str:
            cynogs_str = loci_row.iloc[0].get('cynogs_mid', '') if len(loci_row) > 0 else ''
        if pd.isna(cynogs_str) or not cynogs_str:
            for other_genome in sorted(all_genomes):
                if other_genome == genome:
                    continue
                results.append({
                    'region_id': f"SINGLE__{other_genome}__{lid}",
                    'region_type': 'singleton',
                    'block_id': sing_block_id,
                    'family': family,
                    'genome': other_genome,
                    'scaffold': '',
                    'status': 'N_D',
                    'n_tandem': 0,
                    'target_gene_ids': '',
                    'region_start_rank': '',
                    'region_end_rank': '',
                })
            continue

        singleton_cynogs = set(cynogs_str.split(','))
        singleton_cynogs.discard('')

        for other_genome in sorted(all_genomes):
            if other_genome == genome:
                continue

            other_scaffold = find_scaffold_for_genome(other_genome, chr_name, genome_scaffolds)
            sing_region = None
            if not other_scaffold:
                any_result = find_region_any_scaffold(
                    other_genome, singleton_cynogs, cynog_positions, genome_scaffolds
                )
                if any_result is not None:
                    other_scaffold, _rs, _re = any_result
                    sing_region = (_rs, _re)
            if not other_scaffold:
                results.append({
                    'region_id': f"SINGLE__{other_genome}__{lid}",
                    'region_type': 'singleton',
                    'block_id': sing_block_id,
                    'family': family,
                    'genome': other_genome,
                    'scaffold': '',
                    'status': 'N_D',
                    'n_tandem': 0,
                    'target_gene_ids': '',
                    'region_start_rank': '',
                    'region_end_rank': '',
                })
                continue

            if sing_region is not None:
                region = sing_region
            else:
                region = find_region_via_flanking_cynogs(
                    other_genome, other_scaffold, singleton_cynogs, cynog_positions
                )

            if not region:
                results.append({
                    'region_id': f"SINGLE__{other_genome}__{lid}",
                    'region_type': 'singleton',
                    'block_id': sing_block_id,
                    'family': family,
                    'genome': other_genome,
                    'scaffold': other_scaffold,
                    'status': 'N_D',
                    'n_tandem': 0,
                    'target_gene_ids': '',
                    'region_start_rank': '',
                    'region_end_rank': '',
                })
                continue

            search_start, search_end = region
            hits = scan_region_for_targets(
                other_genome, other_scaffold, search_start, search_end, family,
                genes_by_scaffold, target_members, exclude_pids=all_member_pids
            )

            if hits:
                has_confirmed = any(
                    tier_map.get(pid, 2) == 1 or pid in dbcan_confirmed
                    for pid in hits
                )
                if has_confirmed:
                    status = 'PRESENT'
                else:
                    status = 'UNCONFIRMED'
                n_tandem = len(hits)
                target_ids = ','.join(hits)
            else:
                status = 'EMPTY'
                n_tandem = 0
                target_ids = ''

            results.append({
                'region_id': f"SINGLE__{other_genome}__{lid}",
                'region_type': 'singleton',
                'block_id': sing_block_id,
                'family': family,
                'genome': other_genome,
                'scaffold': other_scaffold,
                'status': status,
                'n_tandem': n_tandem,
                'target_gene_ids': target_ids,
                'region_start_rank': search_start,
                'region_end_rank': search_end,
            })

    print(f"\nSingletons: {n_chr_singletons} on chromosomes (searched), "
          f"{n_frag_singletons} on fragments (owner only)")

    # === OUTPUT ===
    result_df = pd.DataFrame(results)
    result_df = result_df.sort_values(['family', 'block_id', 'genome'])

    output_file = OUTPUT_DIR / 'region_status.tsv'
    result_df.to_csv(output_file, sep='\t', index=False)
    print(f"\nSaved {len(result_df)} region statuses to {output_file}")

    # Summary
    print("\n" + "=" * 60)
    print("STATUS SUMMARY")
    print("=" * 60)

    for status in ['PRESENT', 'UNCONFIRMED', 'EMPTY', 'N_D']:
        count = (result_df['status'] == status).sum()
        print(f"  {status}: {count}")

    block_results = result_df[result_df['region_type'] == 'block']
    print(f"\nPer block:")
    for block_id in sorted(block_results['block_id'].unique()):
        bdf = block_results[block_results['block_id'] == block_id]
        family = bdf['family'].iloc[0]
        n_present = (bdf['status'] == 'PRESENT').sum()
        n_empty = (bdf['status'] == 'EMPTY').sum()
        n_nd = (bdf['status'] == 'N_D').sum()
        print(f"  {block_id} ({family}): {n_present} present, {n_empty} empty, {n_nd} n/d")

    sing_results = result_df[result_df['region_type'] == 'singleton']
    n_sing = len(sing_results)
    n_sing_present = (sing_results['status'] == 'PRESENT').sum()
    n_sing_empty = (sing_results['status'] == 'EMPTY').sum()
    n_sing_nd = (sing_results['status'] == 'N_D').sum()
    print(f"\nSingletons: {n_sing} total ({n_sing_present} PRESENT, "
          f"{n_sing_empty} EMPTY, {n_sing_nd} N_D)")


if __name__ == '__main__':
    main()
