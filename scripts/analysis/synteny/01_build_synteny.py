#!/usr/bin/env python3
"""
Build PCWDE target-centered synteny analysis using CynOG assignments.

Uses flanking gene CynOGs from chromosome_cynog_dict to identify syntenic loci
across genomes. Runs parameter sweep over window sizes to assign confidence.

Inputs:
  - chromosome_cynog_dict.tsv: genome-wide ordered gene list with CynOG assignments
  - membership_table.tsv: protein -> CynOG mapping
  - gca_to_species.tsv: accession -> species name mapping

Outputs:
  - synteny_pairs_{narrow,mid,wide}.tsv: per-window pairwise results
  - synteny_pairs_merged.tsv: merged with confidence tiers
  - synteny_blocks.tsv: block assignments with status and modal_chr
  - pcwde_loci.tsv: per-locus CynOG signatures
"""

import pandas as pd
import numpy as np
import os
from pathlib import Path
from itertools import combinations
from collections import Counter
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

# === CONFIGURATION ===

BASE = Path(os.environ['PCWDE_PROJECT_ROOT']).resolve()
SCRIPT_DIR = Path(__file__).parent
CHR_CYNOG_DICT = SCRIPT_DIR.parent / 'chromosome_cynog_dict.tsv'
MEMBERSHIP = BASE / 'good_data/catalogs/cynog/catalog/membership_table.tsv'
GCA_SPECIES = BASE / 'genomes/pancynipoid_helixer/gca_to_species.tsv'
OUTPUT_DIR = SCRIPT_DIR / 'outputs'
OUTPUT_DIR.mkdir(exist_ok=True)

# Target PCWDE families (CynOG equivalents)
TARGET_FAMILIES = ['GH5_2', 'PL1_8', 'PL4_4', 'CBM13']   # family names; downstream code matches against this

# Window sizes: (name, n_upstream, n_downstream)
WINDOWS = [
    ('narrow', 5, 5),
    ('mid', 15, 15),
    ('wide', 50, 50),
]

JACCARD_THRESHOLD = 0.15

# dbCAN domtblout for Tier 2 confirmation
DBCAN_FILE = SCRIPT_DIR / 'pcwde_dbcan.domtblout'
DBCAN_EVALUE = 1e-5

# Path to dbCAN-defined family hit TSVs (canonical PCWDE definition).
# These files list every protein the dbCAN HMM scan flagged for each family,
# regardless of which CynOGv3 cluster the family catalog assigned them to.
DBCAN_HITS_DIR = BASE / 'analysis/pcwde_effectors/pcwde_paper/data/dbcan_scan'

# Path to the protein -> BK chromosome lookup. Has bk_chr resolved via JCVI
# synteny for non-BK genomes (Dryocosmus kuriphilus, etc.), so it covers
# scaffolds that don't match BK's CM02xxxx naming pattern.
PROTEIN_TO_CHR_FILE = BASE / 'good_data/catalogs/cynog/catalog/protein_to_chromosome.tsv'

# Family names — the canonical "what is this gene" axis. The targets
# DataFrame stores the family name in its cynog_id column for compat with
# downstream code, but downstream code's "cynog_id == fam" is now matching
# against these names, not actual CynOGv3 cluster IDs.
TARGET_FAMILY_NAMES = ['GH5_2', 'PL1_8', 'PL4_4', 'CBM13']

# CynOG -> expected CAZy family prefix (for dbCAN confirmation; legacy)
CYNOG_TO_EXPECTED_CAZY = {
    'CynOGv3_000234': 'GH5',
    'CynOGv3_000148': 'PL1',
    'CynOGv3_000168': 'PL4',
    'CynOGv3_010310': 'CBM',
}

# Adjacent target proteins with <= this many intervening genes form one locus
MAX_LOCUS_GAP = 5

# Corrected-source runs supply accepted targets absent from gene order here.
# An unset path preserves the frozen historical input and output contract.
MISSING_ORDER_FILE = os.environ.get('PCWDE_MISSING_ORDER_TSV')
RECONSTRUCTED_ORDER_FILE = os.environ.get('PCWDE_RECONSTRUCTED_ORDER_TSV')
EXCLUDED_GENOMES = set(filter(None, os.environ.get('PCWDE_EXCLUDED_GENOMES', '').split(',')))
# A production caller may supply a source-derived accepted-copy eligibility set.
# Empty by default, preserving the original native builder contract.
INELIGIBLE_ACCEPTED_IDS = set()

# Scaffold -> chromosome mapping (BK reference-based RagTag scaffolding)
SCAFFOLD_TO_CHR = {
    'CM021338.1': 'chr1',
    'CM021339.1': 'chr2',
    'CM021340.1': 'chr3',
    'CM021341.1': 'chr4',
    'CM021342.1': 'chr5',
    'CM021343.1': 'chr6',
    'CM021344.1': 'chr7',
    'CM021345.1': 'chr8',
    'CM021346.1': 'chr9',
    'CM021347.1': 'chr10',
}


_SCAFFOLD_TO_BK_CHR = None

def _load_scaffold_to_bk_chr():
    """Build {scaffold_name: 'chrN'} from protein_to_chromosome.tsv.

    For each scaffold, take the modal bk_chr value across all proteins on it
    (in case a single scaffold has both placed and unplaced rows). Return
    'fragment' for scaffolds where the modal value is 'unplaced' or empty.
    Both bare ('CM021343.1') and RagTag ('CM021343.1_RagTag') variants are
    populated so callers don't need to normalize.
    """
    global _SCAFFOLD_TO_BK_CHR
    print(f"  Loading scaffold->bk_chr lookup from {PROTEIN_TO_CHR_FILE.name}...")
    df = pd.read_csv(PROTEIN_TO_CHR_FILE, sep='\t', dtype={'bk_chr': str})
    # Modal bk_chr per scaffold
    modes = df.groupby('scaffold')['bk_chr'].agg(
        lambda s: s.mode().iloc[0] if len(s.mode()) > 0 else None
    )
    table = {}
    for scaffold, bk in modes.items():
        if bk is None or bk == 'unplaced' or pd.isna(bk):
            chrom = 'fragment'
        else:
            chrom = f'chr{bk}'
        table[scaffold] = chrom
        # Also populate the bare variant if this is a RagTag-suffixed entry
        if scaffold.endswith('_RagTag'):
            table[scaffold[:-len('_RagTag')]] = chrom
    _SCAFFOLD_TO_BK_CHR = table
    n_placed = sum(1 for v in table.values() if v != 'fragment')
    print(f"  Loaded {len(table)} scaffold mappings ({n_placed} placed on a chromosome)")
    return table

def scaffold_to_chromosome(scaffold):
    """Map scaffold name to BK chromosome via protein_to_chromosome.tsv.

    Falls back to 'fragment' for unplaced scaffolds (small NODE_*_length_*
    contigs and the like).
    """
    global _SCAFFOLD_TO_BK_CHR
    if _SCAFFOLD_TO_BK_CHR is None:
        _load_scaffold_to_bk_chr()
    if scaffold in _SCAFFOLD_TO_BK_CHR:
        return _SCAFFOLD_TO_BK_CHR[scaffold]
    base = scaffold.replace('_RagTag', '')
    return _SCAFFOLD_TO_BK_CHR.get(base, 'fragment')


import re as _re

def load_dbcan_target_proteins():
    """Load dbCAN-defined PCWDE target protein sets per family.

    Reads data/dbcan_scan/{family}_hits.tsv files. Each row's protein_id is
    formatted '{genome_key}|{protein_id}.{version}'. We strip the genome
    prefix and the .N version suffix to obtain the canonical protein_id used
    by chromosome_cynog_dict and protein_to_chromosome.tsv.

    Returns:
      targets_per_family: {family_name: set(protein_id_no_ver)}
      pid_to_family: {protein_id_no_ver: family_name}
    """
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
            # Strip "{genome}|" prefix
            if '|' in raw:
                _, pid_with_ver = raw.split('|', 1)
            else:
                pid_with_ver = raw
            # Strip version suffix (.1, .2, ...)
            pid = _re.sub(r'(?<=\d)\.\d+$', '', pid_with_ver)
            pids.add(pid)
            pid_to_family[pid] = fam
        targets_per_family[fam] = pids
        print(f"  {fam}: {len(pids)} dbCAN-defined target proteins")
    return targets_per_family, pid_to_family


def load_data():
    """Load all required data files.

    Targets are now defined by the dbCAN scan, NOT by CynOGv3 cluster
    membership. The targets DataFrame stores the family name (e.g. 'GH5_2')
    in its cynog_id column so downstream code that does
    `targets['cynog_id'] == fam` keeps working — but the column semantics
    have shifted from "actual CynOGv3 cluster ID" to "family axis label."
    The chromosome_cynog_dict that the find_loci index uses still has actual
    CynOGv3 IDs in its cynog_id column for flanking gene comparisons.
    """
    print("Loading data...")

    # Chromosome CynOG dict: ordered genes per chromosome with CynOG assignments.
    # This stores REAL CynOGv3 IDs in its cynog_id column — used for flanking
    # gene comparisons (the orthology axis). Target identification no longer
    # filters this table by cynog_id.
    chd = pd.read_csv(CHR_CYNOG_DICT, sep='\t')
    if EXCLUDED_GENOMES:
        absent = EXCLUDED_GENOMES - set(chd['genome'])
        if absent:
            raise ValueError(f'excluded genomes absent from order: {sorted(absent)}')
        chd = chd[~chd['genome'].isin(EXCLUDED_GENOMES)].copy()
    print(f"  chromosome_cynog_dict: {len(chd):,} genes across {chd['genome'].nunique()} genomes")

    # dbCAN-defined target proteins per family
    targets_per_family, pid_to_family = load_dbcan_target_proteins()
    all_target_pids = set().union(*targets_per_family.values())

    # Cross-reference target protein IDs against chd to keep only those we
    # can place on a scaffold (i.e. that have a position for flanking lookup).
    # Targets that aren't in chd at all are skipped — these are mostly small
    # unplaced contigs that wouldn't have meaningful flanking context anyway.
    chd_indexed = chd.set_index('protein_id', drop=False)
    in_chd = [pid for pid in all_target_pids if pid in chd_indexed.index]
    missing_chd = len(all_target_pids) - len(in_chd)
    print(f"  Targets in chromosome_cynog_dict: {len(in_chd)} ({missing_chd} missing — small fragments / version mismatches)")

    target_rows = []
    for pid in in_chd:
        row = chd_indexed.loc[pid]
        if isinstance(row, pd.DataFrame):
            row = row.iloc[0]
        target_rows.append({
            'protein_id': pid,
            'cynog_id': pid_to_family[pid],   # family name, not real CynOG ID
            'genome': row['genome'],
            'chrom': row['chrom'],
            'position': row['position'],
            'real_cynog_id': row['cynog_id'],  # actual CynOG, kept for diagnostics
        })
    targets = pd.DataFrame(target_rows)
    if MISSING_ORDER_FILE:
        missing = pd.read_csv(MISSING_ORDER_FILE, sep='\t', dtype=str).fillna('')
        if set(missing['genome']) & EXCLUDED_GENOMES:
            raise ValueError('missing-order input includes an excluded genome')
        required = {'protein_id', 'family', 'genome', 'scaffold', 'start_bp',
                    'end_bp', 'chromosome', 'reason'}
        if not required.issubset(missing.columns):
            raise ValueError(f'missing-order input lacks {sorted(required - set(missing.columns))}')
        if missing.duplicated(['protein_id', 'scaffold']).any():
            raise ValueError('duplicate missing-order protein/scaffold segment')
        if missing.groupby('protein_id')[['family', 'genome']].nunique().gt(1).any().any():
            raise ValueError('inconsistent multi-segment target identity')
        missing = missing.loc[~missing.protein_id.isin(INELIGIBLE_ACCEPTED_IDS)].copy()
        if set(missing['protein_id']) & set(targets['protein_id']):
            raise ValueError('missing-order input includes a joined target')
        for row in missing.itertuples(index=False):
            if row.family not in TARGET_FAMILIES or row.protein_id not in targets_per_family[row.family]:
                raise ValueError(f'not an accepted family target: {row.protein_id}')
            if not row.genome or not row.scaffold or not row.reason:
                raise ValueError(f'missing source identity or reason: {row.protein_id}')
            if int(row.start_bp) < 1 or int(row.end_bp) < int(row.start_bp):
                raise ValueError(f'invalid source coordinates: {row.protein_id}')
            if row.chromosome not in {'unplaced', *(f'chr{i}' for i in range(1, 11))}:
                raise ValueError(f'invalid BK placement: {row.protein_id}')
            if row.protein_id in {r['protein_id'] for r in target_rows}:
                continue  # one accepted target; its other scaffold spans remain in the order
            target_rows.append({
                'protein_id': row.protein_id, 'cynog_id': row.family,
                'genome': row.genome, 'chrom': row.scaffold,
                'position': np.nan, 'real_cynog_id': np.nan,
                'start_bp': int(row.start_bp), 'end_bp': int(row.end_bp),
                'chromosome': row.chromosome, 'reason': row.reason,
            })
        targets = pd.DataFrame(target_rows)

    print(f"  Target proteins (placeable): {len(targets)}")
    for fam in TARGET_FAMILIES:
        n = (targets['cynog_id'] == fam).sum()
        print(f"    {fam}: {n}")

    # GCA -> species mapping
    gca_df = pd.read_csv(GCA_SPECIES, sep='\t')
    acc_to_species = dict(zip(gca_df['accession'], gca_df['species'].str.replace(' ', '_')))
    print(f"  Species mapping: {len(acc_to_species)} entries")

    return chd, targets, acc_to_species


def build_gene_index(chd):
    """Index chromosome_cynog_dict by (genome, chrom) for fast lookup.

    Returns dict: (genome, chrom) -> list of (position, protein_id, cynog_id)
                  sorted by position.
    """
    print("\nBuilding gene index...")
    gene_index = {}
    for (genome, chrom), group in chd.groupby(['genome', 'chrom']):
        genes = list(zip(group['position'], group['protein_id'], group['cynog_id']))
        genes.sort(key=lambda x: x[0])
        gene_index[(genome, chrom)] = genes
    print(f"  Indexed {len(gene_index)} genome-scaffold combinations")
    return gene_index


def inject_reconstructed_order(chd):
    """Use full source-annotation order on affected scaffolds, with an index gate.

    The optional file is produced by the historical v3 admission/order rule.
    Existing indexed proteins must retain their relative order. Affected
    scaffolds below the landmark filter remain available for locus grouping.
    """
    if not RECONSTRUCTED_ORDER_FILE:
        return chd
    order = pd.read_csv(RECONSTRUCTED_ORDER_FILE, sep='\t', dtype=str).fillna('')
    required = {'genome', 'scaffold', 'position', 'protein_id'}
    if not required.issubset(order.columns):
        raise ValueError(f'reconstructed order lacks {sorted(required - set(order.columns))}')
    if order.duplicated(['genome', 'scaffold', 'position']).any():
        raise ValueError('duplicate reconstructed scaffold position')
    if order.duplicated(['genome', 'scaffold', 'protein_id']).any():
        raise ValueError('duplicate reconstructed scaffold protein')
    order['position'] = pd.to_numeric(order['position'], errors='raise').astype(int)
    if (order['position'] < 0).any():
        raise ValueError('negative reconstructed position')
    affected = set(zip(order['genome'], order['scaffold']))
    old = chd[chd[['genome', 'chrom']].apply(tuple, axis=1).isin(affected)]
    for key, group in old.groupby(['genome', 'chrom']):
        original = group.sort_values('position')['protein_id'].astype(str).tolist()
        rebuilt = order[(order['genome'] == key[0]) & (order['scaffold'] == key[1])].sort_values('position')
        rebuilt_ids = rebuilt['protein_id'].tolist()
        old_ids = set(original)
        # AE/AQP deliberately replace the prior native annotation. Their
        # historical protein IDs need not occur in the corrected source GFF;
        # the source-specific admission/coordinate gate is upstream.
        replacement_source = key[0] in {'GCA_035041075.1', 'GCA_035222095.1'}
        if not replacement_source and [pid for pid in rebuilt_ids if pid in old_ids] != original:
            raise ValueError(f'reconstructed order disagrees with existing index: {key}')
    # Retain CynOG labels for historical entries. New admitted models may lack
    # a CynOG assignment in this paper-layer export; they still occupy a rank.
    labels = chd.drop_duplicates(['genome', 'protein_id']).set_index(['genome', 'protein_id'])['cynog_id']
    replacement = order.rename(columns={'scaffold': 'chrom'})[['genome', 'chrom', 'position', 'protein_id']].copy()
    replacement['cynog_id'] = [labels.get((g, p), np.nan) for g, p in zip(replacement['genome'], replacement['protein_id'])]
    retained = chd[~chd[['genome', 'chrom']].apply(tuple, axis=1).isin(affected)]
    return pd.concat([retained, replacement], ignore_index=True)


def find_loci(gene_index, targets, acc_to_species):
    """Identify loci by grouping adjacent target proteins on the same scaffold."""
    print("\nIdentifying loci...")

    if RECONSTRUCTED_ORDER_FILE:
        return find_loci_segment_aware(gene_index, targets, acc_to_species)

    target_set = set(targets['protein_id'])
    target_family = dict(zip(targets['protein_id'], targets['cynog_id']))

    loci = []

    for (genome, chrom), genes in gene_index.items():
        target_positions = []
        for i, (pos, pid, cynog_id) in enumerate(genes):
            if pid in target_set:
                target_positions.append((i, pos, pid, target_family[pid]))

        if not target_positions:
            continue

        current_group = [target_positions[0]]
        for tp in target_positions[1:]:
            prev_idx = current_group[-1][0]
            curr_idx = tp[0]
            gap = curr_idx - prev_idx - 1

            if gap <= MAX_LOCUS_GAP and tp[3] == current_group[0][3]:
                current_group.append(tp)
            else:
                loci.append(_make_locus(current_group, genome, chrom, acc_to_species))
                current_group = [tp]

        loci.append(_make_locus(current_group, genome, chrom, acc_to_species))


    loci_df = pd.DataFrame(loci)
    print(f"  Found {len(loci_df)} loci")
    for fam in TARGET_FAMILIES:
        n = (loci_df['family'] == fam).sum()
        print(f"    {fam}: {n} loci")
    return loci_df


def find_loci_segment_aware(gene_index, targets, acc_to_species):
    """Group target units using every reconstructed scaffold segment.

    A protein can occupy multiple source scaffolds. Edges arise only within a
    scaffold; connected proteins form one locus and each protein is counted once.
    Cross-scaffold articulation is audited before this producer is used.
    """
    identity = {(r.genome, r.protein_id): r.cynog_id for r in targets.itertuples()}
    parent = {key: key for key in identity}

    def root(key):
        while parent[key] != key:
            parent[key] = parent[parent[key]]
            key = parent[key]
        return key

    def join(a, b):
        a, b = root(a), root(b)
        if a != b:
            parent[max(a, b)] = min(a, b)

    occurrences = {}
    for (genome, chrom), genes in gene_index.items():
        by_family = {}
        for i, (pos, pid, _) in enumerate(genes):
            key = (genome, pid)
            if key in identity:
                by_family.setdefault(identity[key], []).append((i, pos, key))
                occurrences.setdefault(key, []).append((genome, chrom, i, pos))
        for entries in by_family.values():
            for left, right in zip(entries, entries[1:]):
                if right[0] - left[0] - 1 <= MAX_LOCUS_GAP:
                    join(left[2], right[2])

    absent = set(identity) - set(occurrences)
    if absent:
        raise ValueError(f'accepted targets absent from ordered index: {len(absent)}')
    components = {}
    for key in identity:
        components.setdefault(root(key), []).append(key)
    loci = []
    for members in components.values():
        families = {identity[key] for key in members}
        if len(families) != 1:
            raise ValueError('cross-family locus component')
        # Anchor to the scaffold representing most distinct component members.
        counts = Counter((g, c) for key in members for g, c, _, _ in occurrences[key])
        genome, chrom = sorted(counts, key=lambda k: (-counts[k], k))[0]
        group = []
        for key in members:
            for g, c, i, pos in occurrences[key]:
                if (g, c) == (genome, chrom):
                    group.append((i, pos, key[1], identity[key]))
        group.sort()
        locus = _make_locus(group, genome, chrom, acc_to_species)
        locus['n_genes'] = len(members)
        locus['protein_ids'] = ','.join(sorted(key[1] for key in members))
        loci.append(locus)
    loci_df = pd.DataFrame(loci).sort_values('locus_id').reset_index(drop=True)
    if sum(loci_df['n_genes']) != len(identity):
        raise ValueError('target-unit locus accounting is not bijective')
    print(f"  Found {len(loci_df)} loci from {len(identity)} accepted targets")
    return loci_df


def _make_locus(group, genome, chrom, acc_to_species):
    """Create a locus dict from a group of adjacent target positions."""
    indices = [g[0] for g in group]
    positions = [g[1] for g in group]
    protein_ids = [g[2] for g in group]
    family = group[0][3]
    species = acc_to_species.get(genome, genome)
    chromosome = scaffold_to_chromosome(chrom)

    locus_id = f"{genome}_{chrom}_pos{min(positions)}"
    if len(positions) > 1:
        locus_id += f"_to_{max(positions)}"

    return {
        'locus_id': locus_id,
        'family': family,
        'genome': genome,
        'scaffold': chrom,
        'chromosome': chromosome,
        'species': species,
        'n_genes': len(group),
        'protein_ids': ','.join(protein_ids),
        'gene_indices': indices,
        'min_pos': min(positions),
        'max_pos': max(positions),
    }


def build_cynog_signatures(loci_df, gene_index):
    """Build CynOG signatures for each locus at each window size."""
    print("\nBuilding CynOG signatures...")

    signatures = {}

    for _, locus in loci_df.iterrows():
        lid = locus['locus_id']
        genome = locus['genome']
        chrom = locus['scaffold']
        indices = locus['gene_indices']

        genes = gene_index.get((genome, chrom), [])
        if not genes:
            signatures[lid] = {w[0]: set() for w in WINDOWS}
            continue

        min_idx = min(indices)
        max_idx = max(indices)
        locus_indices = set(range(min_idx, max_idx + 1))

        signatures[lid] = {}

        for window_name, n_up, n_down in WINDOWS:
            cynogs = set()

            # Upstream
            count = 0
            i = min_idx - 1
            while i >= 0 and count < n_up:
                if i not in locus_indices:
                    cynog_id = genes[i][2]
                    if pd.notna(cynog_id) and cynog_id != '':
                        cynogs.add(cynog_id)
                    count += 1
                i -= 1

            # Downstream
            count = 0
            i = max_idx + 1
            while i < len(genes) and count < n_down:
                if i not in locus_indices:
                    cynog_id = genes[i][2]
                    if pd.notna(cynog_id) and cynog_id != '':
                        cynogs.add(cynog_id)
                    count += 1
                i += 1

            signatures[lid][window_name] = cynogs

    print(f"  Built signatures for {len(signatures)} loci")
    for window_name, _, _ in WINDOWS:
        sizes = [len(signatures[l][window_name]) for l in signatures]
        if sizes:
            print(f"    {window_name}: median {np.median(sizes):.0f} CynOGs, range {min(sizes)}-{max(sizes)}")

    return signatures


def compute_pairwise_synteny(signatures, loci_df, threshold):
    """Compute pairwise Jaccard similarity between loci of the same family."""
    print("\nComputing pairwise synteny...")

    locus_family = dict(zip(loci_df['locus_id'], loci_df['family']))
    results = {w[0]: [] for w in WINDOWS}

    for family in TARGET_FAMILIES:
        fam_loci = [lid for lid, fam in locus_family.items() if fam == family]
        n_pairs = len(fam_loci) * (len(fam_loci) - 1) // 2
        print(f"  {family}: {len(fam_loci)} loci -> {n_pairs} pairs")

        for i, (l1, l2) in enumerate(combinations(fam_loci, 2)):
            if i % 5000 == 0 and i > 0:
                print(f"    Processed {i}/{n_pairs} pairs...")

            for window_name, _, _ in WINDOWS:
                s1 = signatures[l1][window_name]
                s2 = signatures[l2][window_name]

                if len(s1) == 0 or len(s2) == 0:
                    continue

                intersection = len(s1 & s2)
                union = len(s1 | s2)
                jaccard = intersection / union if union > 0 else 0

                if jaccard >= threshold:
                    results[window_name].append({
                        'locus_1': l1,
                        'locus_2': l2,
                        'jaccard': jaccard,
                        'shared_cynogs': intersection,
                        'union_cynogs': union,
                    })

    for window_name, _, _ in WINDOWS:
        print(f"    {window_name}: {len(results[window_name])} syntenic pairs (J >= {threshold})")

    return results


def merge_results(results_by_window, loci_df):
    """Merge results across windows and assign confidence tiers."""
    print("\nMerging results and assigning confidence...")

    loci_meta = loci_df.set_index('locus_id').to_dict('index')

    pairs_by_window = {}
    for window_name, pairs in results_by_window.items():
        pairs_by_window[window_name] = {
            (p['locus_1'], p['locus_2']): p for p in pairs
        }

    all_pairs = set()
    for pairs in pairs_by_window.values():
        all_pairs.update(pairs.keys())

    print(f"  Total unique syntenic pairs: {len(all_pairs)}")

    merged = []
    for pair in all_pairs:
        l1, l2 = pair

        in_narrow = pair in pairs_by_window.get('narrow', {})
        in_mid = pair in pairs_by_window.get('mid', {})
        in_wide = pair in pairs_by_window.get('wide', {})

        if in_narrow:
            confidence = 'high'
        elif in_mid:
            confidence = 'medium'
        else:
            confidence = 'low'

        if in_narrow:
            best = pairs_by_window['narrow'][pair]
        elif in_mid:
            best = pairs_by_window['mid'][pair]
        else:
            best = pairs_by_window['wide'][pair]

        meta1 = loci_meta.get(l1, {})
        meta2 = loci_meta.get(l2, {})

        merged.append({
            'locus_1': l1,
            'locus_2': l2,
            'confidence': confidence,
            'in_narrow': in_narrow,
            'in_mid': in_mid,
            'in_wide': in_wide,
            'jaccard': best['jaccard'],
            'shared_cynogs': best['shared_cynogs'],
            'family_1': meta1.get('family', ''),
            'family_2': meta2.get('family', ''),
            'species_1': meta1.get('species', ''),
            'species_2': meta2.get('species', ''),
            'chr_1': meta1.get('chromosome', ''),
            'chr_2': meta2.get('chromosome', ''),
        })

    df = pd.DataFrame(merged)

    if len(df) > 0:
        print(f"  High confidence: {(df['confidence'] == 'high').sum()}")
        print(f"  Medium confidence: {(df['confidence'] == 'medium').sum()}")
        print(f"  Low confidence: {(df['confidence'] == 'low').sum()}")

    return df


def identify_blocks(signatures, loci_df):
    """Identify synteny blocks using average-linkage clustering on narrow Jaccard."""
    print("\nIdentifying synteny blocks (average-linkage clustering)...")

    loci_meta = loci_df.set_index('locus_id').to_dict('index')
    block_rows = []
    block_counter = 0

    for family in TARGET_FAMILIES:
        fam_loci = loci_df[loci_df['family'] == family]
        locus_ids = list(fam_loci['locus_id'])
        n = len(locus_ids)

        if n < 2:
            for lid in locus_ids:
                meta = loci_meta.get(lid, {})
                block_rows.append(_singleton_row(lid, meta))
            continue

        # Build narrow-window Jaccard distance matrix (1 - J)
        dist_matrix = np.ones((n, n))
        np.fill_diagonal(dist_matrix, 0)

        for i, l1 in enumerate(locus_ids):
            s1 = signatures[l1].get('narrow', set())
            if not s1:
                continue
            for j in range(i + 1, n):
                l2 = locus_ids[j]
                s2 = signatures[l2].get('narrow', set())
                if not s2:
                    continue
                inter = len(s1 & s2)
                union = len(s1 | s2)
                if union > 0:
                    jacc = inter / union
                    dist_matrix[i, j] = 1 - jacc
                    dist_matrix[j, i] = 1 - jacc

        # Average-linkage hierarchical clustering
        condensed = squareform(dist_matrix)
        Z = linkage(condensed, method='average')
        labels = fcluster(Z, t=1 - JACCARD_THRESHOLD, criterion='distance')

        clusters = {}
        for lid, label in zip(locus_ids, labels):
            clusters.setdefault(label, []).append(lid)

        multi_blocks = sorted(
            [c for c in clusters.values() if len(c) >= 2],
            key=len, reverse=True
        )
        single_loci = [c[0] for c in clusters.values() if len(c) == 1]

        print(f"  {family}: {len(multi_blocks)} blocks, {len(single_loci)} singletons")

        for block in multi_blocks:
            chrs = [loci_meta[lid]['chromosome'] for lid in block if lid in loci_meta]
            chr_counts = Counter(c for c in chrs if c != 'fragment')
            modal_chr = chr_counts.most_common(1)[0][0] if chr_counts else 'fragment'

            block_id = f'BLOCK_{block_counter:04d}'
            block_counter += 1

            for lid in block:
                meta = loci_meta.get(lid, {})
                block_rows.append({
                    'block_id': block_id,
                    'block_size': len(block),
                    'locus_id': lid,
                    'family': meta.get('family', ''),
                    'species': meta.get('species', ''),
                    'chromosome': meta.get('chromosome', ''),
                    'status': 'in_block',
                    'modal_chr': modal_chr,
                    'n_genes': meta.get('n_genes', 1),
                    'protein_ids': meta.get('protein_ids', ''),
                })

        for lid in single_loci:
            meta = loci_meta.get(lid, {})
            block_rows.append(_singleton_row(lid, meta))

    return pd.DataFrame(block_rows)


def _singleton_row(lid, meta):
    return {
        'block_id': 'SINGLETON',
        'block_size': 1,
        'locus_id': lid,
        'family': meta.get('family', ''),
        'species': meta.get('species', ''),
        'chromosome': meta.get('chromosome', ''),
        'status': 'true_singleton',
        'modal_chr': meta.get('chromosome', ''),
        'n_genes': meta.get('n_genes', 1),
        'protein_ids': meta.get('protein_ids', ''),
    }


def save_outputs(results_by_window, merged_df, blocks_df, signatures, loci_df):
    """Save all outputs."""
    print("\nSaving outputs...")

    for window_name, pairs in results_by_window.items():
        df = pd.DataFrame(pairs)
        if len(df) > 0:
            df = df.sort_values('jaccard', ascending=False)
        path = OUTPUT_DIR / f'synteny_pairs_{window_name}.tsv'
        df.to_csv(path, sep='\t', index=False)
        print(f"  {path.name}: {len(df)} pairs")

    if len(merged_df) > 0:
        # merge_results() traverses a set of pair keys, so confidence/Jaccard
        # ties inherit hash-randomized process order unless the biological pair
        # identifiers are explicit tie-breakers.
        merged_df = merged_df.sort_values(
            ['confidence', 'jaccard', 'locus_1', 'locus_2'],
            ascending=[True, False, True, True],
        )
    path = OUTPUT_DIR / 'synteny_pairs_merged.tsv'
    merged_df.to_csv(path, sep='\t', index=False)
    print(f"  {path.name}: {len(merged_df)} pairs")

    path = OUTPUT_DIR / 'synteny_blocks.tsv'
    blocks_df.to_csv(path, sep='\t', index=False)
    print(f"  {path.name}: {len(blocks_df)} rows")

    # Locus info with CynOG signatures
    loci_meta = loci_df.set_index('locus_id').to_dict('index')
    loci_out = []
    for locus_id, sigs in signatures.items():
        meta = loci_meta.get(locus_id, {})
        loci_out.append({
            'locus_id': locus_id,
            'family': meta.get('family', ''),
            'genome': meta.get('genome', ''),
            'species': meta.get('species', ''),
            'scaffold': meta.get('scaffold', ''),
            'chromosome': meta.get('chromosome', ''),
            'n_genes': meta.get('n_genes', 1),
            'protein_ids': meta.get('protein_ids', ''),
            'n_cynogs_narrow': len(sigs.get('narrow', set())),
            'n_cynogs_mid': len(sigs.get('mid', set())),
            'n_cynogs_wide': len(sigs.get('wide', set())),
            'cynogs_narrow': ','.join(sorted(sigs.get('narrow', set()))),
            'cynogs_mid': ','.join(sorted(sigs.get('mid', set()))),
            'cynogs_wide': ','.join(sorted(sigs.get('wide', set()))),
        })

    loci_out_df = pd.DataFrame(loci_out)
    path = OUTPUT_DIR / 'pcwde_loci.tsv'
    loci_out_df.to_csv(path, sep='\t', index=False)
    print(f"  {path.name}: {len(loci_out_df)} loci")

    # Summary
    print("\n=== SUMMARY ===")
    for fam in TARGET_FAMILIES:
        fam_loci = loci_df[loci_df['family'] == fam]
        fam_blocks = blocks_df[blocks_df['family'] == fam]
        in_block = fam_blocks[fam_blocks['status'] == 'in_block']
        singletons = fam_blocks[fam_blocks['status'] == 'true_singleton']
        n_blocks = in_block['block_id'].nunique()
        print(f"  {fam}: {len(fam_loci)} loci -> {n_blocks} blocks, {len(singletons)} singletons")


def main():
    chd, targets, acc_to_species = load_data()
    chd = inject_reconstructed_order(chd)
    if MISSING_ORDER_FILE:
        # A missing ordered index cannot establish locus grouping or flanks.
        # Retain the accepted coordinate units for downstream placement and
        # explicit NOT_DETERMINED handling, without inflating locus counts.
        missing = targets[targets['position'].isna()].copy()
        missing.to_csv(OUTPUT_DIR / 'missing_order_targets.tsv', sep='\t', index=False)
    gene_index = build_gene_index(chd)
    loci_df = find_loci(gene_index, targets, acc_to_species)
    signatures = build_cynog_signatures(loci_df, gene_index)
    results_by_window = compute_pairwise_synteny(signatures, loci_df, JACCARD_THRESHOLD)
    merged_df = merge_results(results_by_window, loci_df)
    blocks_df = identify_blocks(signatures, loci_df)
    save_outputs(results_by_window, merged_df, blocks_df, signatures, loci_df)
    print("\nDone!")


if __name__ == '__main__':
    main()
