#!/usr/bin/env python3
"""
Summarize and classify synteny blocks per gene family (CynOG version).

Classifications designed to support per-genome evolutionary history reconstruction:

Position classes:
  - core: in main syntenic block (ancestral location)
  - same_chr_conserved: on modal chr, has conserved flanking CynOGs (near ancestral)
  - same_chr_novel: on modal chr, but no conserved CynOGs (TE insertion, coincidental)
  - diff_chr: on different chromosome than modal

Copy status (relative to genome):
  - ancestral: this is a core block member
  - extra_copy: genome also has core member (this is a duplication)
  - only_copy: genome has no core member (relocated or lineage-specific)
"""

import pandas as pd
from pathlib import Path
from collections import Counter

OUTPUT_DIR = Path(__file__).parent / 'outputs'

MIN_CONSERVED_CYNOGS = 1


def load_cynog_signatures(loci_df):
    """Load flanking CynOG signatures from pcwde_loci.tsv (produced by 01)."""
    signatures = {}
    for _, row in loci_df.iterrows():
        cynogs_str = row.get('cynogs_narrow', '')
        cynogs = set(cynogs_str.split(',')) if pd.notna(cynogs_str) and cynogs_str else set()
        cynogs.discard('')
        signatures[row['locus_id']] = cynogs
    return signatures


def jaccard(set1, set2):
    if len(set1) == 0 or len(set2) == 0:
        return 0.0
    return len(set1 & set2) / len(set1 | set2)


def get_conserved_cynogs(core_locus_ids, signatures, min_freq=0.5):
    """Find CynOGs that appear in >= min_freq of core members."""
    cynog_counts = Counter()
    for lid in core_locus_ids:
        for cynog_id in signatures.get(lid, set()):
            cynog_counts[cynog_id] += 1
    n_members = len(core_locus_ids)
    return {c for c, count in cynog_counts.items() if count >= n_members * min_freq}


def process_family(fam_df, signatures):
    """Classify all loci in a family with history-oriented labels."""
    family = fam_df['family'].iloc[0]
    in_block = fam_df[fam_df['status'] == 'in_block']

    if len(in_block) == 0:
        return [{
            'locus_id': row['locus_id'],
            'position_class': 'unknown',
            'copy_status': 'only_copy',
            'genome_has_core': False,
            'core_block': None,
            'modal_chr': None,
            'j_vs_core': None,
            'n_conserved_cynogs': 0,
            'n_cynogs': len(signatures.get(row['locus_id'], set())),
            'conserved_cynogs_in_family': 0,
        } for _, row in fam_df.iterrows()]

    block_sizes = in_block.groupby('block_id').size()
    main_block_id = block_sizes.idxmax()
    main_block = in_block[in_block['block_id'] == main_block_id]
    modal_chr = main_block['modal_chr'].iloc[0]

    # Check for secondary core block (different chr, substantial size)
    secondary_block_id = None
    secondary_chr = None
    if len(block_sizes) > 1:
        second_id = block_sizes.drop(main_block_id).idxmax()
        second_block = in_block[in_block['block_id'] == second_id]
        if len(second_block) >= 5 and second_block['modal_chr'].iloc[0] != modal_chr:
            secondary_block_id = second_id
            secondary_chr = second_block['modal_chr'].iloc[0]

    core_block_ids = {main_block_id}
    modal_chrs = {modal_chr}
    if secondary_block_id:
        core_block_ids.add(secondary_block_id)
        modal_chrs.add(secondary_chr)

    core_locus_ids = list(in_block[in_block['block_id'].isin(core_block_ids)]['locus_id'])
    conserved_cynogs = get_conserved_cynogs(core_locus_ids, signatures, min_freq=0.5)

    core_all_cynogs = set()
    for lid in core_locus_ids:
        core_all_cynogs.update(signatures.get(lid, set()))

    core_genomes = set(in_block[in_block['block_id'].isin(core_block_ids)]['species'].unique())

    results = []
    for _, row in fam_df.iterrows():
        lid = row['locus_id']
        status = row['status']
        chr_ = row['chromosome']
        species = row['species']
        locus_cynogs = signatures.get(lid, set())

        genome_has_core = species in core_genomes
        j_vs_core = jaccard(locus_cynogs, core_all_cynogs)
        n_cynogs = len(locus_cynogs)
        n_conserved = len(locus_cynogs & conserved_cynogs)

        if status == 'in_block' and row['block_id'] in core_block_ids:
            position_class = 'core'
            copy_status = 'ancestral'
            core_block = row['block_id']
        elif chr_ in modal_chrs:
            if n_conserved >= MIN_CONSERVED_CYNOGS:
                position_class = 'same_chr_conserved'
            else:
                position_class = 'same_chr_novel'
            copy_status = 'extra_copy' if genome_has_core else 'only_copy'
            core_block = None
        else:
            position_class = 'diff_chr'
            copy_status = 'extra_copy' if genome_has_core else 'only_copy'
            core_block = None

        results.append({
            'locus_id': lid,
            'position_class': position_class,
            'copy_status': copy_status,
            'genome_has_core': genome_has_core,
            'core_block': core_block,
            'modal_chr': modal_chr + (f'+{secondary_chr}' if secondary_chr else ''),
            'j_vs_core': round(j_vs_core, 3) if position_class != 'core' else None,
            'n_conserved_cynogs': n_conserved,
            'n_cynogs': n_cynogs,
            'conserved_cynogs_in_family': len(conserved_cynogs),
        })

    return results


def main():
    blocks_df = pd.read_csv(OUTPUT_DIR / 'synteny_blocks.tsv', sep='\t')
    loci_df = pd.read_csv(OUTPUT_DIR / 'pcwde_loci.tsv', sep='\t')

    print("Loading CynOG signatures from pcwde_loci.tsv...")
    signatures = load_cynog_signatures(loci_df)
    print(f"  {sum(len(s) > 0 for s in signatures.values())}/{len(signatures)} loci have CynOG data")

    all_results = []
    for family in sorted(blocks_df['family'].unique()):
        fam_df = blocks_df[blocks_df['family'] == family]
        results = process_family(fam_df, signatures)
        for r in results:
            r['family'] = family
        all_results.extend(results)

    class_df = pd.DataFrame(all_results)
    out_df = blocks_df.merge(class_df, on=['family', 'locus_id'])

    cols = ['family', 'locus_id', 'species', 'chromosome', 'position_class',
            'copy_status', 'genome_has_core', 'j_vs_core', 'n_conserved_cynogs', 'n_cynogs',
            'block_id', 'block_size', 'modal_chr_x', 'status', 'n_genes', 'protein_ids']
    cols = [c for c in cols if c in out_df.columns]
    out_df = out_df[cols].rename(columns={'modal_chr_x': 'block_modal_chr'})
    out_df = out_df.sort_values(['family', 'species', 'chromosome', 'locus_id'])

    outfile = OUTPUT_DIR / 'locus_classifications.tsv'
    out_df.to_csv(outfile, sep='\t', index=False)
    print(f"\nWrote {len(out_df)} loci to {outfile}")

    # === BLOCK SUMMARY ===
    block_results = []
    for family in sorted(out_df['family'].unique()):
        fam = out_df[out_df['family'] == family]

        multi = fam[fam['status'] == 'in_block']
        for block_id in sorted(multi['block_id'].unique()):
            block_df = multi[multi['block_id'] == block_id]
            pos_class = block_df['position_class'].iloc[0]

            block_results.append({
                'family': family,
                'block_id': block_id,
                'block_type': 'syntenic_block',
                'position_class': pos_class,
                'size': len(block_df),
                'chromosome': block_df['block_modal_chr'].iloc[0],
                'n_genomes': block_df['species'].nunique(),
                'genomes': '; '.join(sorted(block_df['species'].unique())),
                'loci': '; '.join(sorted(block_df['locus_id'])),
            })

        singletons = fam[fam['status'] == 'true_singleton']
        for _, row in singletons.iterrows():
            block_results.append({
                'family': family,
                'block_id': row['locus_id'],
                'block_type': 'singleton',
                'position_class': row['position_class'],
                'copy_status': row['copy_status'],
                'size': 1,
                'chromosome': row['chromosome'],
                'n_genomes': 1,
                'genomes': row['species'],
                'loci': row['locus_id'],
                'j_vs_core': row['j_vs_core'],
                'n_conserved_cynogs': row['n_conserved_cynogs'],
            })

    block_summary = pd.DataFrame(block_results)
    block_summary = block_summary.sort_values(['family', 'size', 'block_id'], ascending=[True, False, True])

    block_outfile = OUTPUT_DIR / 'block_summary.tsv'
    block_summary.to_csv(block_outfile, sep='\t', index=False)
    print(f"Wrote {len(block_summary)} entries to {block_outfile}")

    # === PER-GENOME SUMMARY ===
    print("\n" + "=" * 70)
    print("PER-GENOME HISTORY SUMMARY")
    print("=" * 70)

    for family in sorted(out_df['family'].unique()):
        fam = out_df[out_df['family'] == family]
        print(f"\n{family}")
        print("-" * 50)

        for species in sorted(fam['species'].unique()):
            sp_df = fam[fam['species'] == species]

            ancestral = sp_df[sp_df['copy_status'] == 'ancestral']
            extra = sp_df[sp_df['copy_status'] == 'extra_copy']
            only = sp_df[sp_df['copy_status'] == 'only_copy']
            same_conserved = sp_df[sp_df['position_class'] == 'same_chr_conserved']

            if len(ancestral) > 0 and len(extra) == 0:
                pattern = "core_only"
            elif len(ancestral) > 0 and len(extra) > 0:
                pattern = "core+dups"
            elif len(only) > 0 and len(same_conserved) > 0:
                pattern = "near_core"
            elif len(only) > 0:
                pattern = "relocated"
            else:
                pattern = "unknown"

            parts = []
            for pc in ['core', 'same_chr_conserved', 'same_chr_novel', 'diff_chr']:
                n = (sp_df['position_class'] == pc).sum()
                if n > 0:
                    parts.append(f"{n} {pc}")

            print(f"  {species}: [{pattern}] {', '.join(parts)}")


if __name__ == '__main__':
    main()
