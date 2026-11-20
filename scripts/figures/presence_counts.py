#!/usr/bin/env python3
"""PCWDE gene count heatmap on a cladogram — SVG output.

Reads a Hearn-style display cladogram, embeds per-tip counts as NHX, then renders heatmap SVG.
"""

import csv
import re, sys, os
import io
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

from Bio import Phylo

# Public repository renderer entry point.

def required_directory(name):
    value = os.environ.get(name, '')
    if not value:
        raise SystemExit(f'{name} is required')
    path = Path(value).resolve()
    if not path.is_dir():
        raise SystemExit(f'{name} is not a directory: {path}')
    return path


INPUT_DIR = required_directory('PCWDE_PRESENCE_COUNTS_INPUT_DIR')
TREE_FILE = INPUT_DIR / 'cazy_scan_tree.nwk'
HITS_FILE = INPUT_DIR / 'cazy_hits_all.tsv'


# Curated B_lectin/CBM13 HGT family membership (51 proteins, 3 tribes)
BLECTIN_INTERFACE = INPUT_DIR / 'beta_units.tsv'
NEW_GENOME_COUNTS = INPUT_DIR / 'ordinary_family_current_panel_additions.tsv'
# GH1_e72 is supplied by this build's derived CAZy count table in tree_counts.py.


def load_blectin_counts():
    """Per-species count of curated B_lectin/CBM13 HGT family proteins.

    The raw CBM13 column in cazy_hits_all.tsv over-counts: GT27
    N-acetylgalactosaminyltransferases carry Ricin_B_lectin CBM13-like
    domains unrelated to the HGT family. We instead use the 51-protein
    curated set from the blectin_lgt_investigation goal.
    """
    counts = defaultdict(int)
    if not BLECTIN_INTERFACE.is_file():
        raise FileNotFoundError(f'Missing required uniform beta interface: {BLECTIN_INTERFACE}')
    import csv
    with open(BLECTIN_INTERFACE) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            species = row['species'].strip().replace(' ', '_')
            if species:
                counts[species] += 1
    return counts


def load_counts():
    """Load per-genome counts for the four HGT targets.

    - PL1_8, PL4_4, GH5_2: subfamily-level matches from the broad dbCAN scan.
    - CBM13: substituted with curated B_lectin HGT membership
      (see load_blectin_counts); raw CBM13 hits over-count due to GT27
      N-acetylgalactosaminyltransferase contamination.
    - GH1_e72: unique gene/locus counts from uniform competitive full-dbCAN-sub
      v5.2.4 runs across all 27 displayed proteomes.
    """
    counts = defaultdict(lambda: defaultdict(int))
    subfam_targets = {'PL1_8', 'PL4_4', 'GH5_2'}
    with open(HITS_FILE) as f:
        header = f.readline().strip().split('	')
        for line in f:
            fields = line.strip().split('	')
            row = dict(zip(header, fields))
            genome = row['genome']
            subfam = row['cazy_subfamily']
            if subfam in subfam_targets:
                counts[genome][subfam] += 1

    blectin = load_blectin_counts()
    for species, n in blectin.items():
        counts[species]['CBM13'] = n

    # Accepted annotation-backed counts for the nine post-CynOGv3 genomes.
    # These replace missing broad-scan cells. GH1_e72 is populated separately
    # from each genome's eligible canonical annotation interface above.
    with open(NEW_GENOME_COUNTS) as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            species = row['species']
            counts[species][row['family']] = int(row['count'])
    return counts


