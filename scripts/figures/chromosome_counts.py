#!/usr/bin/env python3
"""Combined per-chromosome distribution figure for PL1_8, PL4_4, GH5_2.

Cladogram (Hearn-style display backbone) + heatmap panels showing
per-chromosome copy counts. Gem-tone color gradients per family. Species with
gene family members but no chromosome placement shown as light tint + "?".
"""

import re, sys, os, csv
import xml.etree.ElementTree as ET
from pathlib import Path
from collections import defaultdict

# Public repository renderer entry point.

def required_directory(name):
    value = os.environ.get(name, '')
    if not value:
        raise SystemExit(f'{name} is required')
    path = Path(value).resolve()
    if not path.is_dir():
        raise SystemExit(f'{name} is not a directory: {path}')
    return path


INPUT_DIR = required_directory('PCWDE_CHROMOSOME_COUNTS_INPUT_DIR')
SHARED_INPUT_DIR = required_directory('PCWDE_FIGURE_SHARED_INPUT_DIR')
TREE_FILE = INPUT_DIR / 'species_tree.nwk'
ECOLOGY_STATES = Path(__file__).resolve().parents[2] / 'config' / 'ecology_states.tsv'


# ── Families and their NHX files ──
FAMILIES = ['PL1_8', 'PL4_4', 'GH5_2', 'CBM13']

COPY_TABLE = os.environ.get('PCWDE_COPY_TABLE', '')
GENOME_ROSTER = os.environ.get('PCWDE_GENOME_ROSTER', '')
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from copy_table import read_copies, copy_cells


def load_copy_unit_cells(copy_table, roster):
    """Count copies per species and chromosome from the shared copy table."""
    family_map = {'B-lectin': 'CBM13'}
    display_to_tip = {
        'Ceroptres_kovalevi': 'Ceroptres_masudai',
        'Eschatocerus_acaciae': 'Eschatocerus_acacia',
        'Eubothrus_centaureae': 'Isocolus_centaureae',
    }
    with open(roster) as handle:
        species_list = [display_to_tip.get(r['species'].replace(' ', '_'), r['species'].replace(' ', '_'))
                        for r in csv.DictReader(handle, delimiter='\t')]
    if len(species_list) != 49 or len(set(species_list)) != 49:
        raise ValueError(f'Expected a 49-genome roster: {roster}')
    cells = {(display_to_tip.get(sp.replace(' ', '_'), sp.replace(' ', '_')), family_map.get(fam, fam)): counts
             for (sp, fam), counts in copy_cells(read_copies(copy_table)).items()}
    outside = {sp for sp, _ in cells} - set(species_list)
    if outside:
        raise ValueError(f'Copy table species outside the genome roster: {sorted(outside)}')
    result = {family: {} for family in FAMILIES}
    for family in FAMILIES:
        for species in species_list:
            counts = cells.get((species, family), {})
            total = counts.get('total', 0)
            unplaced = counts.get('unplaced', 0)
            attrs = {'total': total, 'placed_total': total - unplaced, 'unplaced': unplaced}
            attrs.update({ch: counts.get(ch, 0) for ch in (f'chr{i}' for i in range(1, 11))})
            if total and not attrs['placed_total']:
                attrs['annotation'] = 'fragment'
            result[family][species] = attrs
    return result

# Fig. 5 methods: show BK-coordinate chromosome placement from BK RagTag
# pseudomolecules for most genomes; allow JCVI-derived chromosome equivalence
# only for the chromosome-scale EP and DK assemblies. protein_to_chromosome.tsv
# does not preserve per-row provenance, so we enforce the display rule here.


# ── Gem-tone family colors (header + gradient endpoint) ──
FAMILY_COLORS = {
    'PL1_8': (33, 102, 172),   # sapphire #2166AC
    'PL4_4': (178, 24, 43),    # garnet   #B2182B
    'GH5_2': (27, 120, 55),    # emerald  #1B7837
    'CBM13': (217, 140, 0),    # amber    #D98C00 (B-lectin domain HGT family)
}


# Empty genome-locus boxes marked with an asterisk indicate transcript/read
# evidence that makes a strict genome-layer absence unresolved.
ASTERISK_ABSENCES = {
    ('PL4_4', 'Ceroptres_masudai'),
    ('CBM13', 'Druon_quercuslanigerum'),
}

# ── Name fixes (NHX uses truncated names) ──

# ── Species to prune from the tree ──
PRUNE_SPECIES = {'Periclistus_JH2016', 'Callirhytis_spRG'}

ECOLOGY_TRACKS = [
    ('gall_associated', 'Gall', '#8A7967'),
    ('herbivorous', 'Herb.', '#4F4F4F'),
]


def load_ecology_states(tip_names):
    """Load the curated states for the exact accepted Figure 2 panel."""
    required = {'species', 'gall_associated', 'herbivorous'}
    states = {}
    with open(ECOLOGY_STATES) as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f'Missing Figure 2 ecology columns: {sorted(missing)}')
        for row in reader:
            species = row['species'].strip()
            if not species or species in states:
                raise ValueError(f'Blank or duplicate ecology species: {species!r}')
            values = {key: row[key].strip() for key, _, _ in ECOLOGY_TRACKS}
            invalid = {key: value for key, value in values.items()
                       if value not in {'yes', 'no', 'uncertain'}}
            if invalid:
                raise ValueError(f'{species}: invalid ecology states {invalid}')
            states[species] = values
    if set(states) != set(tip_names):
        raise ValueError(
            f'Ecology/Figure 2 species mismatch; missing={sorted(set(tip_names) - set(states))}, '
            f'extra={sorted(set(states) - set(tip_names))}'
        )
    return states

# ── Desired display order (top to bottom) ──
DESIRED_ORDER = [
    'Telenomus_remus',
    # Diplolepidae
    'Pediaspis_aceris', 'Diplolepis_rosae', 'Diplolepis_spinosa',
    # Figitidae
    'Parnips_nigripes', 'Euceroptres_primus',
    'Alloxysta_arcuata', 'Phaenoglyphis_villosa',
    'Ganaspis_brasiliensis',
    'Trybliographa_sp', 'Trichoplasta_sp',
    'Leptopilina_boulardi', 'Leptopilina_clavipes',
    'Leptopilina_drosophilae', 'Leptopilina_heterotoma',
    # Eschatocerus
    'Eschatocerus_acacia',
    # Cynipidae — herb gallers first (flipped), oak gallers last
    'Phanacis_lampsanae', 'Timaspis_phoenixopodos',
    'Hedickiana_levantina', 'Fumariphilus_hypecoi',
    'Isocolus_centaureae', 'Aulacidea_hieracii', 'Aulacidea_tavakolii',
    'Qwaqwaia_scolopiae',
    'Synophrus_hispanicus',
    'Synergus_laeviventris', 'Synergus_magnus',
    'Synergus_coniferae', 'Synergus_itoensis',
    'Synergus_gifuensis', 'Synergus_japonicus', 'Synergus_umbraculus',
    'Periclistus_caninae', 'Diastrophus_rubi', 'Diastrophus_kincaidii',
    'Aylax_minor',
    'Ceroptres_masudai',
    'Protobalandricus_spectabilis',
    'Dryocosmus_kuriphilus',
    'Melikaiella_ostensackeni', 'Neuroterus_valhalla',
    'Amphibolips_quercusostensackenii',
    'Andricus_quercuspetiolicola', 'Andricus_quercusramuli', 'Andricus_curvator',
    'Cynips_quercusfolii',
    'Belonocnema_kinseyi', 'Druon_quercuslanigerum',
    'Acraspis_erinacei',
]

TAXON_LABEL_GROUPS = [
    ('Diplolepididae', [
        'Pediaspis_aceris', 'Diplolepis_rosae', 'Diplolepis_spinosa',
    ]),
    ('Figitidae', [
        'Parnips_nigripes', 'Euceroptres_primus',
        'Alloxysta_arcuata', 'Phaenoglyphis_villosa',
        'Ganaspis_brasiliensis', 'Trybliographa_sp', 'Trichoplasta_sp',
        'Leptopilina_boulardi', 'Leptopilina_clavipes',
        'Leptopilina_drosophilae', 'Leptopilina_heterotoma',
    ]),
    ('Cynipidae', [
        'Eschatocerus_acacia', 'Phanacis_lampsanae',
        'Timaspis_phoenixopodos', 'Hedickiana_levantina',
        'Fumariphilus_hypecoi', 'Aulacidea_hieracii',
        'Aulacidea_tavakolii', 'Isocolus_centaureae', 'Qwaqwaia_scolopiae',
        'Synophrus_hispanicus',
        'Synergus_laeviventris', 'Synergus_magnus', 'Synergus_coniferae',
        'Synergus_itoensis', 'Synergus_gifuensis', 'Synergus_japonicus',
        'Synergus_umbraculus', 'Periclistus_caninae', 'Diastrophus_rubi',
        'Diastrophus_kincaidii', 'Aylax_minor', 'Ceroptres_masudai',
        'Protobalandricus_spectabilis', 'Dryocosmus_kuriphilus',
        'Cynips_quercusfolii', 'Melikaiella_ostensackeni', 'Neuroterus_valhalla',
        'Amphibolips_quercusostensackenii', 'Andricus_quercuspetiolicola',
        'Andricus_quercusramuli', 'Andricus_curvator',
        'Belonocnema_kinseyi', 'Druon_quercuslanigerum', 'Acraspis_erinacei',
    ]),
]


def unresolved_display_state(attrs):
    """Return whether an unplaced annotated copy needs question boxes."""
    placed_total = attrs.get('placed_total')
    if placed_total is None:
        placed_total = sum(
            value for key, value in attrs.items()
            if key.startswith('chr') and isinstance(value, int)
        )
    # Any placement makes the numeral counts authoritative for display.
    if int(placed_total or 0) > 0:
        return False, False
    annotation_state = attrs.get('annotation')
    unresolved_count = int(attrs.get('unplaced', 0) or 0)
    show_questions = (
        annotation_state in ('fragment', 'unplaced')
        or unresolved_count > 0
    )
    return show_questions, False


def get_active_chromosomes(annotations):
    """Return chromosomes that have nonzero counts in at least one species."""
    chroms = set()
    for sp, attrs in annotations.items():
        for k, v in attrs.items():
            if k.startswith('chr') and isinstance(v, int) and v > 0:
                chroms.add(k)
    ordered = sorted(chroms, key=lambda c: int(c.replace('chr', '')))
    return ordered


def parse_newick_cladogram(newick_str):
    """Parse newick, return Bio.Phylo tree."""
    from Bio import Phylo
    import io
    return Phylo.read(io.StringIO(newick_str), 'newick')


# Display-only renames (DB-canonical scientific names; upstream data files
# still use the legacy underscore form so cell-data lookups continue to work).
DISPLAY_EPITHET_RENAME = {
    ('Eschatocerus', 'acacia'): 'acaciae',
}
DISPLAY_BINOMIAL_RENAME = {
    ('Isocolus', 'centaureae'): ('Eubothrus', 'centaureae'),
    ('Ceroptres', 'masudai'): ('Ceroptres', 'kovalevi'),
    ('Bradysia', 'coprophila'): ('Bradysia', 'tilicola'),
}

def short_label(name, abbrevs):
    """Full public species label; *sp.* styling is handled by the SVG renderer."""
    parts = name.split('_')
    if len(parts) >= 2:
        genus = parts[0]
        epithet = ' '.join(parts[1:])
        genus, epithet = DISPLAY_BINOMIAL_RENAME.get((genus, epithet), (genus, epithet))
        epithet = DISPLAY_EPITHET_RENAME.get((genus, epithet), epithet)
        return f"{genus} {epithet}." if epithet == 'sp' else f"{genus} {epithet}"
    return name


def reorder_tree(tree, desired_order):
    """Sort children at each internal node to match desired tip order."""
    order_map = {name: i for i, name in enumerate(desired_order)}

    def min_order(clade):
        if clade.is_terminal():
            return order_map.get(clade.name, 999)
        return min(min_order(c) for c in clade.clades)

    for clade in tree.find_clades(order='postorder'):
        if not clade.is_terminal():
            clade.clades.sort(key=min_order)


def build_display_tree(desired_order):
    """Build the rooted/pruned Figure 5 tree and apply display-only ordering."""
    tree = parse_newick_cladogram(TREE_FILE.read_text().strip())
    tree.root_with_outgroup('Telenomus_remus')
    present_tips = {tip.name for tip in tree.get_terminals()}
    for sp in PRUNE_SPECIES:
        if sp in present_tips:
            tree.prune(sp)
    reorder_tree(tree, desired_order)
    return tree


def lerp_color(intensity, family):
    """Gradient from light gray to family color."""
    hi = FAMILY_COLORS.get(family, (100, 100, 100))
    lo = (235, 235, 235)
    r = int(lo[0] + (hi[0] - lo[0]) * intensity)
    g = int(lo[1] + (hi[1] - lo[1]) * intensity)
    b = int(lo[2] + (hi[2] - lo[2]) * intensity)
    return f'#{r:02x}{g:02x}{b:02x}'


