#!/usr/bin/env python3
from pathlib import Path
import argparse, csv, hashlib, os, re, subprocess, shutil, sys, json
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'release'))
import author_captions  # noqa: E402
AUTHOR_CAPTIONS = author_captions.load()

parser = argparse.ArgumentParser()
parser.add_argument("--base", type=Path, required=True)
parser.add_argument("--out", type=Path, required=True)
parser.add_argument("--built", type=Path, action="append", required=True,
                    help="directory of tables built by a producer; its Table_S{n} files are used instead of --base")
parser.add_argument("--s2a-only", action="store_true", help="render only the expression sample-count panel")
parser.add_argument("--s7", type=Path, required=False,
                    help="build_table_s7.py output: selection_reader_table.tsv and catalytic_reader_table.tsv")
args = parser.parse_args()
BASE = args.base.resolve()
OUT = args.out.resolve()
BUILT = [d.resolve() for d in args.built]
S7 = args.s7.resolve() if args.s7 else None
if not args.s2a_only and S7 is None:
    parser.error("--s7 is required for the full table render")
with (Path(__file__).resolve().parents[2] / 'config/table_numbers.tsv').open(newline='') as handle:
    TABLE_LABELS={r['table_id']:r['label'] for r in csv.DictReader(handle,delimiter='\t')}
SOURCE_IDS={1:'catalog',2:'outgroup',3:'loci',4:'synteny',5:'intron',6:'expression',7:'selection_catalytic',14:'cazy_rank',15:'chromosome_gain_summary',16:'gc3'}
if set(TABLE_LABELS.values())!={f'S{i}' for i in range(1,11)}:raise SystemExit('table numbering map incomplete')
FIG_OLD_IDS={1:'origin_trees_full',2:'accessory_domain_trees',3:'complete_locus_trees',4:'b_lectin_domain_tree',5:'euceroptres_genome'}
with (Path(__file__).resolve().parents[2] / 'config/figure_numbers.tsv').open(newline='') as handle:
    FIG_LABELS={r['figure_id']:r['label'] for r in csv.DictReader(handle,delimiter='\t')}
def remap_refs(text):
    def one(m):
        kind,old=m.group(1),int(m.group(2))
        if kind=='Table':
            stable=SOURCE_IDS.get(old)
            return f'Table {TABLE_LABELS[stable]}' + ('b' if stable == 'expression' else '') if stable else m.group(0)
        stable=FIG_OLD_IDS.get(old)
        return f'Figure {FIG_LABELS[stable]}' if stable else m.group(0)
    return re.sub(r'\b(Table|Figure) S(\d+)\b',one,text)

def table_dir(n):
    """Directory holding Table_S{n}.tsv: a producer's build output, else the deposit."""
    hits = [d for d in BUILT if (d / f"Table_S{n}.tsv").is_file()]
    if len(hits) > 1:
        raise SystemExit(f"Table S{n} built by more than one producer")
    return hits[0] if hits else BASE
if not BASE.is_dir():
    raise SystemExit(f"missing declared table input directory: {BASE}")
if OUT.exists():
    raise SystemExit(f"refusing to overwrite output directory: {OUT}")
OUT.mkdir(parents=True, exist_ok=True)
ORDER = {"GH5_2":0, "PL1_8":1, "PL4_4":2, "B_lectin":3}
with (Path(__file__).resolve().parents[2] / 'config/cazy_lineage_groups.tsv').open(newline='') as handle:
    SPECIES_NAMES = {r['species']: r['source_panel_species'].replace('_', ' ')
                     for r in csv.DictReader(handle, delimiter='\t')}
SPECIES_NAMES.update({'Isocolus_centaureae': 'Eubothrus centaureae',
                      'Eschatocerus_acacia': 'Eschatocerus acaciae'})
SPECIES_COLUMNS = {'species', 'species_display', 'source_species'}
ID_COLUMNS = {'species_id', 'species_code', 'internal_species_code'}
FONT_HASHES = {
    "texgyretermes-regular.otf": "cc3fe7c707b81428d23d54df3eadd9228a2bf6a4d43125d94df56f5f63134659",
    "texgyretermes-bold.otf": "2fb3e952065fa153c7e4e64e04b98b9d79225739b6025aa3f0f0782d299ff61e",
    "texgyretermes-italic.otf": "6dd103a1672e50568cd2f8a706ccd48443d44d7d073a59d2286f4e6f746575d6",
    "texgyretermes-bolditalic.otf": "1bf6af99cb0e26c12951317032d79b96ae009551e59ccf02a5b24f325ecfec87",
}

def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def resolve_font_dir():
    override = os.environ.get("PCWDE_TEX_GYRE_DIR")
    if override:
        font_dir = Path(override)
    else:
        proc = subprocess.run(["kpsewhich", "texgyretermes-regular.otf"], text=True,
                              stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        if proc.returncode or not proc.stdout.strip():
            raise SystemExit("TeX Gyre Termes font resolution failed closed")
        font_dir = Path(proc.stdout.strip()).parent
    for filename, expected in FONT_HASHES.items():
        path = font_dir / filename
        if not path.is_file() or sha256(path) != expected:
            raise SystemExit(f"TeX Gyre Termes font identity failed closed: {filename}")
    return font_dir.resolve()

FONT_DIR = resolve_font_dir()
FONT_DECL = (rf"\setmainfont[Path={FONT_DIR.as_posix()}/,Extension=.otf,"
             r"UprightFont=*-regular,BoldFont=*-bold,ItalicFont=*-italic,"
             r"BoldItalicFont=*-bolditalic]{texgyretermes}")

def esc(value):
    s = str(value)
    replacements = [("\\",r"\textbackslash{}"),("&",r"\&"),("%",r"\%"),("$",r"\$"),("#",r"\#"),("_",r"\_"),("{",r"\{"),("}",r"\}"),("~",r"\textasciitilde{}"),("^",r"\textasciicircum{}")]
    for old,new in replacements: s=s.replace(old,new)
    return s

if not args.s2a_only:
    cazy_rank_source = table_dir(14) / "Table_S14_precision.tsv"
    if not cazy_rank_source.is_file():
        raise SystemExit("missing same-build full-precision CAZy lineage intermediate")
    with cazy_rank_source.open(newline='') as handle:
        precision_rows = list(csv.DictReader(handle, delimiter='\t'))
    cazy_full_percentiles = {(r['lineage_group'], r['family']): r['median_midrank_percentile_full_precision']
                             for r in precision_rows}
    with (table_dir(14) / 'Table_S14.tsv').open(newline='') as handle:
        expected = {(r['lineage_group'], r['family']) for r in csv.DictReader(handle, delimiter='\t')}
    if len(precision_rows) != 15 or len(cazy_full_percentiles) != 15 or set(cazy_full_percentiles) != expected:
        raise SystemExit("missing or duplicate CAZy lineage precision join")
    shutil.copy2(table_dir(4)/"Table_S4.tsv", OUT/f"Table_{TABLE_LABELS['synteny']}.tsv")

    # Table S7 combines the two current scientific summaries into the same
    # section-aware reader-table interface used by the other supplementary tables.
    with (S7 / "selection_reader_table.tsv").open(newline="") as handle:
        selection_rows = list(csv.DictReader(handle, delimiter="\t"))
    with (S7 / "catalytic_reader_table.tsv").open(newline="") as handle:
        catalytic_rows = list(csv.DictReader(handle, delimiter="\t"))
    s7_fields = [
        "section", "family", "analyzed_sequences", "codon_sites",
        "purifying_sites_n_pct", "diversifying_sites_n_pct", "median_omega",
        "significant_tested_terminal_branches_n_N_pct", "site_role",
        "reference_residue", "observable_total", "expected_observable",
        "expected_among_observable_pct",
    ]
    s7_rows = []
    for row in selection_rows:
        s7_rows.append({
            "section": "A", "family": row["gene_family"],
            "analyzed_sequences": row["analyzed_sequences"],
            "codon_sites": row["retained_codon_sites"],
            "purifying_sites_n_pct": row["purifying_display"],
            "diversifying_sites_n_pct": row["diversifying_display"],
            "median_omega": row["median_omega_display"],
            "significant_tested_terminal_branches_n_N_pct": row["terminal_display"],
            "site_role": "", "reference_residue": "", "observable_total": "",
            "expected_observable": "", "expected_among_observable_pct": "",
        })
    for row in catalytic_rows:
        s7_rows.append({
            "section": "B", "family": row["gene_family"],
            "analyzed_sequences": "", "codon_sites": "",
            "purifying_sites_n_pct": "", "diversifying_sites_n_pct": "",
            "median_omega": "", "significant_tested_terminal_branches_n_N_pct": "",
            "site_role": row["site_role"], "reference_residue": row["reference_residue"],
            "observable_total": row["observable_total_display"],
            "expected_observable": row["expected_observable_display"],
            "expected_among_observable_pct": row["expected_observable_pct"],
        })
    s10_labels = {
        'family': 'Gene family', 'analyzed_sequences': 'Analyzed codon alignment sequences (n)',
        'codon_sites': 'Analyzed codon alignment sites (n)',
        'purifying_sites_n_pct': 'Sites under purifying selection (FUBAR), n (%)',
        'diversifying_sites_n_pct': 'Sites under positive selection (FUBAR), n (%)',
        'median_omega': 'Median dN/dS across sites',
        'significant_tested_terminal_branches_n_N_pct': 'Terminal branches with positive selection (aBSREL), n/N (%)',
        'site_role': 'Catalytic site role', 'reference_residue': 'Reference residue',
        'observable_total': 'Scorable / total sequences',
        'expected_observable': 'Conserved / scorable sequences',
        'expected_among_observable_pct': 'Catalytic-residue conservation (%)',
    }
    s10_fields = [s10_labels.get(field, field) for field in s7_fields]
    with (OUT / f"Table_{TABLE_LABELS['selection_catalytic']}.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=s10_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({s10_labels.get(field, field): row[field] for field in s7_fields} for row in s7_rows)

    for n in SOURCE_IDS:
        tdir=table_dir(n)
        src=tdir/f"Table_S{n}.tsv"
        table_label=TABLE_LABELS[SOURCE_IDS[n]] + ('b' if n == 6 else '')
        dst=OUT/f"Table_{table_label}.tsv"
        if n not in (4, 7): shutil.copy2(src,dst)
        if n == 2:
            with src.open(newline='') as handle:
                source_rows = list(csv.DictReader(handle, delimiter='\t'))
            fields = ['Gene family', 'Cynipoid proteins', 'Outgroup proteins',
                      'Median trimmed domain length (aa)', 'Genomes with occurrence',
                      'Outgroup screen hits']
            data_rows = [{
                fields[0]: r['family_id'], fields[1]: r['cynipoid_proteins'],
                fields[2]: r['outgroup_proteins'],
                fields[3]: (format(float(r['median_trimmed_length_aa']), 'g')
                            if r['median_trimmed_length_aa'] else ''),
                fields[4]: r['genomes_with_occurrence'],
                fields[5]: f"bacteria {r['screen_bacteria_proteins']}; eukaryotes {r['screen_eukaryote_proteins']}",
            } for r in source_rows]
            with dst.open('w', newline='') as handle:
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
                writer.writeheader(); writer.writerows(data_rows)
        if n == 4:
            with src.open(newline='') as handle:
                reader = csv.DictReader(handle, delimiter='\t')
                fields = [field for field in reader.fieldnames if field != 'target_positive_patterns']
                data_rows = list(reader)
            with dst.open('w', newline='') as handle:
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
                writer.writeheader()
                writer.writerows({field: row[field] for field in fields} for row in data_rows)
        if n == 5:
            with src.open(newline='') as handle:
                intron_rows = list(csv.DictReader(handle, delimiter='\t'))
            labels = {
                'observation_count': 'copies with this intron (n)',
                'observation_proportion_pct': 'copies with this intron (%)',
                'taxonomic_families': 'cynipoid families with this intron',
            }
            fields = [f for f in intron_rows[0] if f != 'taxonomic_breadth']
            with dst.open('w', newline='') as handle:
                writer = csv.DictWriter(handle, fieldnames=[labels.get(f, f) for f in fields], delimiter='\t', lineterminator='\n')
                writer.writeheader()
                for row in intron_rows:
                    shown = dict(row)
                    for field in ('protein_prevalence_pct', 'observation_proportion_pct'):
                        if shown[field]:
                            shown[field] = f"{float(shown[field]):.1f}".rstrip('0').rstrip('.')
                    writer.writerow({labels.get(f, f): shown[f] for f in fields})
        if n == 16:
            # Reader projection of the unrounded technical GC3 table: headings and
            # display rounding only (checked by check_intron_products.py).
            with src.open(newline='') as handle:
                gc3_rows = list(csv.DictReader(handle, delimiter='\t'))
            gc3_display = [
                ('cohort', 'Gene family', str),
                ('stratum', 'Copy stratum', str),
                ('copies', 'Copies (n)', lambda v: str(int(v))),
                ('species', 'Species (n)', lambda v: str(int(v))),
                ('median_host_gc3_percentile', 'Median transferred-gene percentile in cynipoid genomic background', lambda v: f'{float(v):.1f}'),
                ('species_above_host_median', 'Species with median family-copy percentile above 50 (n)', lambda v: str(int(v))),
                ('median_copy_host_distance_pp', 'Median absolute transferred-gene difference from genomic background (pp)', lambda v: f'{float(v):.1f}'),
                ('median_bacterial_host_distance_pp', 'Median absolute bacterial difference from genomic background (pp)', lambda v: f'{float(v):.1f}'),
                ('paired_wilcoxon_p', 'Paired Wilcoxon P (two-sided; species units)', lambda v: f'{float(v):.1e}'),
                ('median_copy_signed_offset_pp', 'Median signed transferred-gene offset from genomic background (pp)', lambda v: f'{float(v):+.1f}'),
                ('median_bacterial_signed_offset_pp', 'Median signed bacterial offset from genomic background (pp)', lambda v: f'{float(v):+.1f}'),
            ]
            if not gc3_rows or list(gc3_rows[0]) != [c for c, _, _ in gc3_display]:
                raise SystemExit('GC3 technical table schema differs from its reader projection')
            with dst.open('w', newline='') as handle:
                writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
                writer.writerow([label for _, label, _ in gc3_display])
                writer.writerows([fmt(row[c]) for c, _, fmt in gc3_display] for row in gc3_rows)
        if n == 15:
            with src.open(newline='') as handle:
                gain_rows = list(csv.DictReader(handle, delimiter='\t'))
            if len(gain_rows) != 8 or {r['mode'] for r in gain_rows} != {'deltran'}:
                raise SystemExit('Table S7 must contain DELTRAN only')
            fields = [field for field in gain_rows[0] if field != 'mode']
            for row in gain_rows:
                if row['section'] == 'B':
                    categories = {'1': '1', '(1,2]': '>1–2', '(2,4]': '>2–4', '>4': '>4'}
                    row['family_or_stratum'] = categories[row['family_or_stratum']]
                    row['zero_gain_share'] = f"{100 * float(row['zero_gain_share']):.1f}".rstrip('0').rstrip('.')
                    row['gains'] = f"{float(row['gains']):g}"
                if row['section'] == 'A' and row['own_stratum_percentile']:
                    row['own_stratum_percentile'] = f"{float(row['own_stratum_percentile']):g}"
            with dst.open('w', newline='') as handle:
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
                writer.writeheader(); writer.writerows({field: row[field] for field in fields} for row in gain_rows)
        if n != 1:
            with dst.open(newline='') as handle:
                reader = csv.DictReader(handle, delimiter='\t')
                fields = [field for field in reader.fieldnames if field not in ID_COLUMNS]
                data_rows = list(reader)
            for row in data_rows:
                for column in SPECIES_COLUMNS & set(fields):
                    raw = row[column]
                    row[column] = SPECIES_NAMES.get(raw, raw.replace('_', ' '))
            with dst.open('w', newline='') as handle:
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
                writer.writeheader(); writer.writerows({field: row[field] for field in fields} for row in data_rows)
        if n == 1:
            # Use the dedicated Table S1 design and its content-driven pagination.
            env = {**os.environ, "PCWDE_TABLE_S1_INPUT_DIR": str(tdir),
                   "PCWDE_BUILD_ROOT": str(OUT)}
            proc = subprocess.run(
                [sys.executable, "-B", str(Path(__file__).with_name("build_table_s1.py")),
                 "--output", str(OUT / "Table_S1.pdf")],
                env=env, check=True, text=True, stdout=subprocess.PIPE,
            )
            report = json.loads(proc.stdout)
            if report["rows"] != 99 or sum(report["rows_per_page"]) != 99:
                raise SystemExit("Table S1 lost its exact row coverage")
            (OUT / "Table_S1.render.json").write_text(proc.stdout)
            continue
        with dst.open() as h: data=list(csv.DictReader(h,delimiter="\t")); headers=list(data[0])
        key="section" if "section" in headers else ("part" if "part" in headers else None)
        groups=[]
        if key:
            for value in dict.fromkeys(row[key] for row in data): groups.append((value,[row for row in data if row[key]==value]))
        else: groups=[("",data)]
        cap=tdir/f"Table_S{n}_caption.txt"
        caption=cap.read_text().strip() if cap.exists() else f"Table S{n}."
        caption=re.sub(r"^Table S\d+[ab]?\.\s*", "", caption)
        if tdir == BASE: caption=remap_refs(caption)
        author_caption = AUTHOR_CAPTIONS.get(f"table:{SOURCE_IDS.get(n)}")
        landscape=len(headers)>7
        lines=[r"\documentclass[10pt]{article}",r"\usepackage[letterpaper,margin=0.55in]{geometry}",r"\usepackage{fontspec,booktabs,longtable,array,ragged2e,fancyhdr,pdflscape}",FONT_DECL,r"\renewcommand{\arraystretch}{1.08}",r"\setlength{\tabcolsep}{2pt}",r"\setlength{\LTpre}{4pt}",r"\setlength{\LTpost}{4pt}",r"\pagestyle{fancy}",r"\fancyhf{}",rf"\fancyfoot[C]{{\footnotesize Table {table_label} --- Page \thepage}}",r"\newlength{\pcwdeFamilyWidth}",r"\begin{document}"]
        if landscape:
            lines.append(r"\begin{landscape}")
            if len(data) <= 10:
                lines.append(r"\vspace*{0.35in}")
        caption_tex = author_captions.to_tex(author_caption) if author_caption else esc(caption)
        if table_label == TABLE_LABELS['loci']:
            # Species name italics are formatting, applied whether or not the author text marks them.
            caption_tex = re.sub(r'(?:\*|(?<!textit\{))B\. kinseyi\*?', r'\\textit{B. kinseyi}', caption_tex)
        lines += [rf"\section*{{Table {table_label}}}",rf"\small {caption_tex}\par\medskip"]
        for group_label,grp in groups:
            cols=[c for c in headers if c != key and any(row[c] for row in grp)
                  and not (table_label == TABLE_LABELS['synteny'] and c == 'target_positive_patterns')]
            if group_label: lines.append(rf"\textbf{{Part {esc(group_label)}.}}\par\smallskip")
            width=.95/max(1,len(cols)); spec="@{}"+"".join(rf">{{\RaggedRight\arraybackslash}}p{{{width:.3f}\linewidth}}" for _ in cols)+"@{}"
            if table_label == TABLE_LABELS['intron'] and group_label == 'A' and 'phase' in cols:
                widths = {c: (.090 if c == 'phase' else .860 / (len(cols) - 1)) for c in cols}
                spec = "@{}" + "".join(rf">{{\RaggedRight\arraybackslash}}p{{{widths[c]:.3f}\linewidth}}" for c in cols) + "@{}"
            if table_label == TABLE_LABELS['intron'] and group_label == 'B':
                # Short numeric columns donate space to the taxonomic-family list.
                weights = {'reference_position_aa': .055, 'phase': .090,
                           'copies with this intron (n)': .055, 'copies with this intron (%)': .060,
                           'cynipoid families with this intron': .33}
                remaining = .95 - sum(weights.get(c, 0) for c in cols)
                other = [c for c in cols if c not in weights]
                widths = {c: weights.get(c, remaining / max(1, len(other))) for c in cols}
                spec = "@{}" + "".join(rf">{{\RaggedRight\arraybackslash}}p{{{widths[c]:.3f}\linewidth}}" for c in cols) + "@{}"
                # Capture the full table width before longtable changes linewidth per cell.
                lines.append(rf"\setlength{{\pcwdeFamilyWidth}}{{{widths['cynipoid families with this intron']:.3f}\linewidth}}")
            synteny_labels = {
                'family_id': 'Family', 'eligible_target_genes': 'Accepted gene copies',
                'detectable_patterns': 'Synteny blocks',
                'cross_genome_patterns': 'Blocks found in ≥2 genomes',
                'genes_accounted_for': 'Copies assigned to blocks',
                'modal_locus_id': 'Reference chromosome',
                'detected_genomes': 'Genomes with region recovered',
                'occupied_genomes': 'Genomes with target present',
                'empty_genomes': 'Genomes with target absent from recovered region',
                'not_detected_genomes': 'Genomes with region not recovered',
                'blocks_one_genome': 'Blocks found in 1 genome',
                'blocks_2_to_24_genomes': 'Blocks found in 2–24 genomes',
                'blocks_25_plus_genomes': 'Blocks found in ≥25 genomes',
            }
            gain_labels_a = {
                'family_or_stratum': 'Gene family', 'gains': 'Inferred chromosome gains',
                'cynipini_gains': 'Gains in Cynipini', 'other_lineage_gains': 'Gains in other cynipoids',
                'losses': 'Inferred chromosome losses', 'placed_tips_or_native_groups': 'Genomes with chromosome placements',
                'gains_per_placed_tip': 'Gains per genome',
                'all_native_percentile': 'Percentile among all native groups',
                'own_stratum_percentile': 'Percentile among groups with similar copy numbers',
            }
            gain_labels_b = {
                'family_or_stratum': 'Copy-number category (median copies per genome)',
                'gains': 'Median inferred gains', 'placed_tips_or_native_groups': 'Native groups compared',
                'gains_per_placed_tip': 'Median gains per genome',
                'zero_gain_share': 'Groups with no inferred gains (%)',
            }
            labels = synteny_labels if table_label == TABLE_LABELS['synteny'] else (gain_labels_a if group_label == 'A' else gain_labels_b) if table_label == TABLE_LABELS['chromosome_gain_summary'] else {}
            reader_labels = {
                'family_id': 'Gene family', 'median_tpm': 'Median TPM',
                'iqr_tpm': 'TPM interquartile range', 'range_tpm': 'TPM range',
                'species_display': 'Species', 'accepted_target_count': 'Accepted target copies',
                'lineage_group': 'Lineage', 'genomes_n': 'Genomes (n)',
                'protein_prevalence_pct': 'Proteins with introns (%)',
                'reference_position_aa': 'Reference position (aa)',
            }
            if n in (3, 6, 14, 5):
                labels = {**{c: c.replace('_', ' ').capitalize() for c in cols}, **reader_labels, **labels}
            header=" & ".join(rf"\textbf{{{esc(labels.get(c, c.replace('_',' ')))}}}" for c in cols)+r" \\"
            lines += [r"\scriptsize",rf"\begin{{longtable}}{{{spec}}}",r"\toprule",header,r"\midrule\endfirsthead",r"\toprule",header,r"\midrule\endhead"]
            for row in grp:
                values=[]
                for c in cols:
                    raw=remap_refs(row[c]) if tdir == BASE else row[c]
                    # S9-S13 sources carry binary float noise (e.g. 99.08250000000001);
                    # display at most 6 decimals. The TSV keeps the source value.
                    if n >= 9 and re.fullmatch(r"-?\d+\.\d{7,}", raw):
                        raw=f"{round(float(raw),6):.6f}".rstrip("0").rstrip(".")
                    if n in (2, 3, 5, 6, 14, 15):
                        raw = {'B_lectin': 'B-lectin', 'PI_PLC_X': 'PI-PLC-X'}.get(raw, raw)
                    if n == 6 and c not in SPECIES_COLUMNS and not re.search(r'\d', raw):
                        raw = raw.replace('_', ' ')
                    if n == 14 and 'percentile' in c and raw:
                        raw = f"{float(cazy_full_percentiles[row['lineage_group'], row['family']]):.1f}"
                    if n == 15 and c == 'gains_per_placed_tip' and raw:
                        raw = f"{float(raw):.2f}"
                    value=esc(raw)
                    if table_label == TABLE_LABELS['intron'] and group_label == 'B' and c == 'cynipoid families with this intron':
                        lines.append(rf"\settowidth{{\dimen0}}{{{value}}}\ifdim\dimen0>\pcwdeFamilyWidth\PackageError{{pcwde}}{{Intron family cell exceeds column width}}{{}}\fi")
                        value = rf"\mbox{{{value}}}"

                    if (c in SPECIES_COLUMNS or c == 'lineage_group') and re.fullmatch(r'(?:[A-Z][a-z]+|[A-Z]\.) [a-z][a-z.-]+', raw):
                        value=rf'\textit{{{value}}}'
                    if n >= 9:
                        # Long identifiers in S9-S13 (species, genome, reason codes) wrap at underscores.
                        value=value.replace(r"\_", r"\_\allowbreak{}")
                    values.append(value)
                lines.append(" & ".join(values)+r" \\")
            lines += [r"\bottomrule",r"\end{longtable}",r"\normalsize"]
        if landscape:
            lines.extend([r"\clearpage", r"\end{landscape}"])
        lines.append(r"\end{document}")
        tex=OUT/f"Table_{table_label}.tex"; tex.write_text("\n".join(lines)+"\n")
        for _ in range(2):
            proc=subprocess.run(["xelatex","-interaction=nonstopmode","-halt-on-error",tex.name],cwd=OUT,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
            (OUT/f"Table_{table_label}.xelatex.log").write_text(proc.stdout)
            if proc.returncode: raise SystemExit(f"Table S{n} XeLaTeX failed")

# The expression count panel shares the same producer and typeface as S2b.
matrix_source = next((d.parent / "data/expression/sample_count_matrix.tsv" for d in BUILT
                      if (d.parent / "data/expression/sample_count_matrix.tsv").is_file()), None)
if matrix_source is None:
    raise SystemExit("missing same-build expression sample-count matrix")
shutil.copy2(matrix_source, OUT / "Table_S2a.tsv")
with matrix_source.open(newline="") as handle:
    matrix = list(csv.DictReader(handle, delimiter="\t"))
if len(matrix) != 15 or sum(bool(v) for row in matrix for k, v in row.items() if k not in ("species", "source")) != 37:
    raise SystemExit("sample-count panel dimensions changed")
short = ["AG", "VG", "Car", "Whole", "Gut", "SG", "MT", "Whole"]
stage = ["Adult"] * 4 + ["Larva"] * 5
# Letter text width is 8.5 - 2(0.55) = 7.4 in. With @{} outer edges,
# eleven columns contribute twenty 2-pt tabcolsep spaces (40 pt total).
species_width, source_width, context_width = 1.60, 0.78, 0.46
panel_width = species_width + source_width + 8 * context_width + 40 / 72.27
if 7.4 - panel_width < 0.25:
    raise SystemExit(f"Table S2a portrait width margin too small: {7.4 - panel_width:.3f} in")
head = [r"\textbf{Species}", r"\shortstack[l]{\textbf{Source}\\\textit{study /}\\\textit{BioProject}}"] + [rf"\shortstack{{\textbf{{{esc(s)}}}\\{esc(t)}}}" for s, t in zip(stage, short)]
lines = [r"\documentclass[10pt]{article}", r"\usepackage[letterpaper,margin=0.55in]{geometry}",
         r"\usepackage{fontspec,booktabs,array,ragged2e,fancyhdr}", FONT_DECL,
         r"\setlength{\tabcolsep}{2pt}", r"\renewcommand{\arraystretch}{1.22}",
         r"\pagestyle{fancy}", r"\fancyhf{}", r"\fancyfoot[C]{\footnotesize Table S2a --- Page \thepage}",
         r"\begin{document}", r"\section*{Table S2a}",
         r"\small Cells give counts of retained expression profiles by species and context; AG, accessory gland; VG, venom gland; Car, carcass; SG, salivary gland; MT, Malpighian tubule; Whole, whole body; larval Whole includes whole larvae collected from galls.\par\medskip",
         r"\scriptsize", rf"\begin{{tabular}}{{@{{}}>{{\RaggedRight\arraybackslash}}p{{{species_width:.2f}in}}>{{\RaggedRight\arraybackslash}}p{{{source_width:.2f}in}}*{{8}}{{>{{\centering\arraybackslash}}p{{{context_width:.2f}in}}}}@{{}}}}",
         r"\toprule", " & ".join(head) + r" \\", r"\midrule"]
for row in matrix:
    vals = [row["species"], row["source"]] + [row[k] for k in list(row)[2:]]
    cells = [r"\textit{" + esc(vals[0]) + "}"] + [esc(v) for v in vals[1:]]
    lines.append(" & ".join(cells) + r" \\")
lines += [r"\bottomrule", r"\end{tabular}", r"\end{document}"]
tex = OUT / "Table_S2a.tex"
tex.write_text("\n".join(lines) + "\n")
for _ in range(2):
    proc = subprocess.run(["xelatex", "-interaction=nonstopmode", "-halt-on-error", tex.name],
                          cwd=OUT, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                          timeout=120, stdin=subprocess.DEVNULL)
    (OUT / "Table_S2a.xelatex.log").write_text(proc.stdout)
    if proc.returncode:
        raise SystemExit("Table S2a XeLaTeX failed")
