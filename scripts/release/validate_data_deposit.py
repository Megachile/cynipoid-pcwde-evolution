#!/usr/bin/env python3
import csv
import ast
import hashlib
import re
import os
import sys
from pathlib import Path
from openpyxl import load_workbook

sys.path.insert(0, str(Path(__file__).resolve().parent))
from data_root import validate_projected_inputs
from package_archive import validate_archive
from frozen_fasta_ids import declared_raw_sources, validate_raw

def s4_sequence_digest_paths():
    """Read the producer's sole path set without executing its build script."""
    producer = Path(__file__).with_name("build_data_deposit.py")
    tree = ast.parse(producer.read_text(encoding="utf-8"), filename=str(producer))
    matches = [node.value for node in tree.body if isinstance(node, ast.Assign)
               and any(isinstance(target, ast.Name) and target.id == "S4_SEQUENCE_DIGEST_TSVS"
                       for target in node.targets)]
    if len(matches) != 1:
        raise SystemExit("missing unique producer S4_SEQUENCE_DIGEST_TSVS")
    paths = ast.literal_eval(matches[0])
    if not isinstance(paths, set) or not paths or any(not isinstance(path, str) for path in paths):
        raise SystemExit("invalid producer S4_SEQUENCE_DIGEST_TSVS")
    return {"repository_inputs/" + path for path in paths}

def check_public_tsv_headers(root):
    allowed = s4_sequence_digest_paths()
    failures = []
    for path in root.rglob("*.tsv"):
        relative = path.relative_to(root).as_posix()
        if relative == "PROJECTED_INPUTS.tsv":
            # This is the file-level checksum manifest, validated against every
            # projected member; it is not a scientific row-hash table.
            validate_projected_inputs(root)
            continue
        if path.name in {"MANIFEST.tsv", "DATA_DICTIONARY.tsv"}:
            continue
        with path.open(newline="", encoding="utf-8") as handle:
            header = next(csv.reader(handle, delimiter="\t"), [])
        forbidden = [field for field in header if "sha256" in field.casefold()
                     and not (relative in allowed and field == "sequence_sha256")]
        if forbidden:
            failures.append(f"no row-level hash fields: {relative}: {','.join(forbidden)}")
    return failures

if len(sys.argv) == 3 and sys.argv[1] == "--check-public-tsv-headers":
    problems = check_public_tsv_headers(Path(sys.argv[2]))
    if problems:
        raise SystemExit("\n".join(problems))
    print("PUBLIC_TSV_HEADERS_PASS")
    raise SystemExit(0)

REPO = Path(os.environ.get("PCWDE_REPO_ROOT", Path(__file__).resolve().parents[2])).resolve()
BASE = Path(os.environ.get("PCWDE_DEPOSIT_OUT", REPO / "build/data_deposit")).resolve()
ROOT = BASE / "pcwde_supporting_data"
failures = []

def sha256_file(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()

def check(ok, name):
    if not ok:
        failures.append(name)

files = sorted(p for p in ROOT.rglob("*") if p.is_file())
check(not any(p.is_symlink() for p in ROOT.rglob("*")), "no symlinks")
check(not (ROOT / "evidence").exists(), "no internal evidence subtree")
intron_inputs = ROOT / "repository_inputs/analysis/introns"
check(not (intron_inputs / "sources").exists(), "no whole intron GFF or assembly source tree")
check((intron_inputs / "frozen/FROZEN_MANIFEST.tsv").is_file(), "frozen intron layer-one manifest")

with (ROOT / "MANIFEST.tsv").open(newline="", encoding="utf-8") as fh:
    manifest = list(csv.DictReader(fh, delimiter="\t"))
check(len(files) == len(manifest) + 2, "manifest covers every package file except itself and SHA256SUMS")
check(len(manifest) > 65, "refreshed package includes executable repository inputs")
for row in manifest:
    p = ROOT / row["relative_path"]
    check(p.is_file(), f"manifest file exists: {row['relative_path']}")
    if p.is_file():
        check(p.stat().st_size == int(row["bytes"]), f"manifest size: {row['relative_path']}")
        check(sha256_file(p) == row["sha256"], f"manifest hash: {row['relative_path']}")

# The report-only chromosome-gain analysis has the same sealed input route as
# the manuscript figures, including its genome-wide membership and AGP evidence.
with (REPO / "config/chromosome_gain_sources.tsv").open(newline="", encoding="utf-8") as fh:
    chromosome_gain_sources = list(csv.DictReader(fh, delimiter="\t"))
check(len(chromosome_gain_sources) == 59 and len({r["role"] for r in chromosome_gain_sources}) == 59,
      "chromosome-gain source manifest cardinality")
for row in chromosome_gain_sources:
    path = ROOT / "repository_inputs/chromosome_gain" / row["role"]
    check(path.is_file() and not path.is_symlink(), f"chromosome-gain input exists: {row['role']}")
    if path.is_file():
        check(path.stat().st_size == int(row["deposit_bytes"]), f"chromosome-gain input size: {row['role']}")
        check(hashlib.sha256(path.read_bytes()).hexdigest() == row["deposit_sha256"], f"chromosome-gain input hash: {row['role']}")
with (ROOT / "data/chromosome_gain/ecami_scores.tsv").open(newline="", encoding="utf-8") as fh:
    gain_scores = list(csv.DictReader(fh, delimiter="\t"))
check(sum(r["kind"] == "native" for r in gain_scores) >= 10, "eCAMI native eligibility floor")
check(sum(r["kind"] == "focal" for r in gain_scores) == 4, "four focal gain summaries")
check(sum(r["kind"] == "scan" for r in gain_scores) == 3, "three focal eCAMI sensitivities")
for name in ("ecami_tip_states.tsv", "ecami_node_changes.tsv", "background_eligibility.tsv",
             "hand_checks.tsv", "focal_summary.tsv", "copy_strata.tsv", "focal_scan_sensitivity.tsv",
             "focal_component_sizes.tsv", "focal_scan_eligibility.tsv", "focal_scan_max_diagnostic.tsv"):
    check((ROOT / "data/chromosome_gain" / name).is_file(), f"eCAMI deposit output: {name}")
check(not any(token in row["relative_path"] for row in manifest
              for token in ("orthogroup_scores", "all_member_states", "membership_table.tsv")),
      "retired chromosome-gain files absent")

for path in ROOT.rglob("*.tsv"):
    with path.open(newline="", encoding="utf-8") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
        check(bool(header), f"nonempty TSV: {path.name}")
        if header:
            width = len(header)
            check(all(len(row) == width for row in reader), f"rectangular TSV: {path.name}")

raw_sources = declared_raw_sources(REPO)
for path in ROOT.rglob("*.faa"):
    with path.open(encoding="utf-8") as handle:
        ids = [line[1:].strip().split()[0] for line in handle if line.startswith(">")]
    check(bool(ids), f"nonempty FASTA: {path.name}")
    relative = path.relative_to(ROOT).as_posix()
    if relative in raw_sources:
        try:
            digest, scan_relative, scan_digest = raw_sources[relative]
            validate_raw(path, digest, ROOT / scan_relative, scan_digest)
        except ValueError as error:
            failures.append(f"declared raw FASTA ID transform: {relative}: {error}")
    else:
        check(len(ids) == len(set(ids)), f"unique FASTA IDs: {path.name}")

for path in ROOT.rglob("*.treefile"):
    text = path.read_text(encoding="utf-8").strip()
    check(text.startswith("(") and text.endswith(";"), f"Newick framing: {path.name}")

text_suffixes = {".tsv", ".txt", ".md", ".faa", ".fa", ".fasta", ".treefile", ".nwk", ".tree", ".csv", ".r", ".py", ".bash", ".sbatch"}
global_patterns = [
    (r"(?i)claude|codex|chatgpt|large language model|needs_review|pipeline_complete", "no agent/process marks"),
    (r"/(?:" + "carc|easley|users|root" + r")/", "no private absolute paths"),
]
removed_taxon_pattern = (r"(?i)dryocosmus[ _]cannoni|\bDCAN\b", "no removed taxon")
reader_layer_patterns = [
    (r"(?i)author_generated|author-accepted|author_accepted", "no ambiguous/internal authority labels"),
    (r"(?i)\blocked\b", "no internal lock terminology"),
    (r"(?i)\bAU[- _]?(?:test|topology|selection)?\b", "no AU-test route"),
]
for path in files:
    if path.suffix.lower() not in text_suffixes or path.name in {"MANIFEST.tsv", "SHA256SUMS.txt"}:
        continue
    patterns = list(global_patterns)
    relative_text_path = path.relative_to(ROOT).as_posix()
    if not relative_text_path.startswith("repository_inputs/analysis/synteny/"):
        patterns.append(removed_taxon_pattern)
    if not relative_text_path.startswith("repository_inputs/"):
        patterns.extend(reader_layer_patterns)
    found = set()
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            for pattern, label in patterns:
                if label not in found and re.search(pattern, line):
                    found.add(label)
    for label in found:
        check(False, f"{label}: {path.relative_to(ROOT)}")

try:
    validate_projected_inputs(ROOT)
except SystemExit as error:
    failures.append(f"logical-input closure: {error}")

build = BASE.parent
for name, relative in (
    ("Table_S5.tsv", "introns/table_s5/Table_S5.tsv"),
    ("GH5_2_intron_metadata.tsv", "introns/observations/GH5_2_intron_metadata.tsv"),
    ("PL1_8_intron_metadata.tsv", "introns/observations/PL1_8_intron_metadata.tsv"),
    ("PL4_4_intron_metadata.tsv", "introns/observations/PL4_4_intron_metadata.tsv"),
    ("B_lectin_intron_metadata.tsv", "introns/observations/B_lectin_intron_metadata.tsv"),
    ("per_intron_observations.tsv", "introns/observations/per_intron_observations.tsv"),
    ("recurrent_position_summary.tsv", "introns/summary/recurrent_position_summary.tsv"),
    ("pcwde_gc3_per_copy.tsv", "introns/gc3/pcwde_gc3_per_copy.tsv"),
    ("blectin_gc3_per_copy.tsv", "introns/gc3/blectin_gc3_per_copy.tsv"),
    ("pcwde_species_pairs.tsv", "introns/gc3/pcwde_species_pairs.tsv"),
    ("blectin_gc3_comparison.tsv", "introns/gc3/blectin_gc3_comparison.tsv"),
):
    packaged = ROOT / "data/introns" / name
    source = build / relative
    check(packaged.is_file() and source.is_file() and sha256_file(packaged) == sha256_file(source), f"built intron payload: {name}")

with (ROOT / "data/source_panel.tsv").open(newline="", encoding="utf-8") as fh:
    panel = list(csv.DictReader(fh, delimiter="\t"))
may = [r for r in panel if r["species_display"] == "Mayetiola destructor" and r["source_accession_or_identifier"] == "GCA_000149185.1"]
check(len(may) == 1 and may[0]["assembly_busco_complete_pct"] == "95.7" and may[0]["proteome_busco_complete_pct"] == "92.5", "Mayetiola BUSCO 95.7/92.5")

with (ROOT / "data/selection/selection_reader_table.tsv").open(newline="", encoding="utf-8") as fh:
    selection = {r["gene_family"]: r for r in csv.DictReader(fh, delimiter="\t")}
expected = {
    "GH5_2": ("257", "317", "216", "6", "12"),
    "PL1_8": ("300", "331", "212", "8", "26"),
    "PL4_4": ("354", "551", "417", "3", "24"),
}
for fam, values in expected.items():
    row = selection.get(fam, {})
    got = tuple(row.get(k) for k in ("analyzed_sequences", "retained_codon_sites", "purifying_sites", "diversifying_sites", "significant_terminal_branches"))
    check(got == values, f"current selection values: {fam}")

with (ROOT / "data/selection/catalytic_reader_table.tsv").open(newline="", encoding="utf-8") as fh:
    catalytic = list(csv.DictReader(fh, delimiter="\t"))
with (ROOT / "data/selection/catalytic_residue_summary.tsv").open(newline="", encoding="utf-8") as fh:
    catalytic_summary = list(csv.DictReader(fh, delimiter="\t"))
current_sites = [r for r in catalytic_summary if r["alignment_set"] == "de_novo" and r["site_role"] == "catalytic"]
check(len(catalytic) == len(current_sites) == 5, "five current catalytic reader sites")
check([(r["gene_family"], r["reference_residue"], r["expected_residues"], r["observable_sequences"], r["total_sequences"])
       for r in catalytic] ==
      [(r["family"], r["reference_residue"], r["n_expected_residue"], r["n_site_observed"], r["n_total_sequences"])
       for r in current_sites], "catalytic reader matches current per-site intermediate")

expression_path = ROOT / "data/expression/gene_by_sample_tpm.tsv"
with expression_path.open(newline="", encoding="utf-8") as fh:
    expression = list(csv.DictReader(fh, delimiter="\t"))
expected_expression_fields = ["family", "species", "gene_or_protein_id", "parent_gene_or_locus_id", "sample", "study", "life_stage", "tissue", "tpm"]
check(list(expression[0]) == expected_expression_fields, "expression has the approved nine public columns")
check(len(expression) == 4362, "exact 4,362 gene-by-sample expression rows after adult BIP exclusion")
check(len({(r["species"], r["sample"]) for r in expression}) == 158, "exact 158 expression samples")
check(not {"BIP_adult_OVX", "BIP_adult_VGX"} & {r["sample"] for r in expression}, "adult BIP samples excluded")
check({r["study"] for r in expression} == {"This study", "PRJNA1021581", "PRJNA493150", "PRJNA490752"}, "expression uses public study labels")
check(all(re.fullmatch(r"\d+\.\d", r["tpm"]) for r in expression), "expression TPM rounded to one decimal")
check(not any(r["life_stage"] == "larva" and r["tissue"] == "gall" for r in expression), "larval gall label removed")
check(any(r["tissue"] == "whole_body (in gall)" for r in expression), "whole-body-in-gall label present")

synteny_path = ROOT / "data/synteny/conserved_synteny_blocks.xlsx"
wb = load_workbook(synteny_path, read_only=True, data_only=True)
check(len(wb.sheetnames) == 78 and wb.sheetnames[0] == "Index", "synteny workbook has index plus 77 detailed blocks")
check(sum(1 for row in wb["Index"].iter_rows(min_row=2, values_only=True) if row[7]) == 77, "77 cross-genome blocks indexed")
check(not any((ROOT / "data/synteny" / name).exists() for name in ("block_display_crosswalk.tsv", "flanking_gene_maps.tsv", "locus_assignments.tsv")), "internal synteny projections removed")

failures.extend(check_public_tsv_headers(ROOT))

try:
    validate_archive(BASE)
except SystemExit as error:
    failures.append(f"archive: {error}")

metadata = ROOT / "FIGSHARE_METADATA.md"
metadata_text = metadata.read_text(encoding="utf-8") if metadata.is_file() else ""
check("10.6084/m9.figshare.33265659" in metadata_text, "DOI metadata")
check("Version: 1.0.0" in metadata_text and "pcwde_supporting_data.tar.gz" in metadata_text, "version/archive metadata")

# Full transfer-context evidence must accompany the representative display subsets.
from Bio import Phylo, SeqIO
mapping_path = ROOT / "sidecars/origin_context_full_tip_mapping.tsv"
check(mapping_path.is_file(), "complete origin tip mapping present")
origin_mapping = list(csv.DictReader(mapping_path.open(), delimiter="\t")) if mapping_path.is_file() else []
for family, expected_count in (("GH5_2", 427), ("PL1_8", 364), ("PL4_4", 444)):
    tree_path = ROOT / f"trees/origin_context_{family}_full.treefile"
    alignment_path = ROOT / f"alignments/origin_context_{family}_full.aln.faa"
    original_path = ROOT / f"sequences/origin_context_{family}_full.faa"
    method_path = ROOT / f"sidecars/origin_context_{family}_full_method.txt"
    check(all(p.is_file() for p in (tree_path, alignment_path, original_path, method_path)), f"complete origin payload: {family}")
    if not all(p.is_file() for p in (tree_path, alignment_path, original_path, method_path)):
        continue
    aligned = {record.id: str(record.seq) for record in SeqIO.parse(alignment_path, "fasta")}
    original = {record.id: str(record.seq) for record in SeqIO.parse(original_path, "fasta")}
    tips = [node.name for node in Phylo.read(tree_path, "newick").get_terminals()]
    rows = [row for row in origin_mapping if row["family"] == family]
    check(len(aligned) == len(original) == len(tips) == len(set(tips)) == len(rows) == expected_count, f"full origin counts: {family}")
    check(set(aligned) == set(original) == {row["alignment_sequence_id"] for row in rows}, f"original/alignment mapping: {family}")
    check(set(tips) == {row["tree_tip_id"] for row in rows}, f"tree tip mapping: {family}")
    check(all(re.sub(r"[()\[\]=]", "_", row["alignment_sequence_id"]) == row["tree_tip_id"] for row in rows), f"IQ-TREE spelling: {family}")
    normalize = lambda value: value.upper().replace("-", "").replace("*", "")
    check(all(key in original and normalize(value) == normalize(original[key]) for key, value in aligned.items()), f"unaligned/aligned sequence identity: {family}")
    method_text = method_path.read_text()
    check("Model: LG+G4" in method_text and "Ultrafast bootstrap replicates: 1000" in method_text and "Seed:" in method_text, f"origin inference method: {family}")
    for suffix in (".faa", "_domtbl.txt"):
        check((ROOT / "repository_inputs/analysis/domain_preparation" / (family + suffix)).is_file(), f"preprocessing input: {family}{suffix}")

if failures:
    print("FAIL")
    for item in failures:
        print(item)
    raise SystemExit(1)
print(f"PASS: {len(files)} package files; {len(manifest)} checksummed payload files; reader-data, archive, and scientific-value gates passed")
