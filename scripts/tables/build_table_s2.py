#!/usr/bin/env python3
"""Build Table S2 and family_membership.tsv from frozen root-local authority.

Environment: project Python 3 (standard library only).
All inputs and outputs resolve beneath the current release root.
"""

from __future__ import annotations

import csv
import hashlib
import json
import re
import statistics
from collections import Counter, defaultdict
from pathlib import Path
import os
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from copy_table import read_copies, read_tblastn_presence  # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "release"))
from blectin_domain_counts import load_counts, table_s3_tree_clause, verify_table_s3_caption  # noqa: E402


ROOT = Path(__file__).resolve().parents[2]
with (ROOT / "config/figure_numbers.tsv").open(newline="") as handle:
    FIGURE_LABELS = {r["figure_id"]: r["label"] for r in csv.DictReader(handle, delimiter="\t")}
if "b_lectin_domain_tree" not in FIGURE_LABELS:
    raise SystemExit("B-lectin domain tree label missing from figure map")
DATA = Path(os.environ.get("PCWDE_REPOSITORY_INPUTS", ROOT / "figshare_data/repository_inputs")).resolve()
AUTH = DATA / "tables/authority/table_s2"
SYNTENY_AUTH = DATA / "tables/authority/synteny_current"
COPY_TABLE = Path(os.environ["PCWDE_COPY_TABLE"]).resolve()
TBLASTN_PRESENCE = DATA / "shared/tblastn_presence.tsv"
BUILD = Path(os.environ.get("PCWDE_BUILD_ROOT", ROOT / "build"))
STAGING = BUILD / "table_s2/candidate"
# Rendered Figure S6 tree counts, computed once from all tips by blectin_domain_counts.py.
BLECTIN_TREE_COUNTS = BUILD / "blectin_domain_tree_counts.json"
QC = BUILD / "table_s2/qc"
RECEIPTS = BUILD / "table_s2/receipts"

FAMILIES = ("GH5_2", "PL1_8", "PL4_4")
EXPECTED = {
    "GH5_2": (427, 261, 166, 318),
    "PL1_8": (364, 301, 63, 294),
    "PL4_4": (444, 365, 79, 555),
}
FORBIDDEN = ("DCAN", "Dryocosmus cannoni", "Isocolus centaureae", "beta lectin")


from parent_proteins import origin_parent


def blectin_screen_parent(header: str) -> tuple[str, str]:
    """Read the accession field of each pinned B-lectin screen header schema."""
    parts = header.split("|")
    if len(parts) != 3:
        raise SystemExit(f"Unresolved B-lectin screen header: {header}")
    group = parts[0]
    if group in {"bacterium", "plant", "fungus", "metazoa"}:
        accession = parts[1]
        if not accession or any(char.isspace() for char in accession):
            raise SystemExit(f"Unresolved B-lectin screen accession: {header}")
    elif group == "rotifer":
        match = re.fullmatch(r"([A-Z][A-Z0-9_]*\d+\.\d+)(?:_d\d+)?", parts[2])
        if not parts[1] or not match:
            raise SystemExit(f"Unresolved rotifer screen accession: {header}")
        accession = match.group(1)
    else:
        raise SystemExit(f"Unresolved B-lectin screen group: {header}")
    return group, accession


def catalytic_parent(header: str) -> tuple[str, str]:
    """Use accession in the accepted GH30 header; PI-PLC-X has one tip/protein."""
    kind, sep, detail = header.partition("__")
    if not sep:
        raise SystemExit(f"Unresolved catalytic parent: {header}")
    if kind.startswith("GH30_3_F") or kind.startswith("PI_PLC_X_F"):
        return ("focal", header)  # accepted focal sets each have one domain per protein
    if kind.startswith("PI_PLC_X_R"):
        return ("outgroup", header)  # accepted set has no repeated-domain parent
    if kind.startswith("GH30_3_R"):
        matches = re.findall(r"(?:[A-Z]{2,4}_?)?[0-9]{5,}(?:\.[0-9]+)?", detail)
        if len(matches) != 1:
            raise SystemExit(f"Unresolved GH30 parent accession: {header}")
        return ("outgroup", matches[0])
    raise SystemExit(f"Unresolved catalytic parent: {header}")


def catalytic_summary(family, retained_names):
    """Count the pinned Figure S4 domain set, using trimmed amino acids for length."""
    base = DATA / "deposit/base"
    unaligned = dict(fasta(base / f"sequences/catalytic_domains/{family}.unaligned_domains.faa"))
    trimmed = dict(fasta(base / f"alignments/catalytic_domains/{family}.trimmed.faa"))
    if set(unaligned) != set(trimmed):
        raise SystemExit(f"{family}: trimmed and unaligned Figure S4 tips differ")
    focal = [name for name in unaligned if name.startswith(f"{family}_F")]
    outgroup = [name for name in unaligned if name.startswith(f"{family}_R")]
    if len(focal) + len(outgroup) != len(unaligned):
        raise SystemExit(f"{family}: unclassified Figure S4 tip")
    groups = Counter()
    group_proteins = defaultdict(set)
    for name in outgroup:
        group = "bacterium"
        if family == "PI_PLC_X":
            groups["bacterium"] += 1
        else:
            match = re.match(r"GH30_3_R\d+__([a-z]+)_", name)
            if not match:
                raise SystemExit(f"{family}: outgroup tip without taxon group: {name}")
            group = match.group(1)
            groups[group] += 1
        group_proteins[group].add(catalytic_parent(name))
    aliases = {"Isocolus centaureae": "Eubothrus centaureae"}
    focal_species = {aliases.get(name.split("__", 2)[1].replace("_", " "),
                                 name.split("__", 2)[1].replace("_", " ")) for name in focal}
    genomes = len(focal_species & retained_names)
    return dict(family_id=family, counting_unit="domain", focal_count=len(focal),
                outgroup_count=len(outgroup),
                cynipoid_proteins=len({catalytic_parent(name) for name in focal}),
                outgroup_proteins=len({catalytic_parent(name) for name in outgroup}),
                screen_bacteria_proteins=len(group_proteins["bacterium"]),
                screen_eukaryote_proteins=len(set().union(*(v for k, v in group_proteins.items() if k != "bacterium"))),
                median_trimmed_length_aa=statistics.median(len(seq.replace("-", "")) for seq in trimmed.values()),
                genomes_with_occurrence=f"{genomes}/49", phylogeny_unit_count=f"{len(unaligned)} domains",
                screen_bacteria=groups["bacterium"],
                screen_eukaryote=sum(n for taxon, n in groups.items() if taxon != "bacterium"),
                screen_plant=groups["plant"], screen_fungus=groups["fungus"],
                screen_metazoa=groups["metazoa"], screen_rotifer=groups["rotifer"])


def rows(path: Path):
    with path.open(newline="", encoding="utf-8") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def write_tsv(path: Path, fieldnames: list[str], records: list[dict]):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        out = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        out.writeheader()
        out.writerows(records)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_id(prefix: str, value: str) -> str:
    return f"{prefix}-{hashlib.sha256(value.encode()).hexdigest()[:12].upper()}"


def fasta(path: Path):
    header = None
    chunks: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if header is not None:
        yield header, "".join(chunks)


def norm_species(value: str) -> str:
    value = value.replace("_", " ")
    value = re.sub(r"\s+", " ", value).strip()
    return value


def source_group_for(origin_class: str) -> str:
    if origin_class == "cynipoid":
        return "focal_genome"
    if origin_class == "bacterium":
        return "bacterial_outgroup"
    return "eukaryotic_outgroup"


def main():
    source_legacy = list(rows(AUTH / "source_panel_accepted.tsv"))
    ordinary = list(rows(AUTH / "ordinary_accepted_members.tsv"))
    beta_units = list(rows(AUTH / "beta_units.tsv"))
    domains = list(rows(AUTH / "b_lectin_domain_manifest.tsv"))

    db_species = {r["scientific_name"]: r for r in rows(AUTH / "species_metadata.tsv")}

    source_panel: list[dict] = []
    species_to_id: dict[str, str] = {}
    species_to_class: dict[str, str] = {}
    source_identity_to_id: dict[str, str] = {}
    for old_index, old in enumerate(source_legacy, start=1):
        name = old["species"]
        resolved = db_species.get(name)
        accession = old["accession_or_project"]
        if resolved:
            species_id = resolved["species_code"]
            display = resolved["scientific_name"]
        else:
            prefix = "SAMPLE" if old["source_type"] == "transcriptome" else "REF"
            species_id = stable_id(prefix, f"{name}|{accession}")
            display = name
        # Fixed accepted current-name and qualified-sample serialization.
        display = display.replace("Isocolus centaureae", "Eubothrus centaureae")
        display = display.replace("Ceroptres sp. (D. ignotum inquiline)", "Ceroptres sp. ex Druon ignotum")
        display = display.replace("Synergus sp. (B. kinseyi inquiline)", "Synergus sp. ex Belonocnema kinseyi")
        display = display.replace("Synergus sp. (D. ignotum inquiline)", "Synergus sp. ex Druon ignotum")
        display = display.replace("Synergus sp. (D. quercuslanigerum inquiline)", "Synergus sp. ex Druon quercuslanigerum")
        # Some accepted project accessions cover multiple distinct assemblies;
        # include the accepted source description to preserve one source row per record.
        identity = f"{name}|{accession}|{old['source_and_assembly']}|{old_index}"
        source_id = stable_id("SRC", identity)
        source_identity_to_id[identity] = source_id
        species_to_id.setdefault(name, species_id)
        if old["status"] == "retained" and old["source_type"] == "genome":
            species_to_class[name] = "focal_genome"
        elif old["source_type"] == "transcriptome":
            species_to_class[name] = "focal_transcriptome"
        repository = "NCBI" if accession.startswith(("GC", "PRJN")) else "ENA" if accession.startswith("PRJEB") else "author_generated"
        resource_class = {"genome": "genome", "transcriptome": "transcriptome_proteome", "reference_proteome": "reference_proteome"}[old["source_type"]]
        status = {"retained": "included", "excluded": "excluded", "reference outgroup": "considered", "transcriptome input": "considered"}[old["status"]]
        source_panel.append({
            "source_id": source_id, "resource_class": resource_class, "species_id": species_id,
            "species_display": display, "public_repository": repository,
            "public_accession": accession, "annotation_basis": old["source_and_assembly"],
            "assembly_busco_complete_pct": old["assembly_busco_C"],
            "proteome_busco_complete_pct": old["proteome_busco_C"],
            "analysis_panel_status": status,
            "exclusion_reason": old["exclusion_reason"] if status == "excluded" else "",
        })

    # Accepted legacy/header aliases and qualified-sample producer codes.
    # These are serialization crosswalks only; display taxonomy remains current.
    def id_for_source_name(name: str) -> str:
        return species_to_id[name]

    alias_to_accepted_name = {
        "Isocolus centaureae": "Eubothrus centaureae",
        "Eschatocerus acacia": "Eschatocerus acaciae",
        "Ceroptres (DI inquiline)": "Ceroptres sp. (D. ignotum inquiline)",
        "Ceroptres sp. (DI inquiline)": "Ceroptres sp. (D. ignotum inquiline)",
        "Ceroptres sp DI inquiline": "Ceroptres sp. (D. ignotum inquiline)",
        "Ceroptres ex Druon ignotum": "Ceroptres sp. (D. ignotum inquiline)",
        "Synergus (BK7 inquiline)": "Synergus sp. (B. kinseyi inquiline)",
        "Synergus (DI inquiline)": "Synergus sp. (D. ignotum inquiline)",
    }
    for alias, accepted_name in alias_to_accepted_name.items():
        species_to_id[alias] = id_for_source_name(accepted_name)
    producer_code_to_name = {
        "BK_synergus": "Synergus sp. (B. kinseyi inquiline)",
        "DI_ceroptres": "Ceroptres sp. (D. ignotum inquiline)",
        "DI_synergus": "Synergus sp. (D. ignotum inquiline)",
        "DL_synergus": "Synergus sp. (D. quercuslanigerum inquiline)",
    }
    for code, accepted_name in producer_code_to_name.items():
        species_to_id[code] = id_for_source_name(accepted_name)

    membership: list[dict] = []
    seen_membership_keys: set[tuple] = set()

    def add_member(species_id, family, sequence_id, parent, unit, group, occurrence, phylogeny, status="included", reason=""):
        key = (species_id, family, sequence_id, unit, phylogeny)
        if key in seen_membership_keys:
            return
        seen_membership_keys.add(key)
        membership.append({
            "membership_id": stable_id("MEM", "|".join(map(str, key))),
            "species_id": species_id, "family_id": family, "sequence_id": sequence_id,
            "parent_locus_id": parent, "unit_type": unit, "source_group": group,
            "accepted_occurrence": "TRUE" if occurrence else "FALSE", "phylogeny_set": phylogeny,
            "inclusion_status": status, "exclusion_reason": reason,
        })

    # Accepted focal protein/interface set (335 ordinary units). This is an
    # analysis membership interface, not the final-49 genome-occurrence
    # authority; occurrence is added below from the corrected synteny locus
    # detail.
    for r in ordinary:
        species_id = species_to_id.get(r["source_species_code"], r["source_species_code"])
        parent = r["stable_gene_id"] or r["stable_transcript_id"] or r["stable_protein_id"]
        unit = "locus" if r["stable_gene_id"] else "protein"
        group = "focal_transcriptome" if not r["stable_gene_id"] else "focal_genome"
        add_member(species_id, r["family"], r["stable_protein_id"], parent, unit, group, False, "accepted_protein_interface")

    # Origin-tree domain units. Headers themselves are the accepted crosswalk.
    origin_summary = []
    for family in FAMILIES:
        entries = list(fasta(DATA / f"deposit/base/sequences/origin_context_{family}_full.faa"))
        focal = 0
        for header, sequence in entries:
            parts = header.split("|")
            is_legacy_focal = parts[0] == "cynipoid"
            is_ena_focal = parts[0].endswith("_ENA")
            if is_legacy_focal:
                focal += 1
                species_name = norm_species(re.sub(r"_\d+$", "", parts[1]))
                species_id = species_to_id.get(species_name)
                if species_id is None:
                    # Current focal sources are represented deterministically by accepted display identity.
                    species_id = stable_id("REF", species_name)
                    source_panel.append({
                        "source_id": stable_id("SRC", f"origin-focal|{species_name}"),
                        "resource_class": "reference_proteome", "species_id": species_id,
                        "species_display": species_name, "public_repository": "author_generated",
                        "public_accession": parts[1], "annotation_basis": "accepted origin-tree focal source",
                        "assembly_busco_complete_pct": "", "proteome_busco_complete_pct": "",
                        "analysis_panel_status": "considered", "exclusion_reason": "",
                    })
                    species_to_id[species_name] = species_id
                group = "focal_genome"
            elif is_ena_focal:
                focal += 1
                meta_name = next((r["species"] for r in rows(AUTH / "metadata.tsv") if r["protein_id"] == header), "")
                species_name = meta_name or parts[0].removesuffix("_ENA")
                species_id = species_to_id.get(species_name) or stable_id("REF", species_name)
                group = "focal_genome"
            else:
                origin_class, accession, taxon = parts[0], parts[1], parts[2]
                species_name = norm_species(taxon)
                species_id = stable_id("REF", accession)
                group = source_group_for(origin_class)
                identity = f"origin|{accession}"
                if identity not in source_identity_to_id:
                    source_identity_to_id[identity] = stable_id("SRC", identity)
                    source_panel.append({
                        "source_id": source_identity_to_id[identity], "resource_class": "reference_proteome",
                        "species_id": species_id, "species_display": species_name,
                        "public_repository": "NCBI", "public_accession": accession,
                        "annotation_basis": "accepted origin-tree reference sequence", "assembly_busco_complete_pct": "",
                        "proteome_busco_complete_pct": "", "analysis_panel_status": "considered", "exclusion_reason": "",
                    })
            add_member(species_id, family, header, "", "domain", group, False, "origin")
        lengths = [len(seq.replace("-", "")) for _, seq in entries]
        med = statistics.median(lengths)
        total, exp_focal, exp_out, exp_med = EXPECTED[family]
        assert (len(entries), focal, len(entries) - focal, med) == (total, exp_focal, exp_out, exp_med)
        focal_parents = {origin_parent(header) for header, _ in entries
                         if header.startswith("cynipoid|") or header.split("|", 1)[0].endswith("_ENA")}
        outgroup_parents = {origin_parent(header) for header, _ in entries
                            if not (header.startswith("cynipoid|") or header.split("|", 1)[0].endswith("_ENA"))}
        origin_summary.append({"family_id": family, "counting_unit": "domain", "focal_count": focal,
                               "outgroup_count": len(entries)-focal,
                               "cynipoid_proteins": len(focal_parents), "outgroup_proteins": len(outgroup_parents),
                               "median_trimmed_length_aa": str(int(med)),
                               "genomes_with_occurrence": "", "phylogeny_unit_count": f"{len(entries)} domains"})

    # Accepted B-lectin occurrence units, including the explicit 19 tree-ineligible units.
    eligible = {r["biological_unit_id"] for r in domains if r["tip_class"] == "focal"}
    for r in beta_units:
        species_id = species_to_id.get(r["species"]) or stable_id("REF", r["species"])
        eligible_flag = r["biological_unit_id"] in eligible
        add_member(species_id, "B_lectin", r["biological_unit_id"], r["biological_unit_id"], "protein",
                   "focal_genome" if r["source_type"] == "genome" else "focal_transcriptome", False,
                   "none", "included" if eligible_flag else "excluded",
                   "" if eligible_flag else "No qualifying recoverable PF01453 domain; retained for occurrence")
    for r in domains:
        if r["tip_class"] == "focal":
            species_id = species_to_id.get(r["species"]) or stable_id("REF", r["species"])
            group = "focal_genome"
        else:
            species_id = stable_id("REF", r["representative_record_id"])
            group = "bacterial_outgroup"
            identity = f"b_lectin|{r['representative_record_id']}"
            if identity not in source_identity_to_id:
                source_identity_to_id[identity] = stable_id("SRC", identity)
                source_panel.append({
                    "source_id": source_identity_to_id[identity], "resource_class": "reference_proteome",
                    "species_id": species_id, "species_display": r["species"], "public_repository": "NCBI",
                    "public_accession": r["representative_record_id"], "annotation_basis": "accepted B-lectin outgroup sequence",
                    "assembly_busco_complete_pct": "", "proteome_busco_complete_pct": "",
                    "analysis_panel_status": "considered", "exclusion_reason": "",
                })
        add_member(species_id, "B_lectin", r["tip_id"], r["biological_unit_id"], "domain", group, False, "b_lectin_domain")

    # Species-by-family copy inventory: the retained 49-genome panel counted
    # from the shared copy table, the same source as Figures 1 and 2.
    retained = [r for r in source_legacy if r["status"] == "retained" and r["source_type"] == "genome"]
    assert len(retained) == 49 and sum(r["species"] == "Telenomus remus" for r in retained) == 1
    copies = read_copies(COPY_TABLE)
    family_map = {"B-lectin": "B_lectin"}
    copy_counts = Counter()
    retained_names = {r["species"] for r in retained}
    for r in copies:
        family = family_map.get(r["family"], r["family"])
        species_name = norm_species(r["species"])
        assert species_name in retained_names, species_name
        sid = species_to_id[species_name]
        copy_counts[sid, family] += 1
        add_member(sid, family, r["protein_id"], r["locus_id"], "gene_copy",
                   "focal_genome", True, "final_49_copies")
    copy_totals = {f: sum(n for (sid, fam), n in copy_counts.items() if fam == f) for f in (*FAMILIES, "B_lectin")}
    assert copy_totals == {"GH5_2": 187, "PL1_8": 219, "PL4_4": 224, "B_lectin": 46}, copy_totals
    inventory = []
    for r in retained:
        sid = species_to_id[r["species"]]
        inventory.append({"species_id": sid, "species_display": next(x["species_display"] for x in source_panel if x["species_id"] == sid),
                          "gh5_2_copies": copy_counts[sid, "GH5_2"], "pl1_8_copies": copy_counts[sid, "PL1_8"],
                          "pl4_4_copies": copy_counts[sid, "PL4_4"], "b_lectin_copies": copy_counts[sid, "B_lectin"]})

    # Presence (rule b): an accepted gene copy, or a validated tBLASTn locus where no copy exists.
    tblastn = read_tblastn_presence(TBLASTN_PRESENCE)
    genome_counts = {f: sum(x[f.lower()+"_copies"] > 0 or (x["species_display"], f) in tblastn for x in inventory) for f in FAMILIES}
    b_genomes = sum(x["b_lectin_copies"] > 0 or (x["species_display"], "B-lectin") in tblastn for x in inventory)
    summary = origin_summary + [{"family_id": "B_lectin", "counting_unit": "protein", "focal_count": 67,
                                 "outgroup_count": 1, "cynipoid_proteins": 67, "outgroup_proteins": 1,
                                 "median_trimmed_length_aa": "",
                                 "genomes_with_occurrence": f"{b_genomes}/49", "phylogeny_unit_count": "69 domains (67 focal; 2 outgroup)"}]
    for r in summary:
        if r["family_id"] in genome_counts:
            r["genomes_with_occurrence"] = f"{genome_counts[r['family_id']]}/49"
    domain_groups = Counter((r["family_id"], r["source_group"]) for r in membership if r["unit_type"] == "domain")
    for r in summary:
        family = r["family_id"]
        if family in FAMILIES:
            r["screen_bacteria"] = domain_groups[family, "bacterial_outgroup"]
            r["screen_eukaryote"] = domain_groups[family, "eukaryotic_outgroup"]
            r["screen_bacteria_proteins"] = len({origin_parent(header) for header, _ in fasta(DATA / f"deposit/base/sequences/origin_context_{family}_full.faa") if header.startswith("bacterium|")})
            r["screen_eukaryote_proteins"] = r["outgroup_proteins"] - r["screen_bacteria_proteins"]
            if r["screen_bacteria"] + r["screen_eukaryote"] != int(r["outgroup_count"]):
                raise SystemExit(f"{family}: screen groups do not sum to outgroup domains")
    b_screen_headers = [name for name, _ in fasta(DATA / "supplement_facts/b_lectin_outgroup_screen_domains.faa") if "|" in name]
    b_screen = Counter(name.split("|", 1)[0] for name in b_screen_headers)
    if set(b_screen) - {"cynipoid", "bacterium", "plant", "fungus", "metazoa", "rotifer"}:
        raise SystemExit(f"B-lectin screen has unexpected taxon groups: {sorted(b_screen)}")
    summary[-1]["screen_bacteria"] = b_screen["bacterium"]
    summary[-1]["screen_eukaryote"] = sum(b_screen[g] for g in ("plant", "fungus", "metazoa", "rotifer"))
    b_screen_proteins = defaultdict(set)
    for name in b_screen_headers:
        group, accession = blectin_screen_parent(name)
        b_screen_proteins[group].add(accession)
    summary[-1]["screen_bacteria_proteins"] = len(b_screen_proteins["bacterium"])
    summary[-1]["screen_eukaryote_proteins"] = sum(len(b_screen_proteins[group]) for group in ("plant", "fungus", "metazoa", "rotifer"))
    for taxon in ("plant", "fungus", "metazoa", "rotifer"):
        summary[-1][f"screen_{taxon}"] = b_screen[taxon]
    summary += [catalytic_summary(family, retained_names) for family in ("GH30_3", "PI_PLC_X")]
    if [(r["family_id"], int(r["focal_count"]), int(r["outgroup_count"])) for r in summary[-2:]] != [
            ("GH30_3", 16, 188), ("PI_PLC_X", 7, 35)]:
        raise SystemExit("Figure S4 catalytic domain counts differ from pinned sets")

    source_panel.sort(key=lambda r: (r["species_display"], r["public_accession"], r["source_id"]))
    membership.sort(key=lambda r: (r["family_id"], r["species_id"], r["unit_type"], r["sequence_id"]))
    inventory.sort(key=lambda r: r["species_display"])
    write_tsv(STAGING / "data/source_panel.tsv", list(source_panel[0]), source_panel)
    write_tsv(STAGING / "data/family_membership.tsv", list(membership[0]), membership)
    base_fields = ["family_id", "counting_unit", "focal_count", "outgroup_count", "median_trimmed_length_aa", "genomes_with_occurrence", "phylogeny_unit_count"]
    protein_fields = ["cynipoid_proteins", "outgroup_proteins", "screen_bacteria_proteins", "screen_eukaryote_proteins"]
    summary_fields = [*base_fields, "screen_bacteria", "screen_eukaryote", "screen_plant", "screen_fungus", "screen_metazoa", "screen_rotifer", *protein_fields]
    write_tsv(STAGING / "data/domain_summary.tsv", summary_fields, summary)
    table_fields = [*base_fields, "outgroup_screen_hits_by_group", *protein_fields]
    table = [{**{k: r[k] for k in base_fields},
              "outgroup_screen_hits_by_group": f"bacteria {r['screen_bacteria']}; eukaryotes {r['screen_eukaryote']}",
              **{k: r[k] for k in protein_fields}}
             for r in summary]
    # The full 49-genome copy inventory is retained as data, outside the reader table.
    assert len({r["species_display"] for r in inventory}) == len(inventory)
    write_tsv(STAGING / "data/accepted_copy_counts.tsv", list(inventory[0]), inventory)
    write_tsv(STAGING / "supplementary_tables/Table_S2.tsv", table_fields, table)
    with (ROOT / "config/table_numbers.tsv").open(newline="") as handle:
        table_labels = {r["table_id"]: r["label"] for r in csv.DictReader(handle, delimiter="\t")}
    tree_counts = load_counts(BLECTIN_TREE_COUNTS)
    caption = (f"Table {table_labels['outgroup']}. PCWDE family datasets and taxonomic occurrence. Counts are distinct source-qualified parent proteins; the median length measures trimmed domains. Domain-tree tips "
               "are distinguished from accepted protein/locus occurrence units. Outgroup screen hits report "
               "bacterial and eukaryotic proteins in the pinned screen sets; they need not equal tree outgroup proteins. "
               "The 49-genome gene-copy inventory is in the data deposit. B-lectin occurrence "
               "comprises 67 focal loci/proteins, of which 48 yielded 67 qualifying PF01453 domains; 19 remained "
               "tree-ineligible. One bacterial source contributed two domains, giving 69 domains; "
               f"{table_s3_tree_clause(tree_counts, FIGURE_LABELS['b_lectin_domain_tree'])}.\n")
    verify_table_s3_caption(caption, tree_counts, FIGURE_LABELS["b_lectin_domain_tree"])
    cap_path = STAGING / "supplementary_tables/Table_S2_caption.txt"
    cap_path.write_text(caption, encoding="utf-8")

    # Cross-file and public-content QC.
    source_ids = {r["species_id"] for r in source_panel}
    missing_fk = sorted({r["species_id"] for r in membership} - source_ids)
    assert not missing_fk, missing_fk
    assert len(beta_units) == 67
    assert sum(r["tip_class"] == "focal" for r in domains) == 67
    assert sum(r["tip_class"] != "focal" for r in domains) == 2
    outputs = [STAGING / "data/source_panel.tsv", STAGING / "data/family_membership.tsv",
               STAGING / "data/domain_summary.tsv", STAGING / "data/accepted_copy_counts.tsv",
               STAGING / "supplementary_tables/Table_S2.tsv", cap_path]
    public_text = "\n".join(p.read_text(encoding="utf-8") for p in outputs)
    assert not any(token.lower() in public_text.lower() for token in FORBIDDEN)
    QC.mkdir(parents=True, exist_ok=True)
    qc_records = [
        {"check": "origin_counts", "value": "427;364;444", "status": "PASS"},
        {"check": "origin_partition", "value": "261+166;301+63;365+79", "status": "PASS"},
        {"check": "b_lectin_accounting", "value": "67 units;48 eligible;19 ineligible;67+2 domains", "status": "PASS"},
        {"check": "final_genome_panel", "value": "49 (48 cynipoids + Telenomus remus)", "status": "PASS"},
        {"check": "final_panel_copy_totals", "value": "676=187+219+224+46", "status": "PASS"},
        {"check": "final_panel_occurrence_authority", "value": "built copy table (scripts/members/build_copy_table.py from members/)", "status": "PASS"},
        {"check": "membership_rows", "value": str(len(membership)), "status": "PASS"},
        {"check": "species_fk_missing", "value": "0", "status": "PASS"},
        {"check": "forbidden_tokens", "value": "0", "status": "PASS"},
    ]
    write_tsv(QC / "Table_S2_validation.tsv", ["check", "value", "status"], qc_records)
    RECEIPTS.mkdir(parents=True, exist_ok=True)
    receipt = {
        "artifact_id": "ART-TABLE-S2", "status": "staged_qc_pass",
        "species_metadata_sha256": sha256(AUTH / "species_metadata.tsv"),
        "inputs": ({p.name: sha256(p) for p in sorted(AUTH.iterdir()) if p.is_file()} |
                   {"built:members/copy_placements.tsv": sha256(COPY_TABLE),
                    "built:blectin_domain_tree_counts.json": sha256(BLECTIN_TREE_COUNTS),
                    "synteny_current/final_49_roster.tsv": sha256(SYNTENY_AUTH / "final_49_roster.tsv"),
                    **{f"deposit/base/{kind}/catalytic_domains/{family}.{suffix}.faa":
                       sha256(DATA / f"deposit/base/{kind}/catalytic_domains/{family}.{suffix}.faa")
                       for kind, suffix in (("sequences", "unaligned_domains"), ("alignments", "trimmed"))
                       for family in ("GH30_3", "PI_PLC_X")},
                    "supplement_facts/b_lectin_outgroup_screen_domains.faa":
                    sha256(DATA / "supplement_facts/b_lectin_outgroup_screen_domains.faa")}),
        "outputs": {str(p.relative_to(STAGING)): sha256(p) for p in outputs},
        "counts": {"source_rows": len(source_panel), "membership_rows": len(membership), "table_rows": len(table)},
    }
    (RECEIPTS / "ART-TABLE-S2.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(receipt["counts"], sort_keys=True))


if __name__ == "__main__":
    main()
