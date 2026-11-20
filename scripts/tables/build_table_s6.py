#!/usr/bin/env python3
"""Build supplementary intron and expression tables from the supplied tables.

Requires Python 3 and XeLaTeX for PDF output.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import re
import statistics
import subprocess
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
with (ROOT / "config/table_numbers.tsv").open(newline="") as handle:
    TABLE_LABELS = {r["table_id"]:r["label"] for r in csv.DictReader(handle,delimiter="\t")}
DATA = Path(os.environ.get("PCWDE_REPOSITORY_INPUTS", ROOT / "figshare_data/repository_inputs")).resolve()
WORK = Path(os.environ.get("PCWDE_BUILD_ROOT", ROOT / "build")) / "tables_s5_s6"
AUTH = DATA / "tables/authority/tables_s5_s6_current"
BASELINE = AUTH / "baseline"
RETAINED_CROSSWALK_FIELDS = ["full_sample_id", "original_sample_id", "sample_id", "species_code", "species_id", "species_display"]
INTRON = AUTH / "intron_inputs"
TARGET_EXPRESSION = AUTH / "target_transcript_expression.tsv"
AUTH_S2 = DATA / "tables/authority/table_s2"

OUT = WORK / "candidate"
TABLES = OUT / "supplementary_tables"
DATA_I = OUT / "data" / "introns"
DATA_E = OUT / "data" / "expression"
EVIDENCE = WORK / "evidence"
INGRESS = WORK / "ingress"

FAMILIES = ("GH5_2", "PL1_8", "PL4_4", "B_lectin")
SOURCE_FAMILY = {"B_lectin": "CBM13"}
FAMILY_DISPLAY = {"B_lectin": "B-lectin"}
VISIBLE_TABLE_NUMBERS = {"S5": TABLE_LABELS["intron"], "S6": TABLE_LABELS["expression"]}
TARGET_SEQUENCE = "DI_ceroptres_TRINITY_DN781_c0_g1_i11.p1"
LOCKED_TARGET_TPM = {
    "RX02": "1.9", "RX03": "45.0", "RX04": "1.2",
    "WX02": "0.7", "WX04": "29.4", "WX05": "66.0",
}
CAPTION_S6 = (
    f"Table {TABLE_LABELS['expression']}b. Annotation- and catalog-based family TPM summaries for the 158-profile analysis set. Values are sums across eligible family members within a sample and are grouped by species, life stage, tissue, and family; IQR is Q1-Q3. "
    "Annotation sources vary by taxon; Belonocnema kinseyi values use the NCBI genome annotation. "
    "TPM summaries are reported to the nearest 0.1 TPM; whole_body (in gall) denotes whole-larva material collected from galls. "
    "Raw-only conflicts are retained as qualitative QC context and do not replace the reported numeric values. Full sample-family and target-sample values, stable public sample IDs, metadata, and available mapping QC are provided in the accompanying data files."
)
SAMPLE_COUNT_COLUMNS = [
    ("adult", "accessory_gland", "adult_accessory_gland"),
    ("adult", "venom_gland", "adult_venom_gland"),
    ("adult", "carcass", "adult_carcass"),
    ("adult", "whole_body", "adult_whole_body"),
    ("larva", "gut", "larval_gut"),
    ("larva", "salivary_gland", "larval_salivary_gland"),
    ("larva", "malpighian_tubule", "larval_malpighian_tubule"),
    ("larva", "whole_body", "larval_whole_body"),
]
SAMPLE_COUNT_FIELDS = ["species", "source"] + [c for _, _, c in SAMPLE_COUNT_COLUMNS]


def sample_count_matrix(sample_family, metadata, crosswalk, table_rows, source_panel):
    """Count retained profile IDs once per species and context, independent of gene rows."""
    display = {r["sample_id"]: r["species_display"] for r in crosswalk}
    if len(display) != len(crosswalk) or set(display) != {r["sample_id"] for r in metadata}:
        raise AssertionError("retained sample/display crosswalk is incomplete")
    contexts = {(stage, tissue): column for stage, tissue, column in SAMPLE_COUNT_COLUMNS}
    by_cell = defaultdict(set)
    for row in metadata:
        context = (row["life_stage"], "whole_body" if row["life_stage"] == "larva" and row["tissue"] == "gall" else row["tissue"])
        if context not in contexts:
            raise AssertionError(f"unmapped expression context: {context}")
        by_cell[(display[row["sample_id"]], context)].add(row["sample_id"])
    if sum(map(len, by_cell.values())) != len(metadata):
        raise AssertionError("a profile was assigned to multiple matrix cells")
    panel_alias = {
        "Ceroptres sp. (D. ignotum inquiline)": "Ceroptres sp. ex Druon ignotum",
        "Synergus sp. (B. kinseyi inquiline)": "Synergus sp. ex Belonocnema kinseyi",
        "Synergus sp. (D. ignotum inquiline)": "Synergus sp. ex Druon ignotum",
        "Synergus sp. (D. quercuslanigerum inquiline)": "Synergus sp. ex Druon quercuslanigerum",
    }
    panel_names = {panel_alias.get(r["species"], TAXON_FIX.get(r["species"], r["species"])): r["accession_or_project"]
                   for r in source_panel if r["source_type"] == "transcriptome"}
    source = {}
    for species in set(display.values()):
        accession = panel_names.get(species, "")
        if not accession:
            if species != "Belonocnema kinseyi":
                raise AssertionError(f"unresolved expression source: {species}")
            source[species] = "This study"
        elif species == "Dryocosmus quercuspalustris" or accession.startswith("PRJNA150"):
            source[species] = "This study"
        elif "PRJNA" in accession or "GSE" in accession:
            source[species] = re.search(r"(?:PRJNA|GSE)\d+", accession).group()
        else:
            raise AssertionError(f"unresolved expression source: {species}: {accession}")
    result = []
    for species in sorted(source):
        row = {"species": species, "source": source[species]}
        row.update({column: str(len(by_cell.get((species, (stage, tissue)), ()))) or ""
                    for stage, tissue, column in SAMPLE_COUNT_COLUMNS})
        for column in (c for _, _, c in SAMPLE_COUNT_COLUMNS):
            if row[column] == "0":
                row[column] = ""
        result.append(row)
    summary_keys = set()
    summary_counts = defaultdict(int)
    for item in table_rows:
        context = (item["life_stage"],
                   "whole_body" if item["tissue"] == "whole_body (in gall)" else item["tissue"])
        key = (item["species"], *context)
        summary_keys.add(key)
        original_key = (item["species"], item["life_stage"], item["tissue"])
        if original_key in summary_counts and summary_counts[original_key] != int(item["sample_count"]):
            raise AssertionError(f"Table S2b family count mismatch: {original_key}")
        summary_counts[original_key] = int(item["sample_count"])
    merged_counts = defaultdict(int)
    for (species, stage, tissue), count in summary_counts.items():
        merged_counts[(species, stage, "whole_body" if tissue == "whole_body (in gall)" else tissue)] += count
    for key, count in merged_counts.items():
        if len(by_cell.get((key[0], (key[1], key[2])), ())) != count:
            raise AssertionError(f"Table S2b sample_count mismatch: {key}")
    if summary_keys != {(species, *context) for species, context in by_cell}:
        raise AssertionError("Table S2b has a context absent from the sample matrix")
    if len(result) != 15 or sum(map(len, by_cell.values())) != 158:
        raise AssertionError("Table S2a matrix dimensions changed")
    return result
TAXON_FIX = {
    "Isocolus centaureae": "Eubothrus centaureae",
    "Eschatocerus acacia": "Eschatocerus acaciae",
    "Synergus ex BK": "Synergus sp. ex Belonocnema kinseyi",
    "Synergus ex DI": "Synergus sp. ex Druon ignotum",
    "Synergus ex DL": "Synergus sp. ex Druon quercuslanigerum",
    "Ceroptres ex DI": "Ceroptres sp. ex Druon ignotum",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fields: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def norm_id(value: str) -> str:
    return re.sub(r"\.1$", "", value.strip())


def display_family(value: str) -> str:
    return FAMILY_DISPLAY.get(value, value)


def clean_taxon(value: str) -> str:
    return TAXON_FIX.get(value.strip(), value.strip())


def taxon_id(value: str) -> str:
    slug = re.sub(r"[^a-z0-9]+", "-", value.lower()).strip("-")
    return f"TAXON-{slug}"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def percentile(values: list[float], p: float) -> float:
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    position = (len(ordered) - 1) * p
    lo, hi = math.floor(position), math.ceil(position)
    if lo == hi:
        return ordered[lo]
    return ordered[lo] + (ordered[hi] - ordered[lo]) * (position - lo)


def fmt_num(value: float) -> str:
    """Preserve the accepted compact formatting used by Table S6."""
    return f"{value:.3f}".rstrip("0").rstrip(".")


def fmt_tpm(value: float) -> str:
    """Apply the author-approved one-decimal display format to Table S7."""
    return f"{value:.1f}"


def tex_escape(value: object) -> str:
    text = str(value)
    replacements = {
        "\\": r"\textbackslash{}", "&": r"\&", "%": r"\%", "$": r"\$",
        "#": r"\#", "_": r"\_", "{": r"\{", "}": r"\}", "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
    }
    return "".join(replacements.get(char, char) for char in text)


def build_intron() -> tuple[list[dict], list[dict], list[dict], dict]:
    metadata = read_tsv(AUTH_S2 / "metadata.tsv")
    accepted = defaultdict(dict)
    for row in metadata:
        if row["family"] in FAMILIES[:3] and row["source"] == "genome":
            accepted[row["family"]][norm_id(row["protein_id"])] = row

    b_members = read_tsv(AUTH_S2 / "b_lectin_accepted_members.tsv")
    b_units = read_tsv(AUTH_S2 / "beta_units.tsv")
    b_alias = {}
    for row in b_members:
        b_alias[norm_id(row["representative_record_id"])] = row
        b_alias[norm_id(row["biological_unit_id"])] = row
    for row in b_units:
        member = b_alias.get(norm_id(row["biological_unit_id"]))
        if member:
            for key in ("source_record_id", "full_sequence_id", "biological_unit_id"):
                if row.get(key):
                    b_alias[norm_id(row[key])] = member

    public_map: dict[tuple[str, str], dict] = {}
    eligible_rows: list[dict] = []
    raw_by_family: dict[str, list[dict]] = {}
    producer_counts = {}
    excluded_counts = {}
    for family in FAMILIES:
        source_family = SOURCE_FAMILY.get(family, family)
        raw = read_tsv(INTRON / "extracted" / f"{source_family}_intron_metadata.tsv")
        raw_by_family[family] = raw
        raw_ids = sorted({norm_id(r["protein_id"]) for r in raw})
        for pid in raw_ids:
            matching = [r for r in raw if norm_id(r["protein_id"]) == pid]
            exemplar = matching[0]
            if family == "B_lectin":
                member = b_alias.get(pid)
                if not member:
                    continue
                public_id = member["biological_unit_id"]
                species = clean_taxon(member["species"])
            else:
                member = accepted[family].get(pid)
                if not member:
                    continue
                public_id = member["protein_id"]
                species = clean_taxon(member["species"])
            public_map[(family, pid)] = {"locus_id": public_id, "species": species}
            eligible_rows.append({
                "family_id": family,
                "locus_id": public_id,
                "species": species,
                "producer_protein_id": exemplar["protein_id"],
                "intron_count": max(int(float(r["intron_count"])) for r in matching),
                "eligibility_status": "eligible_current_membership_and_cds_recovery",
            })
        producer_counts[family] = len(raw_ids)
        excluded_counts[family] = len(raw_ids) - sum(1 for f, _ in public_map if f == family)

    # The latest mapped-position files retain positive rows only. Missing maps are legitimate NA.
    mapped = {}
    locus_classes = {}
    for family in FAMILIES:
        sf = SOURCE_FAMILY.get(family, family)
        for row in read_tsv(INTRON / f"{sf}_intron_positions.tsv"):
            key = (family, norm_id(row["protein_id"]), int(float(row["intron_number"])))
            mapped[key] = row.get("ref_aa_pos", "")
            locus_classes[key] = row.get("locus_class", "")

    observations = []
    for family in FAMILIES:
        for row in raw_by_family[family]:
            if int(float(row["intron_count"])) <= 0:
                continue
            pid = norm_id(row["protein_id"])
            info = public_map.get((family, pid))
            if not info:
                continue
            number = int(float(row["intron_number"]))
            key = (family, pid, number)
            ref = mapped.get(key, "")
            ref = "" if ref in ("", None, "-1", "-1.0") else str(int(float(ref)))
            observations.append({
                "intron_id": "", "locus_id": info["locus_id"], "family_id": family,
                "raw_cds_position_nt": "", "reference_position_aa": ref,
                "phase": int(float(row["phase_before"])),
                "species": info["species"], "locus_distribution": locus_classes.get(key, "") or "unresolved",
                "source_aa_position_approximate": int(float(row["aa_position"])),
                "donor_dinucleotide": row["donor_dinuc"], "acceptor_dinucleotide": row["acceptor_dinuc"],
            })
    observations.sort(key=lambda r: (FAMILIES.index(r["family_id"]), r["locus_id"], int(r["intron_id"].split("-")[-1]) if r["intron_id"] else 0, r["source_aa_position_approximate"]))
    for index, row in enumerate(observations, 1):
        row["intron_id"] = f"INTRON-{index:05d}"

    eligible_count = Counter(r["family_id"] for r in eligible_rows)
    positive_loci = defaultdict(set)
    total_introns = Counter()
    for row in observations:
        positive_loci[row["family_id"]].add(row["locus_id"])
        total_introns[row["family_id"]] += 1
    part_a = []
    for family in FAMILIES:
        denom = eligible_count[family]
        pos = len(positive_loci[family])
        part_a.append({
            "family_id": family, "eligible_proteins": denom, "proteins_with_introns": pos,
            "total_introns": total_introns[family],
            "protein_prevalence_pct": fmt_num(100 * pos / denom) if denom else "",
        })

    # Re-cluster mapped observations after current-membership filtering. Mode-seeking
    # follows the accepted producer implementation; recurrent means >=5 loci.
    recurrent = []
    for family in FAMILIES:
        items = [r for r in observations if r["family_id"] == family and r["reference_position_aa"]]
        alive = list(range(len(items)))
        clusters = []
        while alive:
            best = max(alive, key=lambda seed: sum(abs(int(items[j]["reference_position_aa"]) - int(items[seed]["reference_position_aa"])) <= 3 for j in alive))
            members = [j for j in alive if abs(int(items[j]["reference_position_aa"]) - int(items[best]["reference_position_aa"])) <= 3]
            center_positions = sorted(int(items[j]["reference_position_aa"]) for j in members)
            center = center_positions[len(center_positions) // 2]
            clusters.append((center, members))
            member_set = set(members)
            alive = [j for j in alive if j not in member_set]
        for center, members in sorted(clusters):
            chosen = {}
            for j in members:
                chosen.setdefault(items[j]["locus_id"], items[j])
            rows = list(chosen.values())
            if len(rows) < 5:
                continue
            phases = Counter(str(r["phase"]) for r in rows if r["phase"] in (0, 1, 2))
            loci = Counter(r["locus_distribution"] for r in rows)
            recurrent.append({
                "position_id": "", "family_id": family, "reference_position_aa": center,
                "intron_count": len(rows), "eligible_protein_count": eligible_count[family],
                "phase_summary": ";".join(f"{p}:{phases[p]}" for p in ("0", "1", "2")),
                "taxonomic_breadth": f"{len({r['species'] for r in rows})} species",
                "locus_distribution": "; ".join(f"{k}:{v}" for k, v in sorted(loci.items())),
            })
    recurrent.sort(key=lambda r: (FAMILIES.index(r["family_id"]), int(r["reference_position_aa"])))
    for index, row in enumerate(recurrent, 1):
        row["position_id"] = f"INTPOS-{index:03d}"

    table_rows = []
    for row in part_a:
        table_rows.append({"part": "A", **row})
    for row in recurrent:
        modes = Counter()
        for token in row["phase_summary"].split(";"):
            phase, count = token.split(":")
            modes[phase] = int(count)
        max_count = max(modes.values()) if modes else 0
        top = [p for p, n in modes.items() if n == max_count]
        table_rows.append({
            "part": "B", "family_id": row["family_id"], "reference_position_aa": row["reference_position_aa"],
            "observation_count": row["intron_count"],
            "observation_proportion_pct": fmt_num(100 * row["intron_count"] / row["eligible_protein_count"]),
            "phase": top[0] if len(top) == 1 else "mixed", "taxonomic_breadth": row["taxonomic_breadth"],
            "locus_distribution": row["locus_distribution"],
        })

    qc = {
        "producer_unique_proteins": producer_counts,
        "excluded_not_in_current_membership": excluded_counts,
        "eligible_proteins": dict(eligible_count),
        "proteins_with_introns": {f: len(positive_loci[f]) for f in FAMILIES},
        "positive_intron_observations": dict(total_introns),
        "mapped_observations": {f: sum(r["family_id"] == f and bool(r["reference_position_aa"]) for r in observations) for f in FAMILIES},
        "recurrent_positions": {f: sum(r["family_id"] == f for r in recurrent) for f in FAMILIES},
        "b_lectin_alias_closure": sum(1 for r in eligible_rows if r["family_id"] == "B_lectin"),
    }
    return observations, recurrent, table_rows, {"eligible_rows": eligible_rows, "qc": qc, "part_a": part_a}


def build_expression() -> tuple[list[dict], list[dict], list[dict], list[dict], list[dict], dict]:
    roster = read_tsv(BASELINE / "source_sample_roster.tsv")
    sample_ids = {r["full_sample_id"]: r["sample_id"] for r in roster}
    if len(sample_ids) != 158 or len(set(sample_ids.values())) != 158 or any(r["species_code"] == "DCAN" for r in roster):
        raise AssertionError("public 158-profile, zero-DCAN roster failed")

    metadata = read_tsv(BASELINE / "sample_metadata_and_mapping_qc.tsv")
    sample_family = read_tsv(BASELINE / "sample_family_tpm.tsv")
    gene_sample = read_tsv(BASELINE / "gene_sample_tpm.tsv")
    display_rows = read_tsv(BASELINE / "expression_sample_public_id_crosswalk.tsv")
    display_by_sample = {r["sample_id"]: r["species_display"] for r in display_rows}
    if len(metadata) != 158 or set(display_by_sample) != {r["sample_id"] for r in metadata} or set(display_by_sample) != set(sample_ids.values()):
        raise AssertionError("accepted public sample metadata/crosswalk closure failed")
    if len(sample_family) != 632 or len(gene_sample) != 4356:
        raise AssertionError("public expression baseline row contract failed")
    exclusions = read_tsv(ROOT / "config/excluded_expression_samples.tsv")
    excluded_full = {r["full_sample_id"] for r in exclusions}
    excluded_ids = set()
    if len(exclusions) != 2 or len(excluded_full) != 2 or not excluded_full.isdisjoint(sample_ids):
        raise AssertionError("excluded expression samples re-entered public roster")
    crosswalk = [{"full_sample_id": r["full_sample_id"],
                  "original_sample_id": r["original_sample_id"],
                  "sample_id": sample_ids[r["full_sample_id"]],
                  "species_code": r["species_code"],
                  "species_id": next(m["species_id"] for m in metadata if m["sample_id"] == sample_ids[r["full_sample_id"]]),
                  "species_display": display_by_sample[sample_ids[r["full_sample_id"]]]}
                 for r in roster if r["full_sample_id"] not in excluded_full]
    if len(crosswalk) != 158 or len({r["sample_id"] for r in crosswalk}) != 158:
        raise AssertionError("retained expression crosswalk is not one-to-one")
    roster = [r for r in roster if r["full_sample_id"] not in excluded_full]
    metadata = [r for r in metadata if r["sample_id"] not in excluded_ids]
    sample_family = [r for r in sample_family if r["sample_id"] not in excluded_ids]
    gene_sample = [r for r in gene_sample if r["sample_id"] not in excluded_ids]
    display_by_sample = {k: v for k, v in display_by_sample.items() if k not in excluded_ids}

    target_rows = [
        r for r in read_tsv(TARGET_EXPRESSION)
        if r["protein_id"] == TARGET_SEQUENCE
        and r["species_code"] == "DI_ceroptres"
        and r["full_sample_id"].startswith("DI_ceroptres__ceroptres:")
    ]
    target_by_original = {r["original_sample_id"]: r for r in target_rows}
    if len(target_rows) != 6 or len(target_by_original) != 6:
        raise AssertionError("target transcript authority must resolve to six unique Ceroptres samples")
    if {sample: target_by_original[sample]["tpm"] for sample in target_by_original} != LOCKED_TARGET_TPM:
        raise AssertionError("target transcript TPM receipt does not match locked values")

    target_public = {}
    for row in target_rows:
        sid = sample_ids.get(row["full_sample_id"])
        if not sid:
            raise AssertionError(f"target sample absent from locked roster: {row['full_sample_id']}")
        target_public[sid] = row
    if len(target_public) != 6:
        raise AssertionError("target transcript did not map one-to-one onto six public samples")

    family_by_key = {(r["sample_id"], r["family_id"]): r for r in sample_family}
    if len(family_by_key) != len(sample_family):
        raise AssertionError("duplicate accepted sample-family keys")
    for sid, target_row in target_public.items():
        family_row = family_by_key[(sid, "B_lectin")]
        family_row["family_tpm"] = fmt_num(float(family_row["family_tpm"]) + float(target_row["tpm"]))
        family_row["quantification_basis"] = "helixer_plus_accepted_trinity_target"

    next_index = max(int(r["gene_sample_id"].split("-")[-1]) for r in gene_sample) + 1
    transcript_id = TARGET_SEQUENCE.removesuffix(".p1")
    for sid, target_row in sorted(target_public.items()):
        gene_sample.append({
            "gene_sample_id": f"GENEEXP-{next_index:05d}",
            "sample_id": sid,
            "stable_target_sequence_id": TARGET_SEQUENCE,
            "parent_gene_or_locus_id": "",
            "transcript_sequence_id": transcript_id,
            "tpm": target_row["tpm"],
        })
        next_index += 1

    target_keys = [(r["sample_id"], r["stable_target_sequence_id"]) for r in gene_sample]
    if len(target_keys) != len(set(target_keys)):
        raise AssertionError("duplicate (sample_id, stable_target_sequence_id) expression keys")

    by_group = defaultdict(list)
    meta_by_sample = {r["sample_id"]: r for r in metadata}
    for row in sample_family:
        meta = meta_by_sample[row["sample_id"]]
        tissue = "whole_body (in gall)" if meta["life_stage"] == "larva" and meta["tissue"] == "gall" else meta["tissue"]
        key = (display_by_sample[row["sample_id"]], meta["life_stage"], tissue, row["family_id"])
        by_group[key].append(float(row["family_tpm"]))
    table_rows = []
    for key, values in sorted(by_group.items(), key=lambda kv: (kv[0][0], kv[0][1], kv[0][2], FAMILIES.index(kv[0][3]))):
        species, stage, tissue, family = key
        table_rows.append({
            "species": species, "life_stage": stage, "tissue": tissue, "family_id": family,
            "sample_count": len(values), "median_tpm": fmt_tpm(statistics.median(values)),
            "iqr_tpm": f"{fmt_tpm(percentile(values, .25))}-{fmt_tpm(percentile(values, .75))}",
            "range_tpm": f"{fmt_tpm(min(values))}-{fmt_tpm(max(values))}",
        })

    qc = {
        "profiles": len(metadata), "sample_family_rows": len(sample_family), "gene_sample_rows": len(gene_sample),
        "families_per_profile": dict(Counter(r["sample_id"] for r in sample_family)),
        "mapped_qc_rows": sum(bool(r["mapping_rate_pct"]) for r in metadata),
        "flagged_mapping_rows": sum(bool(r["qc_note"]) for r in metadata),
        "unique_expression_targets": len({r["stable_target_sequence_id"] for r in gene_sample}),
        "target_transcript_rows_added": len(target_public),
        "target_public_samples": sorted(target_public),
        "duplicate_sample_target_keys": len(target_keys) - len(set(target_keys)),
        "table_groups": len(table_rows),
    }
    return sample_family, gene_sample, metadata, table_rows, crosswalk, qc


def render_pdf(table_no: str, title: str, caption: str, sections: list[tuple[str, list[str], list[str], list[dict]]]) -> Path:
    tex_path = TABLES / f"Table_{table_no}.tex"
    pdf_path = TABLES / f"Table_{table_no}.pdf"
    visible_no = VISIBLE_TABLE_NUMBERS[table_no]
    visible_title = title.replace(f"Table {table_no}.", f"Table {visible_no}.", 1)
    visible_caption = caption.replace(f"Table {table_no}.", f"Table {visible_no}.", 1)
    lines = [
        r"\documentclass[10pt]{article}", r"\usepackage[margin=0.45in,landscape]{geometry}",
        r"\usepackage{longtable,booktabs,array,xcolor}",
        r"\setlength{\tabcolsep}{3pt}", r"\renewcommand{\arraystretch}{1.12}",
        r"\pagestyle{plain}", r"\begin{document}",
        rf"\section*{{{tex_escape(visible_title)}}}", rf"\small {tex_escape(visible_caption)}\par\medskip",
    ]
    for heading, headers, widths, rows in sections:
        colspec = "@{}" + "".join(rf">{{\raggedright\arraybackslash}}p{{{w}}}" for w in widths) + "@{}"
        lines.extend([rf"\subsection*{{{tex_escape(heading)}}}", rf"\scriptsize\begin{{longtable}}{{{colspec}}}", r"\toprule"])
        header = " & ".join(rf"\textbf{{{tex_escape(h)}}}" for h in headers) + r" \\ \midrule"
        lines.extend([header, r"\endfirsthead", r"\toprule", header, r"\endhead", r"\bottomrule", r"\endfoot"])
        for row in rows:
            lines.append(" & ".join(tex_escape(row.get(h, "")) for h in headers) + r" \\")
        lines.append(r"\end{longtable}")
    lines.append(r"\end{document}")
    tex_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    tex_env = os.environ.copy()
    tex_env["SOURCE_DATE_EPOCH"] = "1786665600"
    tex_env["FORCE_SOURCE_DATE"] = "1"
    for _ in range(2):
        subprocess.run(
            ["xelatex", "-interaction=nonstopmode", "-halt-on-error", tex_path.name],
            cwd=TABLES, check=True, stdout=subprocess.DEVNULL, env=tex_env,
        )
    if not pdf_path.exists() or pdf_path.stat().st_size == 0:
        raise AssertionError(f"missing PDF {pdf_path}")
    return pdf_path


def main() -> None:
    for directory in (TABLES, DATA_I, DATA_E, EVIDENCE, INGRESS):
        directory.mkdir(parents=True, exist_ok=True)

    observations, recurrent, table_s5, intron_aux = build_intron()
    current_intron = DATA / "tables/introns_current"
    current_table_s5 = DATA / "tables/accepted/Table_S5.tsv"
    table_s5 = read_tsv(current_table_s5)
    intron_aux["part_a"] = [row for row in table_s5 if row["part"] == "A"]
    current_metadata_proteins = set()
    s5_fields = ["part", "family_id", "eligible_proteins", "proteins_with_introns", "total_introns", "protein_prevalence_pct",
                 "reference_position_aa", "observation_count", "observation_proportion_pct", "phase", "taxonomic_breadth", "locus_distribution"]
    for family in FAMILIES[:3]:
        source = current_intron / f"{family}_intron_metadata.tsv"
        destination = DATA_I / source.name
        destination.write_bytes(source.read_bytes())
        current_metadata_proteins.update((family, row["protein_id"]) for row in read_tsv(source))
    (DATA_I / "recurrent_position_summary.tsv").write_bytes(
        (current_intron / "recurrent_position_summary.tsv").read_bytes())
    write_tsv(EVIDENCE / "intron_eligibility.tsv",
              ["family_id", "locus_id", "species", "producer_protein_id", "intron_count", "eligibility_status"], intron_aux["eligible_rows"])
    write_tsv(TABLES / "Table_S5.tsv", s5_fields, table_s5)

    sample_family, gene_sample, sample_meta, table_s6, crosswalk, expression_qc = build_expression()
    write_tsv(DATA_E / "sample_family_tpm.tsv", ["sample_family_id", "sample_id", "family_id", "family_tpm", "quantification_basis"], sample_family)
    write_tsv(DATA_E / "gene_sample_tpm.tsv",
              ["gene_sample_id", "sample_id", "stable_target_sequence_id", "parent_gene_or_locus_id", "transcript_sequence_id", "tpm"],
              gene_sample)
    write_tsv(DATA_E / "sample_metadata_and_mapping_qc.tsv",
              ["sample_id", "species_id", "life_stage", "tissue", "mapping_rate_pct", "included_in_160_profile_scope", "qc_note"], sample_meta)
    retained_sample_ids = {r["sample_id"] for r in sample_meta}
    write_tsv(EVIDENCE / "expression_sample_public_id_crosswalk.tsv", ["sample_id", "species_id", "species_display"],
              [{k: r[k] for k in ("sample_id", "species_id", "species_display")} for r in crosswalk])
    write_tsv(DATA_E / "retained_sample_crosswalk.tsv", RETAINED_CROSSWALK_FIELDS, crosswalk)
    write_tsv(TABLES / "Table_S6.tsv", ["species", "life_stage", "tissue", "family_id", "sample_count", "median_tpm", "iqr_tpm", "range_tpm"], table_s6)

    caption_s5 = (
        f"Table {TABLE_LABELS['intron']}. Intron prevalence and recurrent positions in accepted current family members with recoverable genomic CDS models. "
        "Part A distinguishes eligible proteins, proteins with at least one intron, and individual intron observations. Part B reports positions represented in at least five proteins after mapping full-length proteins to the family reference and clustering coordinates within three amino acids. "
        "Reference positions are aligned amino-acid coordinates; phase is the coding-frame phase at the intron boundary. Source records retain approximate amino-acid rather than exact raw CDS-boundary coordinates, so unavailable raw nucleotide coordinates and unmapped reference positions are NA in the full data."
    )
    caption_s6 = CAPTION_S6
    (TABLES / "Table_S5_caption.txt").write_text(caption_s5 + "\n", encoding="utf-8")
    (TABLES / "Table_S6_caption.txt").write_text(caption_s6 + "\n", encoding="utf-8")

    s5a = [{"Family": display_family(r["family_id"]), "Eligible proteins": r["eligible_proteins"],
            "Proteins with introns": r["proteins_with_introns"], "Introns observed": r["total_introns"],
            "Protein prevalence (%)": r["protein_prevalence_pct"]} for r in intron_aux["part_a"]]
    s5b = [{"Family": display_family(r["family_id"]), "Reference position": r["reference_position_aa"],
            "Observations": r["observation_count"], "Proportion (%)": r["observation_proportion_pct"],
            "Phase": r["phase"], "Taxonomic breadth": r["taxonomic_breadth"], "Locus distribution": r["locus_distribution"]}
           for r in table_s5 if r["part"] == "B"]
    s6pdf = [{"Species": r["species"], "Life stage": r["life_stage"],
              "Tissue": r["tissue"] if r["tissue"] == "whole_body (in gall)" else r["tissue"].replace("_", " "),
              "Family": display_family(r["family_id"]), "Samples": r["sample_count"], "Median TPM": r["median_tpm"],
              "TPM IQR": r["iqr_tpm"], "TPM range": r["range_tpm"]} for r in table_s6]
    render_pdf("S5", f"Table {TABLE_LABELS['intron']}. Intron prevalence and recurrent positions", caption_s5, [
        ("Part A. Protein prevalence", ["Family", "Eligible proteins", "Proteins with introns", "Introns observed", "Protein prevalence (%)"],
         ["1.2in", "1.1in", "1.35in", "1.1in", "1.35in"], s5a),
        ("Part B. Recurrent reference positions", ["Family", "Reference position", "Observations", "Proportion (%)", "Phase", "Taxonomic breadth", "Locus distribution"],
         ["0.85in", "1.05in", "0.8in", "0.9in", "0.55in", "1.1in", "3.8in"], s5b),
    ])
    render_pdf("S6", f"Table {TABLE_LABELS['expression']}. Expression by biological context", caption_s6, [
        ("Expression summaries", ["Species", "Life stage", "Tissue", "Family", "Samples", "Median TPM", "TPM IQR", "TPM range"],
         ["1.65in", "0.7in", "1.55in", "0.65in", "0.5in", "0.75in", "1.15in", "1.35in"], s6pdf),
    ])

    qc = {"intron": intron_aux["qc"], "expression": expression_qc}
    # Contract gates.
    assert intron_aux["qc"]["b_lectin_alias_closure"] == 23
    assert len(sample_family) == 632
    assert all(n == 4 for n in Counter(r["sample_id"] for r in sample_family).values())
    assert len({r["sample_id"] for r in sample_meta}) == 158
    assert all(float(r["family_tpm"]) >= 0 for r in sample_family)
    assert all(float(r["tpm"]) >= 0 for r in gene_sample)
    assert len(gene_sample) == 4362
    assert len({(r["sample_id"], r["stable_target_sequence_id"]) for r in gene_sample}) == len(gene_sample)
    assert len(current_metadata_proteins) == 605
    assert sha256(DATA_I / "recurrent_position_summary.tsv") == sha256(current_intron / "recurrent_position_summary.tsv")
    assert sha256(TABLES / "Table_S5.tsv") == sha256(current_table_s5)
    assert not any(r["life_stage"] == "larva" and r["tissue"] == "gall" for r in table_s6)
    assert any(r["life_stage"] == "larva" and r["tissue"] == "whole_body (in gall)" for r in table_s6)
    ceroptres = [r for r in table_s6 if r["species"] == "Ceroptres sp. ex Druon ignotum" and r["family_id"] == "B_lectin"]
    assert len(ceroptres) == 1 and ceroptres[0]["sample_count"] == 6 and ceroptres[0]["range_tpm"] == "0.7-66.0"
    one_decimal = re.compile(r"^\d+\.\d$")
    interval = re.compile(r"^\d+\.\d-\d+\.\d$")
    assert all(one_decimal.fullmatch(r["median_tpm"]) for r in table_s6)
    assert all(interval.fullmatch(r["iqr_tpm"]) and interval.fullmatch(r["range_tpm"]) for r in table_s6)
    (EVIDENCE / "QC.json").write_text(json.dumps(qc, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    public_files = sorted(p for p in OUT.rglob("*") if p.is_file() and p.name != "MANIFEST.tsv" and p.suffix not in {".aux", ".log", ".out", ".tex"})
    forbidden = (b"DCAN", b"Dryocosmus cannoni", b"/" + b"carc/", b"/" + b"easley/", b"needs_review", b"pipeline_complete")
    for path in public_files:
        payload = path.read_bytes()
        if path.suffix != ".pdf" and any(token in payload for token in forbidden):
            raise AssertionError(f"forbidden public token in {path}")
    manifest = [{"relative_path": p.relative_to(OUT).as_posix(), "bytes": p.stat().st_size, "sha256": sha256(p),
                 "rows": sum(1 for _ in p.open(encoding="utf-8")) - 1 if p.suffix == ".tsv" else ""} for p in public_files]
    assert all(row["relative_path"] != "MANIFEST.tsv" for row in manifest)
    write_tsv(OUT / "MANIFEST.tsv", ["relative_path", "bytes", "sha256", "rows"], manifest)
    # Add manifest after its content is frozen; checksum list covers delivered artifacts except itself.
    write_tsv(INGRESS / "INGRESS_PROPOSAL.tsv", ["requirement_id", "candidate_relative_path", "sha256", "rows", "source_receipt_id"], [
        {"requirement_id": "REQ-INTRON-POSITIONS", "candidate_relative_path": "data/introns/recurrent_position_summary.tsv", "sha256": sha256(DATA_I / "recurrent_position_summary.tsv"), "rows": len(read_tsv(DATA_I / "recurrent_position_summary.tsv")), "source_receipt_id": "complete_current_roster_20260827"},
        *[{"requirement_id": f"REQ-INTRON-METADATA-{family}", "candidate_relative_path": f"data/introns/{family}_intron_metadata.tsv", "sha256": sha256(DATA_I / f"{family}_intron_metadata.tsv"), "rows": len(read_tsv(DATA_I / f"{family}_intron_metadata.tsv")), "source_receipt_id": "complete_current_roster_20260827"} for family in FAMILIES[:3]],
        {"requirement_id": "REQ-EXPRESSION-FAMILY-TPM", "candidate_relative_path": "data/expression/sample_family_tpm.tsv", "sha256": sha256(DATA_E / "sample_family_tpm.tsv"), "rows": len(sample_family), "source_receipt_id": "accepted_160_profile_baseline_plus_ceroptres_trinity_target"},
        {"requirement_id": "REQ-EXPRESSION-GENE-TPM", "candidate_relative_path": "data/expression/gene_sample_tpm.tsv", "sha256": sha256(DATA_E / "gene_sample_tpm.tsv"), "rows": len(gene_sample), "source_receipt_id": "accepted_160_profile_baseline_plus_six_row_target_receipt"},
        {"requirement_id": "REQ-EXPRESSION-SAMPLE-QC", "candidate_relative_path": "data/expression/sample_metadata_and_mapping_qc.tsv", "sha256": sha256(DATA_E / "sample_metadata_and_mapping_qc.tsv"), "rows": len(sample_meta), "source_receipt_id": "accepted_final_handoff_160_profile_metadata"},
    ])
    print(json.dumps(qc, sort_keys=True))


if __name__ == "__main__":
    main()
