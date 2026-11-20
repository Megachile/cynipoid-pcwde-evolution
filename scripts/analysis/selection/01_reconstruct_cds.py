#!/usr/bin/env python3
"""Reconstruct CDS for the locked current PCWDE protein membership.

Genome CDS are extracted from the DB-canonical assembly/GFF pair. Transcriptome
CDS are recovered from the corresponding ORF-mRNA artifact, including archived
rows when the paper intentionally retains a superseded transcriptome source.
Only CDS whose translation exactly matches the current paper protein is kept.
"""

from __future__ import annotations

import csv
import argparse
import re
from collections import defaultdict
from pathlib import Path
from urllib.parse import unquote

from Bio import SeqIO
from Bio.Seq import Seq
from pyfaidx import Fasta


FAMILIES = ("GH5_2", "PL1_8", "PL4_4")
MANIFEST_FIELDS = (
    "source", "source_file", "species", "species_code", "n_targets",
    "assembly_fasta", "annotation_gff", "orf_mrna_fastas",
    "provenance_method", "accepted_source_detail",
)

# These are display-name changes in the final paper metadata, not source
# resolution aliases. Source identity remains the stable (source, source_file)
# pair, while the frozen manifest names remain authoritative for replay audit
# output. Fail closed on every other metadata-to-manifest display difference.
ALLOWED_SPECIES_DISPLAY_ALIASES = {
    ("genome", "Ceroptres_masudai", "Ceroptres kovalevi"): "Ceroptres masudai",
    ("genome", "Eschatocerus_acacia", "Eschatocerus acaciae"): "Eschatocerus acacia",
    ("genome", "Isocolus_centaureae", "Eubothrus centaureae"): "Isocolus centaureae",
}


def read_fasta(path: Path) -> dict[str, str]:
    return {r.id: str(r.seq).upper() for r in SeqIO.parse(path, "fasta")}


def write_fasta(path: Path, records: dict[str, str]) -> None:
    with path.open("w") as handle:
        for name, seq in records.items():
            handle.write(f">{name}\n{seq}\n")


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def translate(cds: str) -> str:
    cds = cds.upper().replace("-", "").replace(" ", "")
    if len(cds) % 3:
        return ""
    aa = str(Seq(cds).translate())
    return aa[:-1] if aa.endswith("*") else aa


def clean_aa(aa: str) -> str:
    return aa.upper().rstrip("*")


def parse_attrs(raw: str) -> dict[str, list[str]]:
    attrs: dict[str, list[str]] = defaultdict(list)
    for item in raw.strip().split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
        elif " " in item:
            key, value = item.split(" ", 1)
            value = value.strip('"')
        else:
            continue
        attrs[key].extend(unquote(v) for v in value.split(","))
    return attrs


def id_forms(value: str) -> set[str]:
    forms = {value}
    for prefix in ("rna-", "cds-", "gene-", "transcript:", "protein:"):
        if value.startswith(prefix):
            forms.add(value[len(prefix):])
    expanded = set(forms)
    for form in forms:
        # Locked tandem/composite family segments use ``<mRNA>_<segment>``
        # identifiers, while Helixer GFF CDS rows retain only the mRNA parent.
        # Preserve the segment identity as the output key, but allow its exact
        # unsuffixed parent as a lookup alias.  Restrict this to an mRNA-style
        # terminal ``.<isoform>_<segment>`` suffix so ordinary underscore IDs
        # are not broadened accidentally.
        expanded.add(re.sub(r"(\.\d+)_\d+$", r"\1", form))
        expanded.add(re.sub(r"\.gene$", "", form))
        expanded.add(re.sub(r"\.t\d+$", "", form))
        expanded.add(re.sub(r"\.p\d+$", "", form))
        expanded.add(re.sub(r"\.\d+$", "", form))
        if "pseudo_chr_" in form:
            expanded.add(form[form.index("pseudo_chr_"):])
        if "TRINITY_" in form:
            core = form[form.index("TRINITY_"):]
            expanded.add(core)
            expanded.add(re.sub(r"\.p\d+$", "", core))
    return {x for x in expanded if x}


def find_cds_in_transcript(nt: str, aa: str) -> set[str]:
    matches: set[str] = set()
    for oriented in (nt.upper(), str(Seq(nt).reverse_complement()).upper()):
        for frame in range(3):
            usable = oriented[frame:]
            usable = usable[:len(usable) - len(usable) % 3]
            translated_seq = str(Seq(usable).translate())
            start = translated_seq.find(aa)
            while start >= 0:
                cds_start = frame + 3 * start
                cds = oriented[cds_start:cds_start + 3 * len(aa)]
                if translate(cds) == aa:
                    matches.add(cds)
                start = translated_seq.find(aa, start + 1)
    return matches


def transcript_matches(paths: list[Path], targets: dict[str, str]) -> tuple[dict[str, set[str]], dict[str, str]]:
    matches: dict[str, set[str]] = defaultdict(set)
    matched_path: dict[str, str] = {}
    form_to_target: dict[str, set[str]] = defaultdict(set)
    for target in targets:
        for form in id_forms(target):
            form_to_target[form].add(target)
    for path in paths:
        for record in SeqIO.parse(path, "fasta"):
            candidates: set[str] = set()
            for form in id_forms(record.id):
                candidates.update(form_to_target.get(form, set()))
            for target in candidates:
                found = find_cds_in_transcript(str(record.seq), targets[target])
                if found:
                    matches[target].update(found)
                    matched_path[target] = str(path)
    return matches, matched_path


def load_source_manifest(path: Path) -> dict[tuple[str, str], dict[str, object]]:
    """Load the frozen source resolution used by the accepted analysis."""
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    if not rows or tuple(reader.fieldnames or ()) != MANIFEST_FIELDS:
        raise ValueError(f"{path}: expected exact columns {list(MANIFEST_FIELDS)}")

    def resolve(value: str) -> Path:
        candidate = Path(value)
        return candidate if candidate.is_absolute() else (path.parent / candidate).resolve()

    result: dict[tuple[str, str], dict[str, object]] = {}
    for row in rows:
        key = (row["source"], row["source_file"])
        if key in result:
            raise ValueError(f"{path}: duplicate source row {key}")
        result[key] = {
            **row,
            "n_targets": int(row["n_targets"]),
            "assembly_fasta": resolve(row["assembly_fasta"]) if row["assembly_fasta"] else None,
            "annotation_gff": resolve(row["annotation_gff"]) if row["annotation_gff"] else None,
            "orf_mrna_fastas": [resolve(value) for value in row["orf_mrna_fastas"].split("|") if value],
        }
    return result


def validate_source_groups(
    source_manifest: dict[tuple[str, str], dict[str, object]],
    groups: dict[tuple[str, str], list[dict]],
) -> list[dict[str, str]]:
    """Validate stable source keys and approved metadata display differences."""
    missing_sources = sorted(set(groups) - set(source_manifest))
    unused_sources = sorted(set(source_manifest) - set(groups))
    if missing_sources or unused_sources:
        raise RuntimeError(
            f"source manifest mismatch: missing={missing_sources}; unused={unused_sources}"
        )

    aliases: list[dict[str, str]] = []
    for key in sorted(groups):
        source, source_file = key
        display_names = sorted({row["species"] for row in groups[key]})
        if len(display_names) != 1:
            raise RuntimeError(
                f"{key}: expected one metadata species display, found {display_names}"
            )
        metadata_species = display_names[0]
        frozen_species = str(source_manifest[key]["species"])
        if metadata_species == frozen_species:
            continue
        approved = ALLOWED_SPECIES_DISPLAY_ALIASES.get(
            (source, source_file, metadata_species)
        )
        if approved != frozen_species:
            raise RuntimeError(
                f"{key}: unapproved metadata-to-manifest species display difference "
                f"{metadata_species!r} -> {frozen_species!r}"
            )
        aliases.append({
            "source": source,
            "source_file": source_file,
            "metadata_species": metadata_species,
            "frozen_audit_species": frozen_species,
        })
    return aliases


def parse_target_cds(gff: Path, targets: list[str]) -> dict[str, dict]:
    form_to_target: dict[str, set[str]] = defaultdict(set)
    for target in targets:
        for form in id_forms(target):
            form_to_target[form].add(target)
    found: dict[str, dict] = {}
    with gff.open(errors="replace") as handle:
        for raw in handle:
            if not raw or raw.startswith("#"):
                continue
            parts = raw.rstrip("\n").split("\t")
            if len(parts) != 9 or parts[2].lower() != "cds":
                continue
            attrs = parse_attrs(parts[8])
            candidates: set[str] = set()
            for key in ("protein_id", "Parent", "ID", "transcript_id"):
                for value in attrs.get(key, []):
                    for form in id_forms(value):
                        candidates.update(form_to_target.get(form, set()))
            if len(candidates) != 1:
                continue
            target = next(iter(candidates))
            entry = found.setdefault(target, {"strand": parts[6], "features": []})
            if entry["strand"] != parts[6]:
                entry["conflict"] = True
            entry["features"].append((parts[0], int(parts[3]), int(parts[4]), parts[7]))
    return found


def extract_genome_cds(
    assembly: Path, gff: Path, targets: list[str], target_aa: dict[str, str], index_name: Path
) -> tuple[dict[str, str], dict[str, str]]:
    features = parse_target_cds(gff, targets)
    index_name.parent.mkdir(parents=True, exist_ok=True)
    index_name.unlink(missing_ok=True)
    # The same source label can resolve to a newer DB assembly or the frozen
    # paper snapshot. Always rebuild the analysis-local index to avoid reusing
    # offsets from a different FASTA after provenance resolution changes.
    genome = Fasta(str(assembly), indexname=str(index_name), rebuild=True)
    cds: dict[str, str] = {}
    errors: dict[str, str] = {}
    try:
        for target in targets:
            entry = features.get(target)
            if not entry:
                errors[target] = "target_not_found_in_gff"
                continue
            if entry.get("conflict"):
                errors[target] = "conflicting_scaffold_or_strand"
                continue
            reverse = entry["strand"] == "-"
            by_seqid: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
            seqid_order: list[str] = []
            for seqid, start, end, phase in entry["features"]:
                if seqid not in by_seqid:
                    seqid_order.append(seqid)
                by_seqid[seqid].append((start, end, phase))
            candidates = []
            block_orders = [seqid_order]
            if len(seqid_order) > 1:
                block_orders.append(list(reversed(seqid_order)))
            for block_order in block_orders:
                chunks = []
                for seqid in block_order:
                    ordered = sorted(by_seqid[seqid], key=lambda x: x[0], reverse=reverse)
                    for start, end, _phase in ordered:
                        seq = str(genome[seqid][start - 1:end]).upper()
                        chunks.append(str(Seq(seq).reverse_complement()) if reverse else seq)
                candidates.append("".join(chunks))
            exact_segments: set[str] = set()
            for candidate in candidates:
                exact_segments.update(find_cds_in_transcript(candidate, target_aa[target]))
            if len(exact_segments) == 1:
                cds[target] = next(iter(exact_segments))
            else:
                cds[target] = candidates[0]
                if len(exact_segments) > 1:
                    errors[target] = "ambiguous_exact_cds_segments"
    finally:
        genome.close()
    return cds, errors


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Recover exact CDS for the accepted PCWDE protein membership."
    )
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument(
        "--protein-dir", required=True, type=Path,
        help="Directory containing GH5_2.faa, PL1_8.faa, and PL4_4.faa.",
    )
    parser.add_argument(
        "--source-manifest", required=True, type=Path,
        help="Frozen genome/transcriptome source-resolution manifest.",
    )
    parser.add_argument(
        "--display-alias-receipt", type=Path,
        help="Optional TSV receipt for approved metadata-to-frozen display differences.",
    )
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=False)
    metadata = list(csv.DictReader(args.metadata.open(), delimiter="\t"))
    proteins = {family: read_fasta(args.protein_dir / f"{family}.faa") for family in FAMILIES}
    source_manifest = load_source_manifest(args.source_manifest)
    rows = [
        row for row in metadata
        if row["family"] in FAMILIES and row["protein_id"] in proteins[row["family"]]
    ]
    groups: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for row in rows:
        groups[(row["source"], row["source_file"])].append(row)
    aliases = validate_source_groups(source_manifest, groups)
    if args.display_alias_receipt:
        write_tsv(
            args.display_alias_receipt, aliases,
            ["source", "source_file", "metadata_species", "frozen_audit_species"],
        )

    retained: dict[str, dict[str, str]] = {family: {} for family in FAMILIES}
    audit: list[dict] = []
    source_audit: list[dict] = []
    for group_number, (group_key, members) in enumerate(sorted(groups.items()), 1):
        source, source_file = group_key
        manifest = source_manifest[group_key]
        species = str(manifest["species"])
        code = str(manifest["species_code"])
        targets = [row["protein_id"] for row in members]
        if int(manifest["n_targets"]) != len(targets):
            raise RuntimeError(
                f"{source_file}: manifest target count {manifest['n_targets']} != {len(targets)}"
            )
        target_aa = {
            row["protein_id"]: clean_aa(proteins[row["family"]][row["protein_id"]])
            for row in members
        }
        recovered: dict[str, set[str]] = defaultdict(set)
        recovery_path: dict[str, str] = {}
        recovery_reason: dict[str, str] = {}
        source_status = "ok"
        source_detail = str(manifest["accepted_source_detail"])
        try:
            if source == "transcriptome":
                paths = list(manifest["orf_mrna_fastas"])
                if not paths or not all(path.is_file() for path in paths):
                    raise FileNotFoundError(f"missing ORF-mRNA input for {source_file}")
                recovered, recovery_path = transcript_matches(paths, target_aa)
            elif source == "genome":
                assembly = manifest["assembly_fasta"]
                gff = manifest["annotation_gff"]
                if not isinstance(assembly, Path) or not isinstance(gff, Path):
                    raise FileNotFoundError(f"missing assembly/GFF manifest values for {source_file}")
                mapping_detail = re.sub(r";method=[^;]+$", "", source_detail)
                extracted, extraction_errors = extract_genome_cds(
                    assembly, gff, targets, target_aa,
                    out / "fai" / f"{group_number:03d}_{source_file}.fai"
                )
                for target, cds in extracted.items():
                    if translate(cds) == target_aa[target]:
                        recovered[target].add(cds)
                        recovery_path[target] = mapping_detail
                    elif target not in extraction_errors:
                        extraction_errors[target] = "translation_mismatch"
                recovery_reason.update(extraction_errors)
            else:
                raise ValueError(f"unsupported source={source}")
        except Exception as exc:
            source_status = "failed"
            source_detail = f"{type(exc).__name__}: {exc}"
            recovery_reason.update({target: "source_failed" for target in targets})

        source_audit.append({
            "source": source, "source_file": source_file, "species": species,
            "species_code": code, "n_targets": len(targets),
            "n_exact_cds": sum(len(v) == 1 for v in recovered.values()),
            "status": source_status, "source_detail": source_detail,
        })
        print(
            f"[{group_number}/{len(groups)}] {source_file}: "
            f"{sum(len(v)==1 for v in recovered.values())}/{len(targets)}"
        )

        for row in members:
            target = row["protein_id"]
            choices = recovered.get(target, set())
            if len(choices) == 1:
                cds = next(iter(choices))
                retained[row["family"]][target] = cds
                status = "mapped_exact_translation"
                reason = "exact_translation"
            elif not choices:
                status = "unmatched"
                reason = recovery_reason.get(target, "no_exact_translation_match")
            else:
                status = "ambiguous_distinct_cds"
                reason = "multiple_distinct_exact_cds"
            audit.append({
                "family": row["family"], "protein_id": target, "species": species,
                "source": source, "source_file": source_file, "species_code": code,
                "aa_length": len(target_aa[target]), "mapping_status": status,
                "mapping_reason": reason, "n_distinct_cds": len(choices),
                "cds_path": recovery_path.get(target, ""),
            })

    for family in FAMILIES:
        write_fasta(out / f"{family}_locked_current_cds.fna", retained[family])
    write_tsv(
        out / "current_cds_mapping.tsv", audit,
        ["family", "protein_id", "species", "source", "source_file", "species_code", "aa_length",
         "mapping_status", "mapping_reason", "n_distinct_cds", "cds_path"],
    )
    write_tsv(
        out / "current_cds_source_audit.tsv", source_audit,
        ["source", "source_file", "species", "species_code", "n_targets", "n_exact_cds", "status", "source_detail"],
    )
    dl_rows = [row for row in source_audit if row["source_file"] == "DL_bkragtag"]
    if len(dl_rows) != 1:
        raise SystemExit(f"DL_bkragtag: expected one source-audit row, found {len(dl_rows)}")
    dl = dl_rows[0]
    if dl["status"] != "ok" or dl["n_targets"] != 30 or dl["n_exact_cds"] != 30:
        raise SystemExit(
            "DL_bkragtag: fail-closed source gate requires 30/30 exact CDS "
            f"with status=ok; observed {dl['n_exact_cds']}/{dl['n_targets']} status={dl['status']}"
        )
    for family in FAMILIES:
        total = sum(row["family"] == family for row in rows)
        print(f"SUMMARY {family}: {len(retained[family])}/{total} exact translation matches")


if __name__ == "__main__":
    main()
