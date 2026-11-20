#!/usr/bin/env python3
"""Build the isolated selection/catalytic table package."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

import fitz


FAMILIES = ("GH5_2", "PL1_8", "PL4_4")
with (Path(__file__).resolve().parents[2] / "config/table_numbers.tsv").open(newline="") as handle:
    SELECTION_LABEL = {r["table_id"]:r["label"] for r in csv.DictReader(handle,delimiter="\t")}["selection_catalytic"]
S1_FIELDS = (
    "resource_class", "species_display", "source_accession",
    "resource_preparation", "panel_status_exclusion_reason",
    "assembly_n50_bp", "assembly_busco_pct", "proteome_busco_pct",
)
SOURCE_TABLE_S1_ALIASES = {
    # Historical selection labels joined to accepted current Table S1 displays.
    "Ceroptres masudai": "Ceroptres kovalevi",
    "Eschatocerus acacia": "Eschatocerus acaciae",
    "Isocolus centaureae": "Eubothrus centaureae",
    "Ceroptres (DI inquiline)": "Ceroptres sp. ex Druon ignotum",
    "Synergus (BK7 inquiline)": "Synergus sp. ex Belonocnema kinseyi",
    "Synergus (DI inquiline)": "Synergus sp. ex Druon ignotum",
    "Synergus sp. (D. quercuslanigerum inquiline)": "Synergus sp. ex Druon quercuslanigerum",
}
PUBLIC_SPECIES_ALIASES = {
    **{key: value for key, value in SOURCE_TABLE_S1_ALIASES.items()
       if "inquiline" not in key},
    "Ceroptres (DI inquiline)": "Ceroptres sp. ex Druon ignotum",
    "Synergus (BK7 inquiline)": "Synergus sp. ex Belonocnema kinseyi",
    "Synergus (DI inquiline)": "Synergus sp. ex Druon ignotum",
    "Synergus sp. (D. quercuslanigerum inquiline)": "Synergus sp. ex Druon quercuslanigerum",
}


def read_tsv(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader), list(reader.fieldnames or [])


def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({field: row.get(field, "") for field in fields} for row in rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def source_kind(resource_class: str) -> str:
    if resource_class == "Genome":
        return "genome"
    if resource_class == "Transcriptome-derived proteome":
        return "transcriptome"
    return "reference"


def tex_escape(value: object) -> str:
    text = str(value)
    replacements = {
        "\\": r"\textbackslash{}", "&": r"\&", "%": r"\%", "$": r"\$",
        "#": r"\#", "_": r"\_", "{": r"\{", "}": r"\}",
        "~": r"\textasciitilde{}", "^": r"\textasciicircum{}",
    }
    return "".join(replacements.get(char, char) for char in text)


def table_tex(headers: list[str], widths: list[str], data: list[list[str]]) -> list[str]:
    colspec = "@{}" + "".join(
        rf">{{\centering\arraybackslash}}p{{{width}}}" for width in widths
    ) + "@{}"
    header = " & ".join(rf"\textbf{{{item}}}" for item in headers) + r" \\ \midrule"
    lines = [rf"\scriptsize\begin{{longtable}}{{{colspec}}}", r"\toprule", header,
             r"\endfirsthead", r"\toprule", header, r"\endhead",
             r"\bottomrule", r"\endfoot"]
    for row in data:
        lines.append(" & ".join(tex_escape(value) for value in row) + r" \\")
    lines.append(r"\end{longtable}")
    return lines


def render_reader_pdf(path: Path, selection_rows: list[dict[str, str]],
                      catalytic_rows: list[dict[str, str]]) -> None:
    selection_data = [[row[key] for key in (
        "gene_family", "analyzed_sequences", "retained_codon_sites", "purifying_display",
        "diversifying_display", "median_omega_display", "terminal_display")]
        for row in selection_rows]
    catalytic_data = [[row[key] for key in (
        "gene_family", "site_role", "reference_residue", "observable_total_display",
        "expected_observable_display", "expected_observable_pct")]
        for row in catalytic_rows]
    lines = [
        r"\documentclass[10pt]{article}",
        r"\usepackage[margin=0.45in,landscape]{geometry}",
        r"\usepackage{longtable,booktabs,array,xcolor}",
        r"\setlength{\tabcolsep}{3pt}", r"\renewcommand{\arraystretch}{1.12}",
        r"\pagestyle{plain}", r"\begin{document}",
        rf"\section*{{Table {SELECTION_LABEL}. Selection and catalytic-residue conservation across three PCWDE families.}}",
        (r"\small Part A summarizes family-level selection and alignment results. Part B reports "
         r"conservation at reference catalytic residues.\par\medskip"),
        r"\subsection*{Part A. Selection and alignment results by family}",
    ]
    lines.extend(table_tex(
        ["Family", "Analyzed sequences", "Codon sites", r"Purifying sites, n (\%)",
         r"Diversifying sites, n (\%)", r"Median $\omega$", r"Significant/tested terminal branches, n/N (\%)"],
        ["0.78in", "1.0in", "0.82in", "1.35in", "1.0in", "0.78in", "2.45in"], selection_data,
    ))
    lines.append(r"\subsection*{Part B. Catalytic-residue conservation by site}")
    lines.extend(table_tex(
        ["Family", "Site / role", "Reference residue", "Observable/total",
         "Expected/observable", r"Expected among observable (\%)"],
        ["0.9in", "1.65in", "1.15in", "1.25in", "1.45in", "1.7in"], catalytic_data,
    ))
    lines.extend([
        r"\footnotesize\textit{Notes.} FUBAR classes use posterior probability $>0.9$. Median $\omega$ is $\beta/\alpha$ for sites with $\alpha>0.001$.\par",
        r"aBSREL significance uses the Holm-corrected $P<0.05$ threshold; n/N is significant/tested terminal branches.\par",
        r"Catalytic percentages use observable residues as the denominator.\par",
        r"\end{document}",
    ])
    with tempfile.TemporaryDirectory(prefix="selection_table_tex_") as temp_name:
        temp = Path(temp_name)
        tex_path = temp / "selection_catalytic_reader_table.tex"
        tex_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        tex_env = os.environ.copy()
        tex_env["SOURCE_DATE_EPOCH"] = "1786665600"
        tex_env["FORCE_SOURCE_DATE"] = "1"
        for _ in range(2):
            subprocess.run(
                ["xelatex", "-interaction=nonstopmode", "-halt-on-error", tex_path.name],
                cwd=temp, check=True, stdout=subprocess.DEVNULL, env=tex_env,
            )
        rendered = temp / "selection_catalytic_reader_table.pdf"
        if not rendered.is_file() or rendered.stat().st_size == 0:
            raise RuntimeError("XeLaTeX did not produce the reader PDF")
        # Strip engine/process metadata while retaining only publication-facing metadata.
        with fitz.open(rendered) as document:
            document.set_metadata({
                "title": f"Table {SELECTION_LABEL}. Selection and catalytic-residue conservation across three PCWDE families.",
                "subject": "Selection and catalytic-residue conservation",
                "author": "",
                "keywords": "",
                "creator": "",
                "producer": "",
                "creationDate": "",
                "modDate": "",
                "trapped": "",
            })
            document.save(path, garbage=4, deflate=True, no_new_id=True)


def render_native_preview(pdf_path: Path, png_path: Path) -> None:
    with fitz.open(pdf_path) as document:
        if len(document) != 1:
            raise RuntimeError("reader PDF unexpectedly changed page count")
        pixmap = document[0].get_pixmap(matrix=fitz.Matrix(2, 2), alpha=False)
        png_path.write_bytes(pixmap.tobytes("png"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--selection", type=Path, required=True)
    parser.add_argument("--catalytic", type=Path, required=True)
    parser.add_argument("--table-s1", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise SystemExit(f"refusing to overwrite existing output directory: {args.out}")
    args.out.mkdir(parents=True)

    source_files = [
        args.selection / "selection_family_summary.tsv",
        args.selection / "fubar_site_results.tsv",
        args.selection / "absrel_branch_results.tsv",
        args.selection / "alignment_qc.tsv",
        args.selection / "selection_membership_lock.tsv",
        args.selection / "current_cds_mapping.tsv",
        args.catalytic / "catalytic_site_conservation_summary.tsv",
        args.catalytic / "catalytic_site_conservation_by_sequence.tsv",
        args.table_s1,
    ]
    for path in source_files:
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(f"missing or empty accepted input: {path}")

    s1_rows, s1_fields = read_tsv(args.table_s1)
    if tuple(s1_fields) != S1_FIELDS or len(s1_rows) != 99:
        raise RuntimeError(f"accepted Table S1 contract changed: fields={s1_fields}, rows={len(s1_rows)}")
    s1_by_key = {(r["species_display"], source_kind(r["resource_class"])): r for r in s1_rows}

    summary, _ = read_tsv(args.selection / "selection_family_summary.tsv")
    catalytic_summary, _ = read_tsv(args.catalytic / "catalytic_site_conservation_summary.tsv")
    catalytic_primary = [r for r in catalytic_summary
                         if r["alignment_set"] == "de_novo" and r["site_role"] == "catalytic"]
    selection_reader_rows = []
    for row in summary:
        purifying_pct = 100 * int(row["purifying"]) / int(row["sites"])
        diversifying_pct = 100 * int(row["positive"]) / int(row["sites"])
        terminal_pct = 100 * int(row["absrel_significant_terminal"]) / int(row["absrel_terminal_tested"])
        selection_reader_rows.append({
            "gene_family": row["family"], "analyzed_sequences": row["selection_sequences"],
            "retained_codon_sites": row["sites"], "purifying_sites": row["purifying"],
            "purifying_percent": f"{purifying_pct:.1f}", "diversifying_sites": row["positive"],
            "diversifying_percent": f"{diversifying_pct:.1f}",
            "median_omega": row["median_omega"], "median_omega_display": f"{float(row['median_omega']):.3f}",
            "significant_terminal_branches": row["absrel_significant_terminal"],
            "terminal_branches_tested": row["absrel_terminal_tested"],
            "significant_terminal_percent": f"{terminal_pct:.1f}",
            "purifying_display": f"{row['purifying']} ({purifying_pct:.1f}%)",
            "diversifying_display": f"{row['positive']} ({diversifying_pct:.1f}%)",
            "terminal_display": f"{row['absrel_significant_terminal']}/{row['absrel_terminal_tested']} ({terminal_pct:.1f}%)",
        })
    catalytic_reader_rows = [{
        "gene_family": row["family"], "site_role": row["site_label"].replace("_", " "),
        "reference_residue": row["reference_residue"],
        "observable_sequences": row["n_site_observed"], "total_sequences": row["n_total_sequences"],
        "expected_residues": row["n_expected_residue"],
        "expected_observable_pct": f"{float(row['pct_expected_among_observed']):.1f}",
        "observable_total_display": f"{row['n_site_observed']}/{row['n_total_sequences']}",
        "expected_observable_display": f"{row['n_expected_residue']}/{row['n_site_observed']}",
    } for row in catalytic_primary]
    if [r["gene_family"] for r in selection_reader_rows] != list(FAMILIES) or len(catalytic_reader_rows) != 5:
        raise RuntimeError("family or catalytic reader dimensions changed")
    write_tsv(args.out / "selection_reader_table.tsv", selection_reader_rows, list(selection_reader_rows[0]))
    write_tsv(args.out / "catalytic_reader_table.tsv", catalytic_reader_rows, list(catalytic_reader_rows[0]))
    pdf_path = args.out / "selection_catalytic_reader_table.pdf"
    render_reader_pdf(pdf_path, selection_reader_rows, catalytic_reader_rows)
    render_native_preview(pdf_path, args.out / "selection_catalytic_reader_table_page_1.png")

    copy_specs = [
        (args.selection / "selection_family_summary.tsv", "selection_family_summary.tsv"),
        (args.selection / "fubar_site_results.tsv", "fubar_site_evidence.tsv"),
        (args.selection / "absrel_branch_results.tsv", "absrel_branch_evidence.tsv"),
        (args.selection / "alignment_qc.tsv", "selection_alignment_qc.tsv"),
    ]
    for source, name in copy_specs:
        shutil.copyfile(source, args.out / name)

    catalytic_fields = list(catalytic_primary[0])
    write_tsv(args.out / "catalytic_residue_summary.tsv", catalytic_primary, catalytic_fields)
    catalytic_detail, catalytic_detail_fields = read_tsv(args.catalytic / "catalytic_site_conservation_by_sequence.tsv")
    catalytic_detail = [r for r in catalytic_detail
                        if r["alignment_set"] == "de_novo" and r["site_role"] == "catalytic"]
    write_tsv(args.out / "catalytic_residue_by_sequence.tsv", catalytic_detail, catalytic_detail_fields)

    members_all, _ = read_tsv(args.selection / "selection_membership_lock.tsv")
    members = [row for row in members_all if row["selection_status"] == "included"]
    member_counts = {family: sum(row["family"] == family for row in members) for family in FAMILIES}
    expected_member_counts = {
        row["family"]: int(row["selection_sequences"]) for row in summary
    }
    if member_counts != expected_member_counts or len(members) != 911:
        raise RuntimeError(
            f"included membership does not match accepted selection terminals: "
            f"{member_counts} != {expected_member_counts}"
        )
    cds_rows, _ = read_tsv(args.selection / "current_cds_mapping.tsv")
    cds_map = {(r["family"], r["protein_id"]): r for r in cds_rows}
    member_fields = [
        "gene_family", "protein_id", "species", "source", "analysis_id", "selection_status",
        "cds_recovery_method", "source_accession", "resource_preparation", "panel_status",
    ]
    public_members = []
    missing_s1 = set()
    for row in members:
        source_species = SOURCE_TABLE_S1_ALIASES.get(row["species"], row["species"])
        public_species = PUBLIC_SPECIES_ALIASES.get(row["species"], row["species"])
        s1 = s1_by_key.get((source_species, row["source"]))
        if s1 is None:
            missing_s1.add((row["species"], source_species, row["source"]))
            s1 = {}
        mapped = cds_map.get((row["family"], row["protein_id"]), {})
        if mapped.get("mapping_status") != "mapped_exact_translation":
            raise RuntimeError(
                f"included selection member lacks exact CDS mapping: "
                f"{row['family']} {row['protein_id']}"
            )
        public_members.append({
            "gene_family": row["family"], "protein_id": row["protein_id"],
            "species": public_species, "source": row["source"], "analysis_id": row["analysis_id"],
            "selection_status": row["selection_status"],
            "cds_recovery_method": mapped.get("mapping_status", ""),
            "source_accession": s1.get("source_accession", ""),
            "resource_preparation": s1.get("resource_preparation", ""),
            "panel_status": s1.get("panel_status_exclusion_reason", ""),
        })
    if missing_s1:
        raise RuntimeError(f"membership rows absent from accepted Table S1: {sorted(missing_s1)}")
    write_tsv(args.out / "selection_membership_and_public_source.tsv", public_members, member_fields)

    authority_dir = args.out / "authority_snapshot"
    authority_dir.mkdir()
    shutil.copyfile(args.table_s1, authority_dir / "Table_S1.tsv")
    if sha256(args.table_s1) != sha256(authority_dir / "Table_S1.tsv"):
        raise RuntimeError("Table S1 snapshot is not byte-identical")

    audit = {
        "decision": "Full source-panel dependency removed; only accepted public Table S1 fields are joined.",
        "accepted_table_s1_fields": list(S1_FIELDS),
        "joined_membership_fields": ["source_accession", "resource_preparation", "panel_status_exclusion_reason"],
        "public_species_serialization_aliases": PUBLIC_SPECIES_ALIASES,
        "frozen_table_s1_join_aliases": SOURCE_TABLE_S1_ALIASES,
        "source_label_discrepancy": (
            "The accepted selection membership retains historical inquiline labels. Those labels "
            "join to accepted current Table S1 displays, and every derived public selection file "
            "uses the established 'Genus sp. ex Host species' policy."
        ),
        "excluded_dependencies": ["paper_source_expression_refresh/outputs", "AU-test artifacts", "unrelated supplement tables"],
        "table_s1_rows": len(s1_rows), "membership_rows": len(public_members),
        "excluded_no_exact_cds": len(members_all) - len(members), "unmatched_membership_sources": 0,
    }
    (args.out / "source_panel_dependency_audit.json").write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    readme = """# Selection and catalytic-residue conservation across three PCWDE families

`selection_catalytic_reader_table.pdf`, `selection_reader_table.tsv`, and `catalytic_reader_table.tsv` provide the concise reader table. The remaining TSV files provide family summaries and elaborated site, branch, catalytic-residue, alignment-QC, and membership evidence. `Table_S1.tsv` supplies the accession, preparation, and panel-status fields used in the membership table.

FUBAR classes use posterior probability > 0.9; purifying and diversifying percentages use retained codon sites for the same family as the denominator. Median omega is beta/alpha for sites with alpha > 0.001. aBSREL significance uses the Holm-corrected P < 0.05 threshold. Catalytic percentages use observable residues as the denominator.
"""
    (args.out / "README.md").write_text(readme, encoding="utf-8")

    input_manifest = [{"path": str(p), "sha256": sha256(p), "size_bytes": str(p.stat().st_size)}
                      for p in source_files]
    write_tsv(args.out / "INPUT_FREEZE.tsv", input_manifest, ["path", "sha256", "size_bytes"])

    artifact_paths = sorted(p for p in args.out.rglob("*") if p.is_file())
    manifest = []
    for path in artifact_paths:
        rel = path.relative_to(args.out).as_posix()
        rows = ""
        if path.suffix == ".tsv":
            rows = str(max(0, sum(1 for _ in path.open(encoding="utf-8")) - 1))
        manifest.append({"file": rel, "sha256": sha256(path), "size_bytes": str(path.stat().st_size), "data_rows": rows})
    write_tsv(args.out / "MANIFEST.tsv", manifest, ["file", "sha256", "size_bytes", "data_rows"])



if __name__ == "__main__":
    main()
