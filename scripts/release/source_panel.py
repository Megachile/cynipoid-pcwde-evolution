"""Public metadata projection of the same-build complete source panel."""
import csv
from pathlib import Path

def build_source_panel(source: Path, source_panel: Path) -> None:
    if not source.is_file():
        raise SystemExit(f"missing same-build source-panel table: {source}")
    source_panel.parent.mkdir(parents=True, exist_ok=True)
    with source.open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    fields = list(rows[0])
    for row in rows:
        if row["species_display"] == "Mayetiola destructor" and row["public_accession"] == "GCA_000149185.1":
            row["assembly_busco_complete_pct"] = "95.7"
        if row["public_repository"] == "author_generated":
            accession = row["public_accession"]
            if "PRJNA" in accession:
                row["public_repository"] = "NCBI"
            elif "GSE" in accession:
                row["public_repository"] = "NCBI GEO"
            elif "USDA PBARC" in accession:
                row["public_repository"] = "USDA PBARC"
            else:
                row["public_repository"] = "Study dataset"
    with source_panel.open("w", newline="", encoding="utf-8") as fh:
        public_repo_index = fields.index("public_repository")
        public_accession_index = fields.index("public_accession")
        fields[public_repo_index] = "source_repository_or_study"
        fields[public_accession_index] = "source_accession_or_identifier"
        for row in rows:
            row["source_repository_or_study"] = row.pop("public_repository")
            row["source_accession_or_identifier"] = row.pop("public_accession")
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)
