"""Project the frozen reader matrix through the shared retained-sample crosswalk."""
import csv
from pathlib import Path

FIELDS = ("full_sample_id", "original_sample_id", "sample_id", "species_code", "species_id", "species_display")


def project_reader_expression(crosswalk_path: Path, source_path: Path, output_path: Path) -> int:
    if not crosswalk_path.is_file():
        raise ValueError(f"missing retained-sample crosswalk: {crosswalk_path}")
    with crosswalk_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != list(FIELDS):
            raise ValueError("retained-sample crosswalk schema mismatch")
        crosswalk = list(reader)
    if len(crosswalk) != 158 or any(not all(row[k] for k in FIELDS) for row in crosswalk):
        raise ValueError("retained-sample crosswalk has missing or unexpected profiles")
    for key in ("full_sample_id", "sample_id"):
        if len({row[key] for row in crosswalk}) != len(crosswalk):
            raise ValueError(f"duplicate retained-sample crosswalk key: {key}")
    keys = [(row["species_display"], row["original_sample_id"]) for row in crosswalk]
    if len(set(keys)) != len(keys):
        raise ValueError("ambiguous species/sample join in retained-sample crosswalk")
    retained = set(keys)
    with source_path.open(newline="", encoding="utf-8") as source:
        reader = csv.DictReader(source, delimiter="\t")
        if not reader.fieldnames or not {"species", "sample"} <= set(reader.fieldnames):
            raise ValueError("reader expression schema lacks species/sample")
        rows = list(reader)
        known = {(row["species"], row["sample"]) for row in rows}
        if not retained <= known:
            raise ValueError("retained crosswalk samples absent from reader matrix")
        selected = [row for row in rows if (row["species"], row["sample"]) in retained]
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with output_path.open("w", newline="", encoding="utf-8") as target:
            writer = csv.DictWriter(target, fieldnames=reader.fieldnames, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(selected)
    return len(selected)
