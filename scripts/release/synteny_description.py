"""Reader description from the table producer's synteny summary."""
import csv


def workbook_description(summary_path):
    with summary_path.open(newline="", encoding="utf-8") as handle:
        facts = {row["metric"]: int(row["value"])
                 for row in csv.DictReader(handle, delimiter="\t") if row["family"] == "all"}
    return (f"Workbook indexing all {facts['family_specific_blocks']} conserved blocks "
            f"and showing ordered flanking genes for each of the "
            f"{facts['cross_genome_detectable_patterns']} blocks detected in multiple genomes.")
