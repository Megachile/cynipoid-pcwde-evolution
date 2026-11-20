#!/usr/bin/env python3
"""Build Table S6 (expression by biological context) and its expression intermediates only.

Uses build_table_s6.build_expression(); runs no intron step and writes no Table S5,
so nothing saved is copied forward. Method: scripts/tables/README.md.

Environment: project Python 3 (standard library only).
Inputs: PCWDE_REPOSITORY_INPUTS (deposit repository_inputs); output under PCWDE_BUILD_ROOT/table_s6.
"""
import json
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import build_table_s6 as s6  # noqa: E402

OUT = Path(os.environ["PCWDE_BUILD_ROOT"]) / "table_s6"
if OUT.exists():
    raise SystemExit(f"refusing existing output: {OUT}")
tables, expression = OUT / "supplementary_tables", OUT / "data/expression"
tables.mkdir(parents=True)
expression.mkdir(parents=True)
sample_family, gene_sample, sample_meta, table_s6, crosswalk, qc = s6.build_expression()
source_panel = s6.read_tsv(s6.AUTH_S2 / "source_panel_accepted.tsv")
sample_counts = s6.sample_count_matrix(sample_family, sample_meta, crosswalk, table_s6, source_panel)
s6.write_tsv(expression / "sample_family_tpm.tsv",
             ["sample_family_id", "sample_id", "family_id", "family_tpm", "quantification_basis"], sample_family)
s6.write_tsv(expression / "gene_sample_tpm.tsv",
             ["gene_sample_id", "sample_id", "stable_target_sequence_id", "parent_gene_or_locus_id", "transcript_sequence_id", "tpm"],
             gene_sample)
s6.write_tsv(expression / "sample_metadata_and_mapping_qc.tsv",
             ["sample_id", "species_id", "life_stage", "tissue", "mapping_rate_pct", "included_in_160_profile_scope", "qc_note"],
             sample_meta)
s6.write_tsv(expression / "retained_sample_crosswalk.tsv", s6.RETAINED_CROSSWALK_FIELDS, crosswalk)
s6.write_tsv(tables / "Table_S6.tsv",
             ["species", "life_stage", "tissue", "family_id", "sample_count", "median_tpm", "iqr_tpm", "range_tpm"], table_s6)
s6.write_tsv(expression / "sample_count_matrix.tsv", s6.SAMPLE_COUNT_FIELDS, sample_counts)
(tables / "Table_S6a_caption.txt").write_text(
    "Table S2a. Cells give the number of retained expression profiles for each species and life-stage/tissue context; blank cells indicate zero.\n",
    encoding="utf-8")
(tables / "Table_S6_caption.txt").write_text(s6.CAPTION_S6 + "\n", encoding="utf-8")
assert len(sample_family) == 632 and len({r["sample_id"] for r in sample_meta}) == 158
assert not {"EXPSAMPLE-019", "EXPSAMPLE-020"} & {r["sample_id"] for r in sample_meta}
assert not any(r["life_stage"] == "larva" and r["tissue"] == "gall" for r in table_s6)
(OUT / "QC.json").write_text(json.dumps(qc, indent=2, sort_keys=True) + "\n", encoding="utf-8")
print(f"TABLE_S6_PASS rows={len(table_s6)} profiles={qc['profiles']} target_rows_added={qc['target_transcript_rows_added']}")
