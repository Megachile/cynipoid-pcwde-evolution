# Accepted source analyses for Tables S9–S13

The original analysis jobs are retained under the pinned source inventory. Code basenames identify historical producers; the external provenance map retains their original locations. Their exact result files are staged into `repository_inputs/tables/authority/supplement_new/`, with SHA256 pins in `config/deposit_sources.tsv`. `scripts/tables/build_tables_s9_s13.py` performs the final deterministic table projection; it does not rerun biological inference.

| Table | Historical producer code (basenames) | Input and method | Parameters |
|---|---|---|---|
| S9 | `analyze.py` | Accepted catalytic exceptions are mapped to the accepted transfer-context trees. Fitch changes are compared with state permutations among callable tips. | 10,000 permutations; per-family and per-site tests. |
| S10 | `null_core.py`, `summarize_null.py`, `analysis_summary.py` | Copy percentiles from fixed matched orthogroup nulls. The table distinguishes principal OG, all mapped OGs, and accepted paper roster. | Species/genome comparisons in source files; 95th percentile cutoff. |
| S11 | `build_heatmap_cells.py` | One row per species × family cell, preserving copy and chromosome occupancy and unresolved states. | 49 species × six families = 294 cells. |
| S12 | `build_amendment3_cases.py`, `compare_l1_splits.py`, `amendment4_t7_scf.py` | AU tests, site concordance and long-branch sensitivity from the accepted result files. | Fixed LG+G4, accepted alignments and IQ-TREE 3.0.1; 17 alternatives. |
| S13 | `source_inputs/di_inquiline_sort_limits/SORT_LIMITS.md`, `source_inputs/di_inquiline_sort_limits/COORDINATOR_ACCEPTANCE_NOTE.md` | Primary read-rule family counts from the accepted report; the acceptance note qualifies the lineage bound. | 34 audited tips; 14 own-bin, zero confirmed leaks, 20 unresolved. Lineage-qualified bound 27/34 from the note (report: 26/34), carried in `DI_LINEAGE_BOUND.tsv`. |

The table captions are factual drafts for Adam. The DI acceptance note changes the interpretation of the lineage-qualified bound from 26/34 to 27/34 without changing the read-rule table.
