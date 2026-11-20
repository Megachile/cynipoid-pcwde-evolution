# Supplementary table producers

## Current reader intermediates and pinned inputs

`config/logical_inputs.tsv` records the SHA256 for each staged input below; the
build verifies those hashes before a producer runs. This section describes
the processed layer between those inputs and the final PDFs.
The pinned input ledger itself has SHA256
`fb3a3e039d4cc45f514115ad5d7ce3d35e01c4f0eb9d5697d756b0c315652ae6`.

* **Table S3:** `scripts/tables/build_table_s2.py` reads the pinned origin
  sequence sets (GH5_2 `bf94586271cde876bfc707bc566fc7aa55770d300bd54194ae0a2906f2f5b2ed`,
  PL1_8 `1d4d60f3caa0cad4aac66573a75d6921ed8e06d0f3891780d8d37393d3a01256`,
  PL4_4 `77c42b56c356258d782365e6de2be74ab475fc984edaecde912beefdd1bfc06d`)
  and the B-lectin screen (`3a27e8885e0b7b705d82a0645ca24181171bcaf6612dd3008ab7008f76cfbfd0`).
  It maps each domain header to a source-qualified parent protein. Origin headers
  use the accession field, except ENA transcript headers, whose third field is
  the protein ID. B-lectin rotifer headers use the third-field accession with
  the domain suffix removed; distinct accessions within a species remain
  distinct proteins. Unknown schemas fail. The outputs are
  `table_s2/candidate/data/family_membership.tsv` and `domain_summary.tsv`;
  `render_tables.py` makes reader Table S3 from their protein counts, while
  `build_numbers.py` consumes domain and protein totals separately.
* **Table S6:** `build_cazy_rank_table.py` reads the accepted native CAZyme
  rank source (individual SHA256 pins in `config/logical_inputs.tsv`) and
  filters its reader intermediate to cynipoid rows. `render_tables.py` writes
  `tables/Table_S6.tsv` for the final PDF and supplement label check.
* **Table S7:** the chromosome-gain intermediate produced from the accepted
  crown comparison is read by `render_tables.py`; only DELTRAN rows enter
  `tables/Table_S7.tsv`. The renderer checks eight rows and one mode, then
  converts zero-gain fractions to percentages for the reader PDF.
* **Table S8:** `build_gc3_table.py` retains unrounded cohort statistics in
  `introns/gc3/Table_S16.tsv`; `render_tables.py` applies named headings and
  declared display rounding to `tables/Table_S8.tsv`. `check_intron_products.py`
  checks each displayed value against the technical source. The intron inputs
  and their individual SHA256 pins are in `config/logical_inputs.tsv` under
  `repository_inputs/analysis/introns/frozen/`.
* **Table S9:** `scripts/analysis/introns/04_build_table_s5.py` reads the
  pinned extraction and recurrent-position summaries from that same frozen
  intron set, plus `config/cazy_lineage_groups.tsv`
  (`4b5500803a229f20a5178bc880bef5d025fc16e1bb54d4a0f2299a49be6954b8`).
  Each position's member species maps to its taxonomic family; missing mappings
  fail. The technical `introns/Table_S5.tsv` retains species counts and the
  family list. `render_tables.py` projects that list and rounded percentages
  into reader `tables/Table_S9.tsv`; `check_intron_products.py` verifies the
  join and display.
* **Table S10:** `render_tables.py` combines the pinned selection and catalytic
  reader summaries (`repository_inputs/tables/authority/selection_current/` and
  `catalytic_current/`, with per-file SHA256 in the input ledger). It writes
  `tables/Table_S10.tsv` with explicit scorable denominators; the supplement
  typesetter consumes that TSV. No selection model is rerun.
* **Other changed intermediate producers:** `build_table_s6.py` and
  `build_table_s6_expression.py` project the pinned 160-profile expression
  roster into the 158-profile `table_s6/data/expression/` and S2a/S2b files;
  `build_tables_s3_s4.py` projects pinned member and synteny inputs into
  `tables_s3_s4/candidate/` for S4/S5; `ecami_report.py` produces the paired
  native rank summary consumed by `build_cazy_rank_table.py` for S6.
  `build_numbers.py` reads these same build intermediates into `numbers.json`,
  which manuscript and supplement token renderers consume.
  `produce_supplement_sections.py` fills only source-bound tokens into section
  text before `build_supplement.py` typesets it. `build_data_deposit.py` copies
  the validated scientific intermediates and source pins into the deposit;
  its manifest and closure checks consume the copy. All source paths and
  individual input hashes are in the pinned input ledger above.

`make manuscript` runs these producers and typesets their outputs (`render_tables.py`).
Inputs come only from the deposit's `repository_inputs/`, pinned in `config/logical_inputs.tsv`.

## Table S6: expression by biological context (`build_table_s6_expression.py`)

Inputs (`tables/authority/tables_s5_s6_current/`):
- `baseline/source_sample_roster.tsv`: the locked 160-profile RNA-seq roster (no DCAN).
- `baseline/sample_family_tpm.tsv`: per-sample family TPM (640 rows = 160 profiles × GH5_2, PL1_8, PL4_4, B-lectin). Each value is the sum of TPM over the family's eligible members in that sample, from annotation- and catalog-based quantification. This is the frozen expression layer.
- `baseline/gene_sample_tpm.tsv`: per-sample TPM for each target sequence (4,392 rows).
- `baseline/sample_metadata_and_mapping_qc.tsv` and `baseline/expression_sample_public_id_crosswalk.tsv`: life stage, tissue, mapping QC and public species names per sample.
- `target_transcript_expression.tsv`: TPM of the accepted Ceroptres B-lectin Trinity transcript `DI_ceroptres_TRINITY_DN781_c0_g1_i11.p1` in its six samples.

Steps (`build_table_s6.build_expression`):
1. **Public sample IDs.** Assign `EXPSAMPLE-001`… in sorted order of the roster's full sample IDs. Check 160 profiles and zero DCAN.
2. **Ceroptres target transcript.**
   - Take its six TPM values (RX02 1.9, RX03 45.0, RX04 1.2, WX02 0.7, WX04 29.4, WX05 66.0; asserted against `LOCKED_TARGET_TPM`).
   - Add each to that sample's B-lectin family TPM, which becomes quantification basis `helixer_plus_accepted_trinity_target`.
   - Append the six rows to the per-sequence table, giving 4,398 rows.
3. **Grouping.** Group the 640 sample-family values by public species name, life stage, tissue and family. Larval samples labelled `gall` are renamed `whole_body (in gall)`.
4. **Summary per group:**
   - sample count;
   - median;
   - IQR as the 25th–75th percentiles, interpolated linearly between order statistics (`percentile`);
   - range min–max.

   Values are shown to 0.1 TPM (`fmt_tpm`).

Outputs (`build/table_s6/`):
- `supplementary_tables/Table_S6.tsv` and its caption (`build_table_s6.CAPTION_S6`), typeset as Table S6;
- `data/expression/` with the adjusted sample-family table, the per-sequence table and the sample metadata;
- `QC.json`.

## Table S2a/S2b sample counts and TPM summary

The expression producer `scripts/tables/build_table_s6_expression.py` calls `sample_count_matrix` in `build_table_s6.py` after the accepted expression roster is filtered by `config/excluded_expression_samples.tsv`. The accepted baseline is `repository_inputs/tables/authority/tables_s5_s6_current/baseline/`: `source_sample_roster.tsv` (SHA256 84967a626f97e2dc3e7b46f7d8e94513a32d9895eeae4f9319e11786b899432b), `sample_metadata_and_mapping_qc.tsv` (5333fb97864b14c50b1bb9190f0b8533a0a38ee3db014bd9975bed68974a74a3), `expression_sample_public_id_crosswalk.tsv` (8e1e20e550ca9a2d0bae45780094f17db61c809e47660fe0326f54def15b12a6), `sample_family_tpm.tsv` (351be3a83351ec857c9bc412cf3e8e19250497daf14bf51bf9e7fb4b3a77baf4), and `gene_sample_tpm.tsv` (13e0a5e18548cf4fe2a124b1d4d4bc663f2b8b208aa29ade36704c98af8e6ac9). Source accessions come from `repository_inputs/tables/authority/table_s2/source_panel_accepted.tsv` (73fb17eea419a70812d160778fc3c3f67c67cdffa142c5e01d8f3ff96a165e90). The exclusion config SHA256 is ef75a5d0e63c8f32a3d4c912b9306108b32c557ff9d8f918e8e62c0600905fcc.

The baseline contains 160 expression profiles. The two existing excluded sample IDs, EXPSAMPLE-019 and EXPSAMPLE-020, are removed by the accepted filter; the panel uses all 158 retained profiles across 15 species and nine life-stage/tissue contexts. Each retained sample ID contributes once to its species/context cell, irrespective of gene or family rows. Zero cells are blank. Inquiline profile IDs can share originating libraries, so 158 is a profile count, not a count of independent biological libraries. Species display names follow the accepted crosswalk; source labels are resolved against the pinned source panel and an unresolved label fails the producer.

The layer-2 intermediate is `table_s6/data/expression/sample_count_matrix.tsv`; its same-build copy is `tables/Table_S2a.tsv`. `scripts/tables/render_tables.py` renders Table S2a from that intermediate. The Table S2b TPM summary is `table_s6/supplementary_tables/Table_S6.tsv` internally and is checked cell by cell against the panel's sample counts. The numbers producer consumes the same intermediate. Duplicate gene rows cannot increase the count; a mismatch with S2b fails.

## Table S7: selection and catalytic sites (`build_table_s7.py`)

Reads `tables/authority/selection_current/` (positive selection, frozen to the 41-genome scope) and `tables/authority/catalytic_current/`. It writes the two reader tables that make up Table S7: `selection_reader_table.tsv` (purifying, diversifying and terminal-branch percentages) and `catalytic_reader_table.tsv` (de novo alignment, catalytic sites), plus supporting files.

## Not built yet
Table S5 (intron table; the corrected-model intron rerun will build it) is still typeset from its saved copy in `tables/accepted/`. The build's closure check lists it as a known gap.

## Table S1 (`build_table_s1_rnaseq.py`, typeset by `build_table_s1.py`)

`build_table_s1_inventory.py` reconstructs 99 resource rows from the pinned
`tables/authority/table_s2/source_panel_accepted.tsv`. It maps source type,
preparation, accession, status and order to the established S1 presentation.
The Periclistus pair is ordered by display name, and four inquiline display
names follow the accepted paper labels. Eleven BUSCO fields missing from the
source panel are supplied by `tables/authority/table_s1_catalog/catalog_busco.tsv`,
a frozen catalog snapshot export. Ten protein-mode scores have matching October
16 reports pinned beside the export. These are corroborating later measurements;
the pre-October S1 fill chain is unrecorded. The Mayetiola 95.7 score is the
catalog's insecta_odb10 assembly-mode result for assembly 114, not its
protein-mode 92.5. Druon annotation 205 supplies 76.2 and has a matching
report; historical annotation 29 also has 76.2, and the old annotation choice
is unrecorded. The export script requires a gateway-materialized host-local
snapshot and records its SHA-256. The historical saved S1 is used only for
external cell-and-order regression, never as a build input.

The RNA-seq join reads the built inventory, the accepted caption and the
built Table S6 (`build/table_s6/supplementary_tables/Table_S6.tsv`). Every Table S1 row and cell is carried unchanged except
`resource_class`: a Genome row whose species has a library in the built Table S6 becomes "Genome + RNA-seq", and every
Transcriptome-derived proteome row becomes "Transcriptome (RNA-seq)". Species join by exact `species_display`; a Table S6
species without a Table S1 row fails the build (`build/table_s1/join_report.tsv`). One caption sentence defines the two
labels. `scripts/release/check_table_s1.py` re-derives both label sets and checks every other cell is unchanged.

## Tables S9–S13 (`build_tables_s9_s13.py`)

The stage step pins the accepted summaries in `repository_inputs/tables/authority/supplement_new/`; `config/deposit_sources.tsv` records every upstream SHA256 and `config/logical_inputs.tsv` pins the staged bytes. The builder is a deterministic projection of the 2026-10-28 uniform preview (`prepare_new_tables.py`), with each TSV checked against its preview SHA256 before rendering. Table S9 uses the accepted 10,000-permutation catalytic result. Table S10 keeps principal orthogroup, all mapped orthogroups, and paper roster as distinct units. Table S11 retains all 294 matrix cells. Table S12 combines the accepted AU, sCF and long-branch files. Table S13 parses the accepted read-rule table (14 own-bin, 0 confirmed leaks, 20 unresolved of 34). The lineage-qualified bound is a separate value: the builder parses the report's `(34 − 8)/34 = 26/34` and the acceptance note's `7/34` strongly attributed, `27/34`, checks both arithmetic and the shared denominator, and writes `build/tables_s9_s13/DI_LINEAGE_BOUND.tsv`. The caption and `numbers.json` (`table_s13_lineage_bound`) both come from that file; `build_numbers.py` fails if the caption source or the typeset caption disagree. Captions are factual drafts for Adam. All five tables use the existing `render_tables.py` style.
The upstream code and method pointers for these five producers are in `upstream/supplement_table_sources/README.md`.

### BIP adult 454 exclusion

The shared Table S6 expression producer reads `config/excluded_expression_samples.tsv` after checking the unchanged 160-profile accepted input. It excludes `BIP_adult_OVX` and `BIP_adult_VGX` because those adult 454 libraries are outside the current larval-only BIP reference. Public EXPSAMPLE IDs are assigned from the full accepted roster before filtering, so all retained IDs and TPM values stay stable. Both Table S6 entry points write the same retained crosswalk (`table_s6/data/expression/retained_sample_crosswalk.tsv`) with full, original and public sample IDs and species identity. The reader deposit requires this crosswalk and joins on species display plus original sample ID; original IDs can repeat across species. The frozen inputs themselves are unchanged.
