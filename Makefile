PYTHON ?= python3
RSCRIPT ?= Rscript
GS ?= gs
BUILD ?= $(CURDIR)/build
DATA_ROOT ?= $(if $(PCWDE_DATA_ROOT),$(PCWDE_DATA_ROOT),$(CURDIR)/figshare_data)
INTRON_SOURCES := config/intron_snapshot_sources.tsv
FROZEN_INTRONS := $(DATA_ROOT)/repository_inputs/analysis/introns/frozen
BACTERIAL_GRID := $(FROZEN_INTRONS)/bacterial_grid.tsv
BLECTIN_BACTERIAL_CDS := $(FROZEN_INTRONS)/blectin_bacterial_cds.fna
export PYTHONDONTWRITEBYTECODE=1

.PHONY: verify manuscript paper intron-preflight intron-gc3 intron-fixture intron-input-check deposit-input-check data-deposit stage-check clean

paper: manuscript

# Retained chromosome-gain analysis against native eCAMI groups.
CHROM_GAIN_INPUTS = $(DATA_ROOT)/repository_inputs
CHROM_GAIN_PANEL = $(CHROM_GAIN_INPUTS)/chromosome_gain/PANEL_GENOMES.tsv
CHROM_GAIN_TREE = $(CHROM_GAIN_INPUTS)/figures/chromosome_counts/species_tree.nwk
CHROM_GAIN_OUT = $(BUILD)/chromosome_gain
CHROM_GAIN_FIGURE_CELLS = $(BUILD)/figures/tree_counts/full_combined_cells.tsv

.PHONY: ecami-chromosome-gain
ecami-chromosome-gain:
	test -f $(BUILD)/cazy_counts/component_gene_units.tsv -a -f $(BUILD)/members/copy_placements.tsv -a -f $(CHROM_GAIN_FIGURE_CELLS)
	mkdir -p $(CHROM_GAIN_OUT)
	$(PYTHON) -B scripts/analysis/chromosome_gain/verify_cells.py --focal $(BUILD)/members/copy_placements.tsv --figure-cells $(CHROM_GAIN_FIGURE_CELLS) --out $(CHROM_GAIN_OUT)/figure1_cell_qc.json
	$(PYTHON) -B scripts/analysis/chromosome_gain/build_states.py --data-root $(DATA_ROOT) --out $(CHROM_GAIN_OUT)/native_protein_states.tsv --protein-genes $(BUILD)/cazy_counts/protein_genes.tsv --scan-manifest config/cazy_scans.tsv --panel $(CHROM_GAIN_PANEL)
	$(PYTHON) -B scripts/analysis/chromosome_gain/ecami_score.py --states $(CHROM_GAIN_OUT)/native_protein_states.tsv --components $(BUILD)/cazy_counts/component_gene_units.tsv --group-mapping $(BUILD)/cazy_counts/group_mapping.tsv --protein-genes $(BUILD)/cazy_counts/protein_genes.tsv --native-counts $(BUILD)/cazy_ranks/native_component_counts.tsv --crosswalk $(BUILD)/cazy_crosswalk/member_crosswalk.tsv --primary-ranks $(BUILD)/cazy_ranks/primary_ecami_ranks.tsv --focal $(BUILD)/members/copy_placements.tsv --panel $(CHROM_GAIN_PANEL) --tree $(CHROM_GAIN_TREE) --figure-code scripts/figures/chromosome_counts.py --expected-focal config/ecami_focal_parity.tsv --out $(CHROM_GAIN_OUT)
	$(PYTHON) -B scripts/analysis/chromosome_gain/handcheck.py --ecami --scores $(CHROM_GAIN_OUT)/ecami_scores.tsv --tips $(CHROM_GAIN_OUT)/ecami_tip_states.tsv --changes $(CHROM_GAIN_OUT)/ecami_node_changes.tsv --tree $(CHROM_GAIN_TREE) --out $(CHROM_GAIN_OUT)/hand_checks.tsv
	$(PYTHON) -B scripts/analysis/chromosome_gain/ecami_report.py --scores $(CHROM_GAIN_OUT)/ecami_scores.tsv --eligibility $(CHROM_GAIN_OUT)/background_eligibility.tsv --scan-eligibility $(CHROM_GAIN_OUT)/focal_scan_eligibility.tsv --component-sizes $(CHROM_GAIN_OUT)/focal_component_sizes.tsv --handchecks $(CHROM_GAIN_OUT)/hand_checks.tsv --table-numbers config/table_numbers.tsv --figure-numbers config/figure_numbers.tsv --out $(CHROM_GAIN_OUT)/report

# Diagnostic only: require an explicit source root and frozen deposit snapshot.
# A successful invocation can report unresolved identities; inspect its TSVs.
.PHONY: cazy-diagnostic
cazy-diagnostic:
	@test -f "$(DATA_ROOT)/repository_inputs/members/member_table.tsv" || { echo 'frozen member_table.tsv is required' >&2; exit 2; }
	$(PYTHON) -B scripts/tables/diagnose_cazy_inputs.py --manifest config/cazy_scans.tsv --paper-sources config/cazy_paper_proteomes.tsv --member-table $(DATA_ROOT)/repository_inputs/members/member_table.tsv --data-root $(DATA_ROOT) --out $(BUILD)/cazy_diagnostic

.PHONY: cazy-science
cazy-science:
	@test -f config/cazy_processed_proteomes.tsv || { echo 'processed scan input pin is required' >&2; exit 2; }
	$(PYTHON) -B scripts/tables/build_cazy_scan_counts.py --manifest config/cazy_scans.tsv --data-root $(DATA_ROOT) --fig1-counts $(DATA_ROOT)/repository_inputs/figures/presence_counts --out $(BUILD)/cazy_counts
	$(PYTHON) -B scripts/tables/reconcile_cazy_members.py --manifest config/cazy_scans.tsv --processed-map config/cazy_processed_proteomes.tsv --member-table $(DATA_ROOT)/repository_inputs/members/member_table.tsv --sequences-dir $(DATA_ROOT)/sequences --data-root $(DATA_ROOT) --out $(BUILD)/cazy_crosswalk
	$(PYTHON) -B scripts/tables/build_cazy_ranks.py --scan-counts $(BUILD)/cazy_counts --member-table $(DATA_ROOT)/repository_inputs/members/member_table.tsv --member-crosswalk $(BUILD)/cazy_crosswalk/member_crosswalk.tsv --manifest config/cazy_scans.tsv --out $(BUILD)/cazy_ranks
	$(PYTHON) -B scripts/tables/build_cazy_deposit.py --manifest config/cazy_scans.tsv --unit-counts $(BUILD)/cazy_counts/unit_counts.tsv --paired-ranks $(BUILD)/cazy_ranks/paired_ranks.tsv --out $(BUILD)/cazy_deposit/all_unit_counts.tsv
	$(PYTHON) -B scripts/tables/build_cazy_lineage_summary.py --lineages config/cazy_lineage_groups.tsv --paired-ranks $(BUILD)/cazy_ranks/paired_ranks.tsv --native-counts $(BUILD)/cazy_ranks/native_counts.tsv --component-ranks $(BUILD)/cazy_ranks/component_ranks.tsv --primary-ranks $(BUILD)/cazy_ranks/primary_ecami_ranks.tsv --native-components $(BUILD)/cazy_ranks/native_component_counts.tsv --group-mapping $(BUILD)/cazy_counts/group_mapping.tsv --primary-precision $(BUILD)/cazy_ranks/primary_ecami_precision.tsv --out $(BUILD)/cazy_summary

verify:
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/release/portable_sources.py --repo $(CURDIR)
	$(PYTHON) -B scripts/release/check_readable_names.py --repo $(CURDIR) --data-root $(DATA_ROOT)
	$(PYTHON) -B scripts/release/test_readable_names.py
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/deposit_inputs/test_public_scope.py
	PCWDE_REPOSITORY_INPUTS=$(DATA_ROOT)/repository_inputs PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/tables/test_sample_count_matrix.py
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/release/test_supplement_page_numbers.py
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/release/test_supplement_cleanup.py
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/release/test_main_text_pointers.py
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/release/test_pdf_glyphs.py
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/release/check_intron_retirement.py
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/release/verify_repository.py --repo $(CURDIR) --data-root $(DATA_ROOT)

data-deposit: verify data-deposit-assemble

.PHONY: data-deposit-assemble
data-deposit-assemble:
	test -s $(BUILD)/table_s6/data/expression/retained_sample_crosswalk.tsv
	PCWDE_REPO_ROOT=$(CURDIR) PCWDE_DATA_ROOT=$(DATA_ROOT) PCWDE_DEPOSIT_OUT=$(BUILD)/data_deposit PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/release/build_data_deposit.py
	$(PYTHON) -B scripts/release/check_deposit_size.py --payload $(BUILD)/data_deposit/pcwde_supporting_data --baseline config/deposit_size_baseline.tsv --receipt $(BUILD)/data_deposit/DEPOSIT_SIZE.json
	PCWDE_REPO_ROOT=$(CURDIR) PCWDE_DATA_ROOT=$(DATA_ROOT) PCWDE_DEPOSIT_OUT=$(BUILD)/data_deposit PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/release/validate_data_deposit.py

deposit-input-check:
	PCWDE_REPO_ROOT=$(CURDIR) PCWDE_DATA_ROOT=$(DATA_ROOT) PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/release/build_data_deposit.py --validate-inputs-only

intron-input-check:
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/analysis/introns/preflight.py --data-root $(DATA_ROOT) --sources $(INTRON_SOURCES) --members $(DATA_ROOT)/repository_inputs/members/member_table.tsv --blectin-models config/blectin_models.tsv --frozen $(FROZEN_INTRONS) --out $(BUILD)/introns/gc3 --validate-inputs-only

intron-fixture:
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/analysis/introns/fixture.py --work $(BUILD)/intron_fixture

intron-preflight: verify intron-fixture
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/analysis/introns/preflight.py --data-root $(DATA_ROOT) --sources $(INTRON_SOURCES) --members $(DATA_ROOT)/repository_inputs/members/member_table.tsv --blectin-models config/blectin_models.tsv --frozen $(FROZEN_INTRONS) --out $(BUILD)/introns/gc3

intron-gc3: intron-preflight
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/analysis/introns/stage_frozen.py --frozen $(FROZEN_INTRONS) --kind member_cds --out $(BUILD)/introns/member_cds
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/analysis/introns/build_member_observations.py --accounting $(BUILD)/introns/member_cds/cohort_accounting.tsv --segments $(BUILD)/introns/member_cds/member_cds_segments.tsv --member-table $(DATA_ROOT)/repository_inputs/members/member_table.tsv --blectin-accounting $(BUILD)/introns/member_cds/blectin_accounting.tsv --blectin-segments $(BUILD)/introns/member_cds/blectin_cds_segments.tsv --sources $(INTRON_SOURCES) --project-root $(DATA_ROOT) --index-dir $(BUILD)/introns/indexes --focal-cds $(FROZEN_INTRONS)/focal_cds.fna --focal-junctions $(FROZEN_INTRONS)/focal_junctions.tsv --out $(BUILD)/introns/observations
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/analysis/introns/stage_frozen.py --frozen $(FROZEN_INTRONS) --kind annotation_proteins --out $(BUILD)/introns/annotation_proteins
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/analysis/introns/stage_frozen.py --frozen $(FROZEN_INTRONS) --kind host_background --out $(BUILD)/introns/host_background
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/analysis/introns/02_map_introns.py --input-dir $(BUILD)/introns/observations --sequence-dir $(BUILD)/introns/annotation_proteins --member-table $(DATA_ROOT)/repository_inputs/members/member_table.tsv --output-dir $(BUILD)/introns/mapped
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/analysis/introns/03_summarize_introns.py --input-dir $(BUILD)/introns/mapped --sequence-dir $(BUILD)/introns/annotation_proteins --output-dir $(BUILD)/introns/summary
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/analysis/introns/04_build_table_s5.py --extraction-dir $(BUILD)/introns/observations --summary-dir $(BUILD)/introns/summary --b-lectin-extraction-dir $(BUILD)/introns/observations --output-dir $(BUILD)/introns/table_s5
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/analysis/introns/summarize_gc3.py --host-dir $(BUILD)/introns/host_background --member-table $(DATA_ROOT)/repository_inputs/members/member_table.tsv --blectin-models config/blectin_models.tsv --focal-composition $(BUILD)/introns/observations/focal_composition.tsv --blectin-composition $(BUILD)/introns/observations/blectin_focal_composition.tsv --bacterial-grid $(BACTERIAL_GRID) --blectin-bacterial-cds $(BLECTIN_BACTERIAL_CDS) --out $(BUILD)/introns/gc3
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/analysis/introns/check_reproduction.py --accounting $(BUILD)/introns/member_cds/cohort_accounting.tsv --observations $(BUILD)/introns/observations/per_intron_observations.tsv --table-s5 $(BUILD)/introns/table_s5/Table_S5.tsv --host-qc $(BUILD)/introns/host_background/QC.tsv --gc3-stats $(BUILD)/introns/gc3/pcwde_gc3_stats.json --focal-composition $(BUILD)/introns/observations/focal_composition.tsv

manuscript: verify
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/release/check_environment.py --rscript $(RSCRIPT)
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/release/test_check_closure.py
	$(PYTHON) -B scripts/release/figure_build_gate.py --repo $(CURDIR) --data-root $(DATA_ROOT) --build $(BUILD)
	mkdir -p $(BUILD)
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/release/check_closure.py --reset --data-root $(DATA_ROOT) --build $(BUILD)
	$(MAKE) s4-architecture-replay PYTHON=$(PYTHON) BUILD=$(BUILD) DATA_ROOT=$(DATA_ROOT)
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/members/build_copy_table.py --data-root $(DATA_ROOT) --out $(BUILD)/members/copy_placements.tsv
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/members/build_locus_tips.py --data-root $(DATA_ROOT) --copy-table $(BUILD)/members/copy_placements.tsv --out $(BUILD)/members/locus_tips.tsv
	$(RSCRIPT) scripts/trees/prune_origin_display_trees.R --complete-dir $(DATA_ROOT)/repository_inputs/deposit/base/trees --tips config/origin_display_tips.tsv --logical-inputs config/logical_inputs.tsv --out $(BUILD)/trees
	$(MAKE) cazy-science PYTHON=$(PYTHON) BUILD=$(BUILD) DATA_ROOT=$(DATA_ROOT)
	$(MAKE) intron-gc3 BUILD=$(BUILD) DATA_ROOT=$(DATA_ROOT) PYTHON=$(PYTHON)
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/analysis/introns/assemble_figure_inputs.py --observations $(BUILD)/introns/observations --summary $(BUILD)/introns/summary --background $(BUILD)/introns/host_background --out $(BUILD)/introns/figure_inputs
	RENDER_COPY_TABLE=$(BUILD)/members/copy_placements.tsv RENDER_LOCUS_TIPS=$(BUILD)/members/locus_tips.tsv RENDER_INTRONS_INPUT_DIR=$(BUILD)/introns/figure_inputs RENDER_REPO=$(CURDIR) RENDER_REPO_BUILD=$(BUILD) RENDER_DATA=$(DATA_ROOT)/repository_inputs/figures RENDER_OUT=$(BUILD)/figures RENDER_PYTHON=$(PYTHON) RENDER_RSCRIPT=$(RSCRIPT) /bin/bash scripts/release/render_all.bash
	$(MAKE) ecami-chromosome-gain PYTHON=$(PYTHON) BUILD=$(BUILD) DATA_ROOT=$(DATA_ROOT)
	$(PYTHON) -B scripts/release/figure_build_gate.py --repo $(CURDIR) --data-root $(DATA_ROOT) --build $(BUILD) --receipt
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/release/figure_numbers.py --build $(BUILD) --data-root $(DATA_ROOT)
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/release/blectin_domain_counts.py --tree $(BUILD)/figures/b_lectin_domain_tree/primary/trees/Figure_S4_B_lectin_domains.rooted.treefile --manifest $(BUILD)/figures/b_lectin_domain_tree/primary/sidecars/Figure_S4_domain_manifest.tsv --out $(BUILD)/blectin_domain_tree_counts.json
	PCWDE_COPY_TABLE=$(BUILD)/members/copy_placements.tsv PCWDE_REPOSITORY_INPUTS=$(DATA_ROOT)/repository_inputs PCWDE_BUILD_ROOT=$(BUILD) PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/tables/build_table_s2.py
	PCWDE_COPY_TABLE=$(BUILD)/members/copy_placements.tsv PCWDE_REPOSITORY_INPUTS=$(DATA_ROOT)/repository_inputs PCWDE_BUILD_ROOT=$(BUILD) PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/tables/build_tables_s3_s4.py
	PCWDE_REPOSITORY_INPUTS=$(DATA_ROOT)/repository_inputs PCWDE_BUILD_ROOT=$(BUILD) PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/tables/build_table_s6_expression.py
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/tables/build_table_s1_inventory.py --source-panel $(DATA_ROOT)/repository_inputs/tables/authority/table_s2/source_panel_accepted.tsv --catalog-pin $(DATA_ROOT)/repository_inputs/tables/authority/table_s1_catalog/catalog_busco.tsv --out $(BUILD)/table_s1_inventory
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/tables/build_table_s1_rnaseq.py --table-s1-dir $(BUILD)/table_s1_inventory --caption $(DATA_ROOT)/repository_inputs/tables/accepted/Table_S1_caption.txt --table-s6 $(BUILD)/table_s6/supplementary_tables/Table_S6.tsv --out $(BUILD)/table_s1
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/analysis/catalytic/subset_current_origin.py --data-root $(DATA_ROOT) --out $(BUILD)/catalytic_current
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/tables/build_table_s7.py --selection $(DATA_ROOT)/repository_inputs/tables/authority/selection_current --catalytic $(BUILD)/catalytic_current --table-s1 $(BUILD)/table_s1_inventory/Table_S1.tsv --out $(BUILD)/table_s7
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/tables/build_tables_s9_s13.py --data-root $(DATA_ROOT) --out $(BUILD)/tables_s9_s13 --analysis $(DATA_ROOT)/repository_inputs/tables/authority/supplement_new/catalytic_clustering/ANALYSIS.json --catalytic-summary $(BUILD)/catalytic_current/catalytic_site_conservation_summary.tsv
	$(PYTHON) -B scripts/tables/build_cazy_rank_table.py --paired-ranks $(BUILD)/cazy_ranks/paired_ranks.tsv --primary-ranks $(BUILD)/cazy_ranks/primary_ecami_ranks.tsv --numbers config/table_numbers.tsv --lineage-groups config/cazy_lineage_groups.tsv --lineage-precision $(BUILD)/cazy_summary/lineage_precision.tsv --out $(BUILD)/cazy_rank_table
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/members/sequence_set_roster.py --data-root $(DATA_ROOT) --build $(BUILD) --out $(BUILD)/sequence_sets
	$(MAKE) intron-paper-products BUILD=$(BUILD) DATA_ROOT=$(DATA_ROOT) PYTHON=$(PYTHON)
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/release/check_origin_key_nodes.py --build $(BUILD)
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/release/check_origin_key_geometry.py --figures $(BUILD)/figures/origin_trees --out $(BUILD)/figures/origin_trees/evidence/key_geometry.json

	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/release/check_table_s1.py --data-root $(DATA_ROOT) --build $(BUILD)
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/release/check_copy_counts.py --data-root $(DATA_ROOT) --copy-table $(BUILD)/members/copy_placements.tsv --locus-tips $(BUILD)/members/locus_tips.tsv --build $(BUILD)
	PCWDE_MANUSCRIPT_FIGURES=$(BUILD)/figures PCWDE_MANUSCRIPT_OUT=$(BUILD)/manuscript PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/release/build_manuscript.py
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/release/test_supplement_template.py --repo $(CURDIR) --numbers $(BUILD)/numbers.json
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/release/test_blectin_domain_counts.py --tree $(BUILD)/figures/b_lectin_domain_tree/primary/trees/Figure_S4_B_lectin_domains.rooted.treefile --manifest $(BUILD)/figures/b_lectin_domain_tree/primary/sidecars/Figure_S4_domain_manifest.tsv --counts $(BUILD)/blectin_domain_tree_counts.json --table-s3-caption $(BUILD)/table_s2/candidate/supplementary_tables/Table_S2_caption.txt --figure-numbers config/figure_numbers.tsv
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/release/produce_supplement_sections.py --repo $(CURDIR) --numbers $(BUILD)/numbers.json --report $(BUILD)/supplement_numbers_gate.txt --out $(BUILD)/supplement_section_text.json
	PCWDE_BUILD_ROOT=$(BUILD) PCWDE_DATA_ROOT=$(DATA_ROOT) PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/release/build_supplement.py
	$(PYTHON) -B scripts/release/check_supplement_labels.py --repo $(CURDIR) --build $(BUILD)
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/release/render_supplement_pages.py --gs $(GS) --pdf $(BUILD)/supplement/supplement.pdf --out $(BUILD)/supplement_pages
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/release/check_closure.py --inputs --data-root $(DATA_ROOT) --build $(BUILD)
	$(MAKE) data-deposit PYTHON=$(PYTHON) BUILD=$(BUILD) DATA_ROOT=$(DATA_ROOT)
	$(MAKE) intron-product-check BUILD=$(BUILD) PYTHON=$(PYTHON)
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/release/check_closure.py --deposit --data-root $(DATA_ROOT) --build $(BUILD)

stage-check:
	test -n "$(BASE)" -a -n "$(WORKSPACE)" || { echo 'stage-check needs BASE and WORKSPACE' >&2; exit 2; }
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/deposit_inputs/stage.py --base $(BASE) --workspace $(WORKSPACE) --out $(BUILD)/staged_deposit/pcwde_supporting_data
	# Compare the restaged root with DATA_ROOT in the same form (raw source pins or
	# projected manifest), then always prove the restaged root projects cleanly.
	if test -e $(DATA_ROOT)/PROJECTED_INPUTS.tsv; then \
	  PCWDE_REPO_ROOT=$(CURDIR) PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/release/build_data_deposit.py --project-inputs-only $(BUILD)/staged_deposit/pcwde_supporting_data/repository_inputs && \
	  PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/deposit_inputs/stage.py --compare $(BUILD)/staged_deposit/pcwde_supporting_data $(DATA_ROOT); \
	else \
	  PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/deposit_inputs/stage.py --compare $(BUILD)/staged_deposit/pcwde_supporting_data $(DATA_ROOT) && \
	  PCWDE_REPO_ROOT=$(CURDIR) PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/release/build_data_deposit.py --project-inputs-only $(BUILD)/staged_deposit/pcwde_supporting_data/repository_inputs; \
	fi

clean:
	@printf '%s\n' 'Generated files are beneath build/; remove that directory explicitly if desired.'

.PHONY: previews
PREVIEW_GS ?= gs
previews:
	@test -d "$(BUILD)" || { echo "missing build tree: $(BUILD)" >&2; exit 2; }
	@command -v "$(PREVIEW_GS)" >/dev/null || { echo "Ghostscript unavailable in build environment: $(PREVIEW_GS)" >&2; exit 2; }
	@test -f "$(BUILD)/supplement/supplement.pdf" || { echo "missing supplement PDF" >&2; exit 2; }
	@mkdir -p "$(BUILD)/previews"; \
	  "$(PREVIEW_GS)" -q -dNOPAUSE -dBATCH -dSAFER -sDEVICE=png16m -r100 -sOutputFile="$(BUILD)/previews/Supplement_%03d.png" "$(BUILD)/supplement/supplement.pdf"

# Targeted visual fixture: same expression producer and S2a renderer as the full paper build.
.PHONY: table-s2a-fixture
table-s2a-fixture:
	@test -d "$(DATA_ROOT)/repository_inputs" -a -d "$(DATA_ROOT)/repository_inputs/tables/accepted" || { echo 'missing pinned repository inputs or table base' >&2; exit 2; }
	@test -s "$(DATA_ROOT)/build_inputs/euceroptres_genome/euceroptres_genome.pdf" || { echo 'missing pinned S5 asset' >&2; exit 2; }
	@command -v xelatex >/dev/null && command -v kpsewhich >/dev/null || { echo 'TeX tools unavailable' >&2; exit 2; }
	@test ! -e "$(BUILD)/table_s6" -a ! -e "$(BUILD)/tables" || { echo 'fixture output exists' >&2; exit 2; }
	PCWDE_REPOSITORY_INPUTS=$(DATA_ROOT)/repository_inputs PCWDE_BUILD_ROOT=$(BUILD) PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/tables/build_table_s6_expression.py
	PCWDE_REPOSITORY_INPUTS=$(DATA_ROOT)/repository_inputs PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/tables/test_sample_count_matrix.py
	@if [ -n "$(S2A_LOCAL_ROOT)" ]; then \
	  test -d "$(S2A_LOCAL_ROOT)" && test ! -e "$(S2A_LOCAL_ROOT)/s2a"; \
	  mkdir -p "$(S2A_LOCAL_ROOT)/s2a/scripts/tables" "$(S2A_LOCAL_ROOT)/s2a/config" "$(S2A_LOCAL_ROOT)/s2a/base" "$(S2A_LOCAL_ROOT)/s2a/table_s6/data/expression" "$(S2A_LOCAL_ROOT)/s2a/table_s6/supplementary_tables" "$(S2A_LOCAL_ROOT)/s2a/fonts"; \
	  cp scripts/tables/render_tables.py "$(S2A_LOCAL_ROOT)/s2a/scripts/tables/"; \
	  cp config/table_numbers.tsv config/figure_numbers.tsv config/cazy_lineage_groups.tsv "$(S2A_LOCAL_ROOT)/s2a/config/"; \
	  cp "$(BUILD)/table_s6/data/expression/sample_count_matrix.tsv" "$(S2A_LOCAL_ROOT)/s2a/table_s6/data/expression/"; \
	  cp "$(BUILD)/table_s6/supplementary_tables/Table_S6a_caption.txt" "$(S2A_LOCAL_ROOT)/s2a/table_s6/supplementary_tables/"; \
	  font_dir=$$(dirname "$$(kpsewhich texgyretermes-regular.otf)"); \
	  for style in regular bold italic bolditalic; do cp "$$font_dir/texgyretermes-$$style.otf" "$(S2A_LOCAL_ROOT)/s2a/fonts/"; done; \
	  (cd "$(S2A_LOCAL_ROOT)/s2a" && sha256sum scripts/tables/render_tables.py config/*.tsv table_s6/data/expression/sample_count_matrix.tsv table_s6/supplementary_tables/Table_S6a_caption.txt fonts/*.otf > staged_sha256.txt && sha256sum -c staged_sha256.txt); \
	  PCWDE_TEX_GYRE_DIR="$(S2A_LOCAL_ROOT)/s2a/fonts" TEXMFVAR="$(S2A_LOCAL_ROOT)/s2a/texmf-var" TEXMFCONFIG="$(S2A_LOCAL_ROOT)/s2a/texmf-config" TEXMFCACHE="$(S2A_LOCAL_ROOT)/s2a/texmf-cache" PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B "$(S2A_LOCAL_ROOT)/s2a/scripts/tables/render_tables.py" --s2a-only --base "$(S2A_LOCAL_ROOT)/s2a/base" --built "$(S2A_LOCAL_ROOT)/s2a/table_s6/supplementary_tables" --out "$(S2A_LOCAL_ROOT)/s2a/tables"; \
	  test -s "$(S2A_LOCAL_ROOT)/s2a/tables/Table_S2a.pdf"; \
	  mkdir -p "$(BUILD)/tables"; cp "$(S2A_LOCAL_ROOT)/s2a/tables/"Table_S2a.* "$(BUILD)/tables/"; cp "$(S2A_LOCAL_ROOT)/s2a/staged_sha256.txt" "$(BUILD)/tables/S2a_staged_sha256.txt"; \
	else \
	  PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/tables/render_tables.py --s2a-only --base $(DATA_ROOT)/repository_inputs/tables/authority --built $(BUILD)/table_s6/supplementary_tables --out $(BUILD)/tables; \
	fi
	@test -s "$(BUILD)/tables/Table_S2a.pdf" -a -s "$(BUILD)/tables/Table_S2a.tsv" -a -s "$(BUILD)/table_s6/supplementary_tables/Table_S6a_caption.txt"

.PHONY: intron-paper-products intron-product-check
intron-paper-products:
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/release/build_gc3_table.py --pcwde $(BUILD)/introns/gc3/pcwde_gc3_stats.json --blectin $(BUILD)/introns/gc3/blectin_gc3_stats.json --out $(BUILD)/introns/gc3
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/tables/render_tables.py --base $(DATA_ROOT)/repository_inputs/tables/authority --built $(BUILD)/table_s2/candidate/supplementary_tables --built $(BUILD)/tables_s3_s4/candidate/supplementary_tables --built $(BUILD)/table_s6/supplementary_tables --built $(BUILD)/table_s1 --s7 $(BUILD)/table_s7 --built $(BUILD)/tables_s9_s13 --built $(BUILD)/cazy_rank_table --built $(CHROM_GAIN_OUT)/report --built $(BUILD)/introns/table_s5 --built $(BUILD)/introns/gc3 --out $(BUILD)/tables
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/release/build_numbers.py --data-root $(DATA_ROOT) --copy-table $(BUILD)/members/copy_placements.tsv --table-s2 $(BUILD)/table_s2/candidate/supplementary_tables/Table_S2.tsv --domain-summary $(BUILD)/table_s2/candidate/data/domain_summary.tsv --origin-key-nodes $(BUILD)/trees/origin_key_nodes.tsv --new-tables $(BUILD)/tables_s9_s13 --rendered-tables $(BUILD)/tables --synteny-table $(BUILD)/tables_s3_s4/candidate/supplementary_tables/Table_S4.tsv --synteny-summary $(BUILD)/tables_s3_s4/synteny_summary.tsv --chromosome-gain-summary $(CHROM_GAIN_OUT)/report/focal_summary.tsv --chromosome-gain-table $(CHROM_GAIN_OUT)/report/Table_S15.tsv --cazy-rank-table $(BUILD)/cazy_rank_table/Table_S14.tsv --intron-summary $(BUILD)/introns/summary/recurrent_position_summary.tsv --intron-observations $(BUILD)/introns/observations/per_intron_observations.tsv --intron-composition $(BUILD)/introns/observations/focal_composition.tsv --intron-comparator $(BUILD)/introns/host_background/whole_genome_comparator.json --table-s5 $(BUILD)/introns/table_s5/Table_S5.tsv --gc3-stats $(BUILD)/introns/gc3/pcwde_gc3_stats.json --blectin-gc3-stats $(BUILD)/introns/gc3/blectin_gc3_stats.json --blectin-gc3-comparison $(BUILD)/introns/gc3/blectin_gc3_comparison.tsv --source-panel $(DATA_ROOT)/repository_inputs/tables/authority/table_s2/source_panel_accepted.tsv --expression-table $(BUILD)/table_s6/supplementary_tables/Table_S6.tsv --expression-sample-count-matrix $(BUILD)/table_s6/data/expression/sample_count_matrix.tsv --family-membership $(BUILD)/table_s2/candidate/data/family_membership.tsv --sequence-sets $(BUILD)/sequence_sets/pcwde_sequence_set_summary.tsv --b-lectin-manifest $(BUILD)/figures/b_lectin_domain_tree/primary/sidecars/Figure_S4_domain_manifest.tsv --b-lectin-tree-counts $(BUILD)/blectin_domain_tree_counts.json --catalytic-domains $(DATA_ROOT)/repository_inputs/deposit/base/sequences/catalytic_domains --out $(BUILD)/numbers.json

intron-product-check:
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/release/check_intron_products.py $(BUILD)

# Bounded contracts use the same producers as paper, without a saved product fallback.
.PHONY: s6-render source-panel-public
s6-render:
	RENDER_RSCRIPT=$(RSCRIPT) $(PYTHON) -B scripts/figures/b_lectin_domain_tree/build_b_lectin_tree_display.py --input-dir $(DATA_ROOT)/repository_inputs/figures/b_lectin_domain_tree --out $(BUILD)/figures/b_lectin_domain_tree/primary --rebuild $(BUILD)/figures/b_lectin_domain_tree/internal_rebuild
source-panel-public:
	$(PYTHON) -B -c 'from pathlib import Path; import sys; sys.path.insert(0,"scripts/release"); from source_panel import build_source_panel; build_source_panel(Path("$(BUILD)/table_s2/candidate/data/source_panel.tsv"),Path("$(BUILD)/public_source_panel.tsv"))'

.PHONY: s4-architecture-replay
s4-architecture-replay:
	test ! -e $(BUILD)/s4_architecture_replay/architecture.tsv
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -B scripts/figures/accessory_domain_trees/build_architecture.py --mapping $(DATA_ROOT)/repository_inputs/analysis/s4_architecture/full_sequence_mapping.tsv --fasta $(DATA_ROOT)/repository_inputs/analysis/s4_architecture/full_sequences.faa --domtbl $(DATA_ROOT)/repository_inputs/analysis/s4_architecture/pfam_hits.domtblout --out $(BUILD)/s4_architecture_replay/architecture.tsv
	cmp $(BUILD)/s4_architecture_replay/architecture.tsv $(DATA_ROOT)/repository_inputs/figures/accessory_domain_trees/architecture.tsv
	@echo S4_ARCHITECTURE_REPLAY_PASS
