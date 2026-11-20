#!/bin/bash
set -euo pipefail

: "${RENDER_REPO:?RENDER_REPO is required}"
: "${RENDER_DATA:?RENDER_DATA is required}"
: "${RENDER_OUT:?RENDER_OUT is required}"
: "${RENDER_PYTHON:?RENDER_PYTHON is required}"
: "${RENDER_RSCRIPT:?RENDER_RSCRIPT is required}"

if [ ! -d "$RENDER_REPO/scripts/figures" ] || [ ! -d "$RENDER_DATA" ]; then
  echo "missing repository or data-root interface" >&2
  exit 2
fi
if [ -e "$RENDER_OUT" ]; then
  echo "refusing existing render output: $RENDER_OUT" >&2
  exit 2
fi
mkdir -p "$RENDER_OUT"
SCRIPTS="$RENDER_REPO/scripts/figures"

: "${RENDER_COPY_TABLE:?RENDER_COPY_TABLE is required (built by scripts/members/build_copy_table.py)}"
COPY_TABLE="$RENDER_COPY_TABLE"
: "${RENDER_LOCUS_TIPS:?RENDER_LOCUS_TIPS is required (built by scripts/members/build_locus_tips.py)}"
GENOME_ROSTER="$RENDER_DATA/../tables/authority/synteny_current/final_49_roster.tsv"
PCWDE_PRESENCE_COUNTS_INPUT_DIR="$RENDER_DATA/presence_counts" \
PCWDE_CHROMOSOME_COUNTS_INPUT_DIR="$RENDER_DATA/chromosome_counts" \
PCWDE_FIGURE_SHARED_INPUT_DIR="$RENDER_DATA/shared" \
PCWDE_COPY_TABLE="$COPY_TABLE" \
PCWDE_TBLASTN_PRESENCE="$RENDER_DATA/../shared/tblastn_presence.tsv" \
PCWDE_GENOME_ROSTER="$GENOME_ROSTER" \
PCWDE_GH1_DERIVED_COUNTS="$RENDER_REPO_BUILD/cazy_counts/derived_gh1e72_gene_counts.tsv" \
PCWDE_TREE_COUNTS_OUTPUT_DIR="$RENDER_OUT/tree_counts" \
"$RENDER_PYTHON" "$SCRIPTS/tree_counts.py"
PCWDE_ORIGIN_TREES_INPUT_DIR="$RENDER_DATA/origin_trees" PCWDE_ORIGIN_DISPLAY_TREE_DIR="$RENDER_REPO_BUILD/trees" PCWDE_ORIGIN_TREES_OUTPUT_DIR="$RENDER_OUT/origin_trees" "$RENDER_RSCRIPT" "$SCRIPTS/origin_trees.R"
"$RENDER_PYTHON" "$SCRIPTS/origin_trees_full/merge_panels.py" \
  --panel "$RENDER_OUT/origin_trees/origin_trees_full_GH5_2_ggtree.pdf" \
  --panel "$RENDER_OUT/origin_trees/origin_trees_full_PL1_8_ggtree.pdf" \
  --panel "$RENDER_OUT/origin_trees/origin_trees_full_PL4_4_ggtree.pdf" \
  --output "$RENDER_OUT/origin_trees/origin_trees_full.pdf" --work-dir "$RENDER_OUT/origin_trees/s1_composer_work"
PCWDE_INTRONS_INPUT_DIR="${RENDER_INTRONS_INPUT_DIR:-$RENDER_DATA/introns}" PCWDE_INTRONS_OUTPUT_DIR="$RENDER_OUT/introns" "$RENDER_PYTHON" "$SCRIPTS/introns.py"
PCWDE_LOCUS_TIP_CROSSWALK="$RENDER_LOCUS_TIPS" PCWDE_LOCUS_TREES_INPUT_DIR="$RENDER_DATA/locus_trees" PCWDE_LOCUS_TREES_COMPLETE_INPUT_DIR="$RENDER_DATA/complete_locus_trees" PCWDE_LOCUS_TREES_OUTPUT_DIR="$RENDER_OUT/locus_trees" "$RENDER_RSCRIPT" "$SCRIPTS/locus_trees.R"
PCWDE_COPY_TABLE="$COPY_TABLE" PCWDE_B_LECTIN_TREE_INTERFACE="$RENDER_DATA/b_lectin_tree/interface.tsv" PCWDE_B_LECTIN_TREE_TREE_FILE="$RENDER_DATA/b_lectin_tree/tree.treefile" PCWDE_B_LECTIN_TREE_ALIGNMENT_FILE="$RENDER_DATA/b_lectin_tree/domain_lane_alignment.tsv" PCWDE_B_LECTIN_TREE_OUTPUT_DIR="$RENDER_OUT/b_lectin_tree" "$RENDER_RSCRIPT" "$SCRIPTS/b_lectin_tree.R"
PCWDE_FIGURE_DATA_ROOT="$RENDER_DATA" PCWDE_EXPRESSION_OUTPUT_DIR="$RENDER_OUT/expression" "$RENDER_PYTHON" "$SCRIPTS/expression.py"
PCWDE_ACCESSORY_DOMAIN_TREES_INPUT_DIR="$RENDER_DATA/accessory_domain_trees" PCWDE_ACCESSORY_DOMAIN_TREES_OUTPUT_DIR="$RENDER_OUT/accessory_domain_trees" "$RENDER_RSCRIPT" "$SCRIPTS/accessory_domain_trees/render_catalytic_trees.R"
"$RENDER_PYTHON" "$SCRIPTS/merge_panels.py" --panel "$RENDER_OUT/accessory_domain_trees/accessory_domain_trees_GH30_3_ggtree.pdf" --panel "$RENDER_OUT/accessory_domain_trees/accessory_domain_trees_PI_PLC_X_ggtree.pdf" --expected-pages 2 --output "$RENDER_OUT/accessory_domain_trees/accessory_domain_trees.pdf"
PCWDE_LOCUS_TIP_CROSSWALK="$RENDER_LOCUS_TIPS" PCWDE_COMPLETE_LOCUS_TREES_INPUT_DIR="$RENDER_DATA/complete_locus_trees" PCWDE_COMPLETE_LOCUS_TREES_OUTPUT_DIR="$RENDER_OUT/complete_locus_trees" "$RENDER_RSCRIPT" "$SCRIPTS/complete_locus_trees.R"
"$RENDER_RSCRIPT" "$SCRIPTS/validate_complete_locus_display_root.R" "$RENDER_DATA/complete_locus_trees" "$RENDER_OUT/complete_locus_trees/evidence"
"$RENDER_PYTHON" "$SCRIPTS/merge_panels.py" --panel "$RENDER_OUT/complete_locus_trees/complete_locus_trees_GH5_2_ggtree.pdf" --panel "$RENDER_OUT/complete_locus_trees/complete_locus_trees_PL1_8_ggtree.pdf" --panel "$RENDER_OUT/complete_locus_trees/complete_locus_trees_PL4_4_ggtree.pdf" --expected-pages 3 --output "$RENDER_OUT/complete_locus_trees/complete_locus_trees.pdf"
"$RENDER_PYTHON" "$SCRIPTS/b_lectin_domain_tree/build_b_lectin_tree_display.py" --input-dir "$RENDER_DATA/b_lectin_domain_tree" --out "$RENDER_OUT/b_lectin_domain_tree/primary" --rebuild "$RENDER_OUT/b_lectin_domain_tree/internal_rebuild"

echo "RENDER_COMPONENT_PASS output=$RENDER_OUT"
