# Revised deposit inputs

`stage.py` produces the revised `repository_inputs/` from the submitted
Figshare deposit and a fixed set of upstream analysis outputs. It is the only
way deposit inputs change; `make` then builds every figure and table from the
result.

Steps, in order:

| Step | What changes | Upstream |
|---|---|---|
| `bp_cb_pruning` | Removes the three PL1_8/PL4_4 and one B-lectin *Bassettia pallida* transcripts that the fixed removal rule assigns to *Callirhytis quercusbatatoides* origin (per-tip values in `removed_cb_origin_bp_tips.tsv`): records dropped from sequences, alignments, sidecars and domain-preparation files (`drop_rows.py`), tips pruned from trees without re-inference (`prune_trees.py`, `config/bp_cb_prune_plan.tsv`) | 2026-10-19 pruning run: removal record, two method notes, Figure S4 inputs |
| `copy_table` | Adds `shared/copy_placements.tsv`, one row per placed locus listing its gene copies | 2026-10-21 placement-conformance table (internal `evidence` column dropped) |
| `computed_tables` | Removes saved Tables S2–S4; `make` computes them | — |
| `locus_trees` | Replaces Figure 4 / S3 inputs; adds `shared/locus_tips.tsv` | 2026-10-23 locus-tree rebuild; 2026-10-28 roster with the same-gene BK venom slot |
| `reader_locus` | Replaces reader-deposit complete-locus alignments, sequences, trees and tip labels | 2026-10-23 locus-tree rebuild |
| `gh1_counts` | Adds Figure 1 GH1_e72 counts | Figure 1 GH1_e72 scan |
| `tblastn_presence` | Adds `shared/tblastn_presence.tsv`: the stated tBLASTn filters applied to every genome/family pair without an accepted member; `admitted_validated_rescue` rows count as present (presence rule b) but carry no gene copy | 2026-10-24 presence-count audit |
| `member_layer` | Adds the frozen member layer `members/member_table.tsv` (one row per gene copy: identity, the accepted placement run's per-locus placement carried on each copy, and for the 46 curated B-lectin units their curated unit id, scaffold and route-B evidence from `beta_units.tsv`) and `members/member_order.tsv` (each three-family copy's positions in the locked ordered gene list: scaffold, rank, position, and coordinates where the locked annotation has them, otherwise the placement run's reconstructed-scaffold inputs, with `coordinate_status`) | Step 2a run (locked final gene order, member attributes); placement run reconstructed-scaffold inputs |
| `synteny_existing_identity` | Adds `members/synteny_existing_identity.tsv`: a per-copy disposition linking accepted members to pre-existing ordered-gene identities. It matches direct IDs, then uniquely identical amino-acid sequences in the genome's pinned annotation proteome. Missing or ambiguous identity never creates a new order row. | Frozen member layer and original gene order; SHA-pinned accepted B-lectin sequences and annotation proteomes in `config/intron_sources.tsv` |
| `corrected_synteny` | Copies the frozen order, affected-scaffold order and species labels; reruns the native backend and accepted recovery with the eligible-copy identity table. Derived membership, block summary and block-by-genome states are checked against `config/synteny_derived_pins.tsv` after the pin-generation build. | `config/synteny_existing_sources.tsv`, the frozen member layer and source-derived identity disposition |
| `derived_copy_table` | Removes `shared/copy_placements.tsv`; `make` rebuilds it from `members/` with `scripts/members/build_copy_table.py` and the locus rule in `scripts/members/locus_rule.py` | — |
| `derived_locus_tips` | Replaces `shared/locus_tips.tsv` with `members/locus_tip_crosswalk.tsv`: the tree-intrinsic columns (display tip, raw tip, source kind, source protein; species only for transcript tips). `make` derives genome tips' species and BK chromosome from the built copy table (`scripts/members/build_locus_tips.py`) | — |
| `retired_figures` | Removes the synteny figure's inputs (`figures/figure4/`); the figure was dropped from the paper | — |
| `stable_names` | Renames figure input folders to stable figure ids (e.g. `figures/figure3/` → `figures/origin_trees/`, `tables/figure5_current/` → `tables/introns_current/`); numbers are applied only at assembly from `config/figure_numbers.tsv` | — |
| `reader_cleanup` | Removes 32 inputs: 26 that no build step reads (legacy chromosome-count overlays, two legacy shared tables, a header-only audit, a style PDF, four legacy synteny tables, the older intron set), the superseded saved Tables S6/S7, and the stale synteny-locus reader tables `deposit/base/data/family_membership.tsv` and `chromosome_placement_evidence.tsv` (the reader deposit ships the build's gene-copy versions) | — |

Every upstream file is pinned by sha256 in `config/deposit_sources.tsv`; a
changed or missing upstream file stops the run. Alignment, tree inference and
placement are not rerun here.

    make stage-check BASE=<submitted data root> WORKSPACE=<analysis workspace> DATA_ROOT=<revised data root>

reruns all steps into the build directory and fails unless every file under
`repository_inputs/` matches the revised deposit byte for byte.
`stage.py --pin DATA_ROOT` updates `config/logical_inputs.tsv` after a deposit change.

## Member layer and the locus rule

Gene copies are grouped into loci in `make`: same-family copies on one scaffold join when at most five genes lie between them in the ordered gene list (every gene counts, member or not), and a copy present on several scaffold segments links its loci (`scripts/members/locus_rule.py`). This is the joining logic of `find_loci_segment_aware`/`_make_locus` in `scripts/analysis/synteny/01_build_synteny.py` as run by placement job 1344052 (2026-10-20), without that module's live chromosome lookup; on the locked inputs both give the same 369 three-family loci (membership, locus ids, anchor scaffolds) as the accepted copy table. The ordered positions come from the locked annotations listed in `config/upstream_locked_inputs.tsv` (step 2a, all matched to what the placement run recorded); the full ordered gene list is not deposited. B-lectin units are curated (route-B pairing rule) and frozen, not regrouped.

### Unplaced *Eubothrus* PL1_8 rescue

The shared copy producer reads the checksum-pinned `shared/tblastn_presence.tsv` and
`shared/fig5_validated_raw_locus_rescue.tsv`. It requires their unique *Eubothrus*
PL1_8 records to agree on `NODE_38372_length_2218_cov_4.820177:1323-346`.
The validated raw locus has no Helixer model or ordered-gene position, so the
producer appends exactly one copy outside the annotation-based locus grouping.
Its BK placement is unplaced under RagTag, and its `provenance` field names both
inputs. The frozen corrected-synteny cohort excludes this raw locus because it
was not available to the ordered-gene/synteny inference; the placement table
includes it as an unplaced copy.
