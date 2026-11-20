# Reproduction graph

The repository separates computational results from display. The paired supporting-data archive
holds the scientific inputs and intermediates; Git holds their consumers. The
local package will be published at the reserved Figshare DOI with the manuscript.

```text
Figshare archive
  └── pcwde_supporting_data/repository_inputs/
        ├── figures/ ───────────────> scripts/figures/* ──> build/figures/
        ├── tables/accepted/ ───────> scripts/tables/* ───> build/tables/
        ├── tables/authority/ ──────> table/analysis validation scripts
        ├── analysis/domain_preparation/ -> preprocess_domains.py -> build/proteins/
        └── deposit/ ───────────────> supporting-data package validation

config/logical_inputs.tsv + repository_inputs/
  └── make verify

verified inputs
  └── make manuscript
        ├── scripts/release/render_all.bash
        ├── scripts/tables/render_tables.py
        ├── scripts/release/build_manuscript.py
        └── scripts/release/build_data_deposit.py
```

## Upstream analytical code

| Analysis | Code | Principal deposited input | Output unit |
|---|---|---|---|
| Domain trimming and fusion filtering | `scripts/analysis/preprocess_domains.py`, `fusion_filter.py` | proteins and HMMER domtblout | trimmed protein FASTA |
| Intron and GC3 analysis | `scripts/analysis/introns/build_*.py`, `02_map_introns.py`, `03_summarize_introns.py`, `04_build_table_s5.py`, `summarize_gc3.py` | frozen per-genome host metrics, focal CDS segments and nucleotide sequences, accepted annotation proteins, splice junctions and bacterial comparators; member and B-lectin rosters; genome GFFs and assemblies stay upstream | Figure 3, Table S5, intron and GC3 deposit tables |
| Catalytic-site conservation | `scripts/analysis/catalytic/rebuild_catalytic_sites.py` | accepted family sequences and references | per-sequence residue calls |
| Selection preparation and parsing | `scripts/analysis/selection/` | proteins, CDS, metadata, HyPhy outputs | codon alignments and selection summaries |
| Final-49 synteny packaging | `scripts/analysis/synteny/` | locus, block, and genome-state tables | public synteny tables |

The paper build does not rerun MAFFT/IQ-TREE inference, HyPhy, genome
annotation, transcriptome assembly, or read quantification. Their accepted
outputs and identifiers are deposited so the reported results and displays can
be inspected without implying that every upstream historical operation is one
portable command.

Make targets and domain preprocessing write beneath `build/` by default.
`BUILD` and preprocessing `--output-dir` provide explicit output overrides;
preprocessing refuses existing or input-overlapping output directories.

## Reader input forms

An unprojected source root is checked against every current projected-input hash in
`config/logical_inputs.tsv`. An extracted `pcwde_supporting_data` archive carries
`PROJECTED_INPUTS.tsv`; `make verify` checks the fixed source-manifest identity,
all 479 source paths, and projected hashes and closure before its scientific
inventory checks. `make stage-check` then compares the extracted projected tree
with a newly staged projected tree. A missing projected manifest selects raw
source validation and cannot validate transformed GFF comments.
