# Bacterial gene transfers and herbivory in cynipoid wasps

Code for *Multiple ancient gene transfers from bacteria facilitated herbivory in
cynipoid wasps*. This repository holds the code that turns the analysis results
into the paper's figures, tables, supplement and manuscript PDF.

All genomes, transcriptomes and sequencing reads used in the paper come from
public archives; accessions are listed in Table S1. The Figshare record
[10.6084/m9.figshare.33265659](https://doi.org/10.6084/m9.figshare.33265659)
is not new primary data. It is a set of intermediate files derived from those
public sources (protein sets, alignments, trees, expression values, synteny
and analysis tables) so the paper can be rebuilt without repeating
genome annotation, read mapping or tree inference.

## What's here

| Folder | Contents |
|---|---|
| `scripts/analysis/` | Analysis code: domain trimming, introns and GC3, catalytic sites, selection, synteny, chromosome gains |
| `scripts/figures/` | One plotting script per figure (table below) |
| `scripts/tables/` | Builders for the supplementary tables |
| `scripts/members/` | Builds the master list of gene copies (species, family, chromosome) that the figures and tables share |
| `scripts/release/` | Assembles the manuscript and supplement PDFs and runs consistency checks |
| `config/` | Settings and lookup tables, including the map from each input file to its Figshare path and checksum (`logical_inputs.tsv`) |
| `authority/` | Manuscript, supplement and caption text that the build typesets |
| `docs/` | Method notes for individual analyses |
| `environment/` | Python and R package versions |

## Where each figure comes from

| Figure | Content | Script (`scripts/figures/`) |
|---|---|---|
| 1, S2 | Species tree with gene counts | `tree_counts.py` |
| 2, S3 | Origin trees with bacterial and eukaryotic outgroups | `origin_trees.R` |
| 3 | Intron positions | `introns.py` |
| 4, S5 | Cynipoid locus trees | `locus_trees.R`, `complete_locus_trees.R` |
| 5 | B-lectin locus tree | `b_lectin_tree.R` |
| 6 | Expression across tissues | `expression.py` |
| S4 | GH30 and PI-PLC-X trees | `accessory_domain_trees/` |
| S6 | B-lectin domain tree | `b_lectin_domain_tree/` |

Figure S1 (*Euceroptres* genome) is a supplied image. Some script and folder
names carry older figure or table numbers; `config/figure_numbers.tsv` and
`config/table_numbers.tsv` give the current numbering.

Trees and selection results were computed once (MAFFT, IQ-TREE, HyPhy) and are
read from the Figshare data; the build does not rerun them.

## Rebuilding the paper

Download and unpack the Figshare archive, then:

```bash
export PCWDE_DATA_ROOT=/path/to/pcwde_supporting_data
make verify       # check the data against the recorded checksums
make manuscript   # build every figure, table, the supplement and the manuscript under build/
```

Needs Python 3, R and XeLaTeX; package versions are in `environment/`.
Details are below.

## License

CC BY 4.0 (`LICENSE`). Cite as in `CITATION.cff`.

---

## Reproduction details (for automated agents and the very thorough)

**Inputs.** `config/logical_inputs.tsv` lists the staged data inputs the build reads, with
its archive-relative path and SHA-256. The archive's `MANIFEST.tsv` and
`SHA256SUMS.txt` cover every packaged file. Verify the download with
`sha256sum -c ARCHIVE.sha256` before unpacking. `make verify` checks
inventories, checksums and the code/data boundary; with an extracted archive it
also checks `PROJECTED_INPUTS.tsv`.

**Environment.** Python and R versions are locked in `environment/python-lock.txt`
and `environment/r-lock.txt`; `make manuscript` checks both first
(`scripts/release/check_environment.py`). Override interpreters with the Make
variables `PYTHON` and `RSCRIPT`; `BUILD` changes the output directory (default
`build/`). `scripts/release/local_env.bash LOCALDIR PYTHON_ENV R_ENV` writes
wrappers with node-local font caches, for shared filesystems. `kpsewhich` and
Ghostscript must be on `PATH`.

**Build order.** `make manuscript` runs, in order: environment and closure
self-tests; the gene-copy table (`scripts/members/build_copy_table.py`,
`build_locus_tips.py`); display pruning of the full origin trees
(`scripts/trees/prune_origin_display_trees.R`, tip lists in
`config/origin_display_tips.tsv`); CAZyme rank analysis
(`make cazy-science`); intron and GC3 summaries (`make intron-gc3`); all figures
(`scripts/release/render_all.bash`); chromosome-gain analysis
(`make ecami-chromosome-gain`); table
builders (`scripts/tables/`); the manuscript (`scripts/release/build_manuscript.py`)
and supplement (`build_supplement.py`); and the supporting-data package
(`make data-deposit`). Numbers quoted in the text are computed once
(`build/numbers.json`) and the build fails if a figure, table or text value
disagrees with them.

**Closure check.** `scripts/release/check_closure.py` requires every file under
`repository_inputs/` to be read by the build or declared in
`config/closure_declared.tsv`, and every deposit file to match
`config/deposit_layout.tsv`. Reads are detected from access times, so the data
root must be on a filesystem that updates atime; the check fails otherwise.
`scripts/release/test_check_closure.py` runs first to confirm that unrelated,
altered and stale files are rejected.

**What is not rerun.** Genome annotation, transcriptome assembly, read
quantification, MAFFT/IQ-TREE inference and HyPhy are upstream; their outputs
are in the archive. The full transfer-context trees (427/364/444 tips for
GH5_2/PL1_8/PL4_4) are in `trees/`, alignments in `alignments/`, unaligned
inputs in `sequences/`; `sidecars/origin_context_full_tip_mapping.tsv` maps
alignment IDs to tree tips. Scripts under `scripts/analysis/` are the code that
produced the deposited intermediates; some large upstream inputs (genome
assemblies, GFFs) are not redistributed, so these are reference
implementations, not a single raw-data-to-paper command.
`docs/methods_code_crosswalk.tsv` links each Methods step to its code and input,
and `docs/READER_DAG.md` shows the input/output graph.

**Domain preprocessing (standalone).**

```bash
python3 -B scripts/analysis/preprocess_domains.py --protein-dir /path/to/domain_inputs --output-dir build/proteins
```

The output directory must be new and outside the input directory. Overlapping
hits are resolved by best independent E-value. On the frozen inputs this gives
426/363/443 entries versus 427/364/444 in the historical alignments; the build
uses the historical alignments and trees. Tests:
`python3 -B scripts/analysis/test_preprocess_domains.py`.

No command here uploads, publishes, pushes or writes to a database.
