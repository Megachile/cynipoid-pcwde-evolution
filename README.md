# PCWDE evolution in cynipoid wasps

This repository contains the analysis and rendering code for *Multiple ancient
gene transfers from bacteria facilitated herbivory in cynipoid wasps*. The
paired supporting-data archive is reserved for Figshare DOI
`10.6084/m9.figshare.33265659` and will be released with the manuscript.

Git contains code, configuration, and the manuscript source. Scientific data,
alignments, trees, analysis tables, and figure inputs are in the paired
supporting-data archive, which is not yet publicly released. The file
`config/logical_inputs.tsv` connects each input used by the code to its
Figshare-relative path and checksum.

## Reproduce the paper outputs

1. Obtain and unpack the paired supporting-data archive (local until publication).
2. Clone this repository.
3. Set `PCWDE_DATA_ROOT` to the unpacked `pcwde_supporting_data` directory.
4. Run:

```bash
make verify
make manuscript
```

`make manuscript` also checks build closure (`scripts/release/check_closure.py`): every file under
`repository_inputs/` must be read by the build or declared in `config/closure_declared.tsv` (frozen member
layer, cached trees, out-of-family inputs, inventories, or a `known_gap` with its planned fix), and every
reader-deposit file must match `config/deposit_layout.tsv`: a byte copy of the input at its mapped path, or a
declared product of a named committed script written during this build. `scripts/release/test_check_closure.py`
runs first and proves an unrelated file, an altered copy and a stale generated file all fail. Reads are detected from file access times
on the extracted data root; the check fails if the filesystem does not track them. Known gaps are printed
in every build log.

`make verify` checks file inventories, SHA-256 checksums, and the Git/Figshare
boundary. `make manuscript` renders Figures 1–6, Figures S1–S6, Tables S1–S10 (including S2a/S2b),
and the manuscript PDF beneath `build/`.

The rendering environment requires Python 3, R, XeLaTeX and `kpsewhich`. Python and R package versions are recorded under
`environment/`. Executable locations can be supplied with the Make variables
`PYTHON` and `RSCRIPT`. `scripts/release/local_env.bash LOCALDIR PYTHON_ENV R_ENV` writes
wrappers with node-local font caches; `make manuscript` first checks both
environments against the locks (`scripts/release/check_environment.py`).

## Analysis scope

The paper outputs are reproducible from the paired archived intermediates. The
repository also preserves the current upstream analysis scripts used for
domain processing, introns, synteny, catalytic sites, and selection. Some large
or externally sourced upstream inputs are not redistributed, so these scripts
are reference implementations rather than a promised raw-data-to-paper
workflow. Expensive MAFFT, IQ-TREE, and HyPhy results are not rerun by the
manuscript target. Figure renderers consume deposited trees, and the synteny
release code consumes the deposited final-49 locus and block tables.

The input/output graph is in `docs/READER_DAG.md`, and method-to-code links and
input availability are in `docs/methods_code_crosswalk.tsv`.

The Make targets write beneath `build/` by default (`BUILD` overrides it).
For standalone domain preprocessing, provide the family protein FASTAs and
original `hmmsearch` domain tables in a separate input directory:

```bash
python3 -B scripts/analysis/preprocess_domains.py --protein-dir /path/to/domain_inputs --output-dir build/proteins
```

`--output-dir` defaults to `build/proteins` relative to the working directory.
It must be a new directory, separate from and outside the input directory;
existing output directories are rejected. All six required input files are
checked before output creation. `PCWDE_PROTEIN_DIR` supplies the input default;
otherwise the script uses `repository_inputs/analysis/domain_preparation`
beneath `PCWDE_DATA_ROOT` (or `figshare_data`). The script preserves its inputs.

Overlapping domain hits are considered in increasing independent domain E-value,
with deterministic coordinate-based tie breaks. A retained hit blocks any other
hit sharing an alignment residue (one-based inclusive coordinates). Retained
hits are returned in sequence order. `--overlap-rule legacy-diagnostic` is only
for comparison with historical start-first preprocessing. On the frozen inputs,
the corrected rule yields 426/363/443 entries rather than 427/364/444 for
GH5_2/PL1_8/PL4_4. These differences have not been adopted as new scientific
inputs: `make manuscript` continues to use the fixed archived alignments and
trees. Details of the full inference evidence are in `docs/DATA_ACCESS.md`.
Run the interface, parser and overlap tests with
`python3 -B scripts/analysis/test_preprocess_domains.py`.

No provided command uploads, publishes, commits, pushes, or modifies a database.
