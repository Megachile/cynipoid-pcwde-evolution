# Figure renderers

`scripts/release/render_all.bash` is the paper-level entry point. It calls the
individual renderers for Figures 1–8 and Figures S1–S4 using inputs from
`repository_inputs/figures/` in the Figshare deposit.

The tree renderers consume deposited Newick trees and metadata; they do not
rerun alignment or phylogenetic inference.

Figure 6 is a 31/33/30-tip display pruning of the same complete locus trees
used by Figure S3. Its renderer receives both the Figure 6 annotation interface
and the Figure S3 tree interface from the single unpacked data root, maps
support from source-tree splits, and applies only the documented display roots.

All paths are supplied by the release controller through environment variables,
and all generated objects are written beneath the selected build directory.
