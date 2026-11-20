# Figure S6 rendering and source-panel projection

Figure S6 is redrawn by `build_b_lectin_tree_display.py` and `render_domains.R`
from the pinned rooted B-lectin domain tree and the 68-row architecture table.
Sequences, alignment, tree, support values, branch lengths and domain coordinates
are retained. The R drawing implementation comes from the accepted pruning
producer; it ladderizes right=FALSE, displays no branch lengths, shows UFBoot
support >=70, and aligns complete protein architectures at the represented
B-lectin domain. No topology or coordinates are recovered from artwork.

The Python wrapper normalizes the established display taxonomy and italicizes
species names, then converts the newly drawn SVG to PDF and PNG. Each render
writes its own terminal order, plot data, expanded domain rectangles and QC.
The historical SVG and QC remain comparison evidence only. They are never
opened by the renderer. `make s6-render` exercises the same producer separately.

The deposit uses every row and column of
`BUILD/table_s2/candidate/data/source_panel.tsv`. `source_panel.py` applies the
existing public metadata projection: Mayetiola destructor GCA_000149185.1 BUSCO
complete becomes 95.7; author_generated repositories become NCBI, NCBI GEO,
USDA PBARC or Study dataset according to their existing identifiers; the public
repository/accession headers become source_repository_or_study and
source_accession_or_identifier. Scientific rows and all other cells retain their
order and values. A missing same-build table is fatal. The saved reader table
is excluded when the base reader directories are copied. `make source-panel-public`
exercises the same transformation independently of deposit packaging.
