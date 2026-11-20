# Current catalytic-domain membership

`make manuscript` runs `scripts/analysis/catalytic/subset_current_origin.py` before
building the catalytic reader table. The producer reads the pinned catalytic
per-sequence residue calls and summary in
`repository_inputs/tables/authority/catalytic_current/` and the current full
origin-domain FASTAs in `repository_inputs/deposit/base/sequences/`. It takes
the exact cynipoid IDs in each current full origin FASTA as membership authority.
The recorded removed-tip sidecar must identify every excluded ID. The producer
fails if the origin and catalytic IDs differ in any other way.

The pinned catalytic calls were read from reference-augmented alignments. The
producer filters their rows and recomputes each site's total, observed, missing,
expected-residue and percentage counts. Each retained residue call and alignment
column is copied unchanged; no alignment is recomputed. The historical source
TSVs remain frozen inputs for provenance, while `build/catalytic_current/` is the
single current intermediate consumed by Table S10 and its deposit products.
The positive-selection input and catalytic-clustering permutation test are
unchanged.

Pinned input SHA256 values for this integration build:

| Input | SHA256 |
|---|---|
| catalytic per-sequence calls | `ef6dab963e49fbdfba266beebcd7f61b159ea3430e62cb5a1f205fd3242a7c09` |
| catalytic old summary | `d1292fef425268aec53eb7f982237e50aa3f091c19de1287aa0248321adb09ee` |
| GH5_2 full origin FASTA | `bf94586271cde876bfc707bc566fc7aa55770d300bd54194ae0a2906f2f5b2ed` |
| PL1_8 full origin FASTA | `1d4d60f3caa0cad4aac66573a75d6921ed8e06d0f3891780d8d37393d3a01256` |
| PL4_4 full origin FASTA | `77c42b56c356258d782365e6de2be74ab475fc984edaecde912beefdd1bfc06d` |
| removed-tip record | `8a1120c74c52773691939ecd449ae70aef025fe61f3d1d92ca62c06bc10a6f13` |
