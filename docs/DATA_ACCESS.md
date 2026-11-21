# Supporting-data access

Download `pcwde_supporting_data.tar.gz` from
Figshare DOI `10.6084/m9.figshare.33265659`.

After downloading the archive and its `ARCHIVE.sha256` sidecar:

```bash
sha256sum -c ARCHIVE.sha256
tar -xzf pcwde_supporting_data.tar.gz
export PCWDE_DATA_ROOT="$PWD/pcwde_supporting_data"
make verify
```

The archive has one top-level directory. `config/logical_inputs.tsv` maps the
inputs used by repository code to paths and SHA-256 checksums. Within the
archive, `MANIFEST.tsv` and `SHA256SUMS.txt` cover every packaged payload.

The archive includes full transfer-context trees with 427 GH5_2,
364 PL1_8 and 444 PL4_4 tips (365 and 446 were inferred; see below). For each family, `trees/origin_context_FAMILY_full.treefile`
is paired with `alignments/origin_context_FAMILY_full.aln.faa` and the original
unaligned MAFFT input `sequences/origin_context_FAMILY_full.faa`. Here “full”
means complete inference membership, not untrimmed protein length. The
retained `origin_context` subset exports contain 186/95/125 tips. They differ
from the current executable Figure S1 trees (170/67/82 tips) and Figure 3 trees
(29/28/29 tips) under `repository_inputs/figures/origin_trees/`. All of these tip sets
are subsets of the complete inference trees.
`sidecars/origin_context_full_tip_mapping.tsv` maps alignment IDs to IQ-TREE
tip spelling; per-family method records give the original version, model,
seed and bootstrap settings. Original sequence bytes are preserved. Alignment
comparison removes alignment gaps and stop symbols.

`repository_inputs/analysis/domain_preparation/` contains frozen family proteins,
the original hmmsearch hit rows restricted to those proteins, and historical
trimmed FASTAs. Preprocessing requires no new HMM search. The corrected overlap
rule produces diagnostic sequence differences; the manuscript targets use the
fixed historical alignments and trees.

Four *Bassettia pallida* transcript-derived sequences were removed after
inference because they are near-identical to *Callirhytis quercusbatatoides*
sequences. `sidecars/removed_cb_origin_bp_tips.tsv` lists each removed tip with
the tree distances and protein identity that met the removal rule. They were
pruned from every deposited and figure-input tree and dropped as rows from the
matching sequences, alignments, sidecars, figure inputs and domain-preparation
inputs, without realignment or column changes. All other topology, branch
lengths and support values are exactly as inferred, so rerunning inference from
the deposited alignments will not reproduce the deposited trees byte-for-byte.

## Portable upstream source aliases

`config/source_aliases.tsv` binds each `source_inputs/` name to a file in the pinned input tree, with its byte count and SHA256. For upstream restaging, set `PCWDE_SOURCE_INPUTS` to that tree. The resolver rejects missing aliases, changed bytes, duplicate identifiers and paths outside the tree. Matching files reuse `repository_inputs/`; original records needed before projection are frozen under `build_inputs/source_inputs/` and do not enlarge the public deposit. Public reference files identified by Table S1 accessions are not redistributed: NCBI annotation GFFs used only to collapse proteins to genes are replaced by compact protein-to-gene index tables, and outgroup reference proteomes are kept only as SHA-256 pins. The *Belonocnema kinseyi* NCBI GFF and RefSeq proteome stay: intron restaging and synteny member identity read the GFF, and accepted *B. kinseyi* members are matched to the proteome by sequence. `config/public_reference_reductions.tsv` lists each accession, removed file and its original SHA-256; `scripts/deposit_inputs/reduce_public_references.py --derive-index GFF OUT [--bk]` regenerates an index from a downloaded GFF. Original site-specific locations and source hashes are retained in the external build provenance map.
