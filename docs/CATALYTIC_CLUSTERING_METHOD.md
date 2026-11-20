# Catalytic clustering input reconciliation

Table S9 uses the accepted 2026-10-14 tree reader, Fitch score, per-site and
per-family seeds, 10,000 permutations, and callability-mask-preserving family
null. The producer is `scripts/analysis/catalytic/clustering/reconcile.py`.
The locked intermediate `LOCKED_ANALYSIS.json` was rerun with the same accepted
full origin trees, metadata, and source-to-tree crosswalk. Its only changed input
was the catalytic call table after the current-origin membership filter used by
Table S10. Removed Bassettia tree tips remain in the tree and receive unknown
states, exactly as the original method treats tips without accepted calls.

Original producer SHA256: reconcile.py
`49406aeeca0cfbcc34100ce20cde1e8e4ca58818d366eaeb54032b9979c34048`;
analyze.py (tree-reader source)
`cf81c12b5bd9c0dfb41e7d96c3e5db9f065488357b1dd5bbb80874674ce8da16`.
Crosswalk SHA256:
`7553b612d27e14dcbaab6e7fdac0f273e9c4c9383d3aef6b5ff6ec5f950a552e`.
Locked call TSV SHA256:
`e3b0cfdc8c69528b54ad5cf74ef1cbf9c3965514a6adbcbe35bf0d26e00d608f`.
The producer's old-call regression matched every original statistical and coverage
field exactly. The accepted JSON adds a later sequence-unit audit field.
The build checks S9 retained/callable counts against S10 expected/observable
counts for each catalytic site and fails on disagreement.
