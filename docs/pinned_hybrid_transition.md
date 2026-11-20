# Public input projection of the pinned intron stage

The intron rerun stages 518 files described by `config/logical_inputs.tsv`. The source manifest is checked before deposit construction and is never rewritten by public projection. The current pinned stage contains public TSV fields and source labels from the earlier projection, with private annotation paths still present in 33 GFF comment lines. It has the previously omitted provenance files absent.

The deposit producer accepts this state only when the source manifest has its committed digest, every staged path, byte count and SHA256 matches that manifest, the public TSV headers and omission state match the declared transition, and the only remaining public hazards are private paths in declared GFF comments. It then replaces those comments, refreshes the figure byte inventory without changing row order or other fields, and rebuilds the repository input inventory. GFF feature lines remain byte-identical. The source manifest still describes the staged inputs; `PROJECTED_INPUTS.tsv` records the separate hashes of the transformed public inputs.

The deposit validator checks the projected manifest and its complete path closure. The stage check calls the same projection producer on its independently staged source inputs before comparing them to the archive. Fully raw projection and fully projected idempotence retain their existing paths; other partial states are rejected.
`PROJECTED_INPUTS.tsv` is integrity metadata with the exact three-column
`relative_path`, `bytes`, `sha256` schema. The deposit validator validates its
paths, hashes, and closure through `validate_projected_inputs` before applying
the general ban on scientific TSV row hashes. No other TSV gains this exception.
