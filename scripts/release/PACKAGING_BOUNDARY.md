# Shared public-input packaging boundary

`data_root.py` verifies the immutable source logical-input manifest and the distinct projected-input manifest. `public_inputs.py` applies the eight declared public transformations and validates the projected repository inventory. The production builder calls it with the fixed source-manifest digest; the tiny test supplies a synthetic digest only at the lower-level function boundary. `package_archive.py` writes the same deterministic tar/gzip member order and metadata as the prior builder and verifies the checksum, member inventory, modes, containment and bytes. The production validator calls this verifier while retaining all scientific checks. The normal `stage.py --compare` requires both roots to have the same form and validates each by that form's contract (projected manifest, or raw source pins) before comparing bytes; `make stage-check` then projects the restaged root and validates its projected manifest.

`test_packaging_command_fixture.py` creates the 518 committed logical paths with tiny synthetic content and structurally meaningful GFF, CDS and figure inventory inputs. It validates source pins, projects, validates the projected manifest, writes and verifies the archive, extracts it, and invokes the real stage comparison command. It tests changed, missing and extra source files, wrong manifest digest, unexpected omission, altered projected checksum and corrupted archive. The focused projection and pinned-hybrid tests retain their separate negative cases, including private GFF feature lines and the fixed production-wrapper digest. This fixture qualifies packaging only. Production scientific validation, real-input A/B fixture and make stage-check remain required on the full data.

The reader verifier selects raw or projected form by the projected manifest.
Raw form requires original source hashes; projected form requires fixed
source-manifest identity, complete source paths, and projected hashes and closure.
The repository-only policy switch uses the same repository checks as the
combined `make verify` entry point, without opening the data root.
