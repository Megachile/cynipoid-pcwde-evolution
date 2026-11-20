# Refreshing supplement prose

Adam writes the supplement prose in the editable author document. The build never types result numbers from that document: it renders the repository template `authority/supplement/supplement_source.md`, in which every result number is a token `{{key}}` or `{{key|format}}` resolved from the build's `numbers.json` (formats: `comma`, `f0`/`f1`/`f2`, `word`/`Word`). `scripts/release/produce_supplement_sections.py` renders the template and fails the build on an unknown or unresolved token, any `[CHECK` text, any numeral that is neither a token nor a row in `config/supplement_literals.tsv`, a literal row that no longer matches the template, and a parameter literal that disagrees with its producer constant (`check_key`). `scripts/release/test_supplement_template.py` holds negative fixtures (a typed stale result must fail). The pinned export the template was derived from is recorded as a class U input in `config/upstream_locked_inputs.tsv` (currently `authority/supplement/supplement_draft_20261118T0358Z.md`); the historical exports are retained only for provenance.

`config/supplement_literals.tsv` lists the numerals that stay literal: method parameters, software versions, dates, citations and fixed protocol values (columns: section, numeral, context, class, check_key). Result numbers do not belong there.

For later author edits:

1. Sync and export the document as Markdown on the droplet (`lib/gdoc_sync.py`); keep the private document identifier out of Git and public artifacts.
2. Run `python3 scripts/supplement/doc_edit_report.py --repo . --new <export.md>`. It compares the new export with the pinned one paragraph by paragraph and lists changed text and every new numeral. It reports only; it does not apply edits.
3. On a new branch, the terminal ports the prose changes into `supplement_source.md` as a candidate: new result numbers become tokens (adding the key to the producer that owns the fact, computed from the same intermediate the tables use), new parameters/versions/dates/citations become literal rows. Pin the new export's SHA256 in place of the old one.
4. Rebuild with `make paper` from a clean exact-SHA checkout, review the rendered supplement and the closure comparison, prepare a Gallery candidate, and merge only after Adam picks it.

The manuscript is rendered from the pinned current draft `authority/manuscript/manuscript_draft_20261118T0422Z.md` (see `authority/manuscript/SOURCE_AUTHORITY.md`); its prose is verbatim and its numbers are not tokenised. `authority/manuscript/manuscript_current.md` stays unchanged as the submission record.

The current refresh adopts the pinned supplement and 1531Z manuscript. Approved producer-backed numerical projections are declared in `config/author_result_routes.json`; exact paragraph gates enforce author words. No pending source decisions remain.
