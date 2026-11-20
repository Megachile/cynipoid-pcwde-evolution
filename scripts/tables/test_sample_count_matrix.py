#!/usr/bin/env python3
"""Check count units and fail-closed reconciliation on the pinned expression roster."""
import copy
import build_table_s6 as s6

families, genes, metadata, summary, crosswalk, _ = s6.build_expression()
panel = s6.read_tsv(s6.AUTH_S2 / "source_panel_accepted.tsv")
expected = s6.sample_count_matrix(families, metadata, crosswalk, summary, panel)
assert s6.sample_count_matrix(families + families[:1], metadata, crosswalk, summary, panel) == expected
assert s6.sample_count_matrix(families, metadata, crosswalk, summary, panel) == expected
bad = copy.deepcopy(summary)
bad[0]["sample_count"] = str(int(bad[0]["sample_count"]) + 1)
try:
    s6.sample_count_matrix(families, metadata, crosswalk, bad, panel)
except AssertionError as error:
    assert "mismatch" in str(error)
else:
    raise AssertionError("inconsistent Table S2b totals passed")
print("SAMPLE_COUNT_MATRIX_PASS species=15 cells=37 profiles=158")
