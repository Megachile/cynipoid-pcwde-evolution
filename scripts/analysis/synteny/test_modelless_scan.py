#!/usr/bin/env python3
"""Interval-boundary and scope checks for the post-detection model-less scan."""
import sys
import tempfile
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from produce_existing_recovery import BETA, modelless_loci, scan_modelless  # noqa: E402

chd = pd.DataFrame({"genome": ["G"] * 4, "chrom": ["s1"] * 4, "position": ["0", "1", "2", "3"],
                    "protein_id": ["a", "b", "c", "d"], "cynog_id": ["", "", "", ""],
                    "start_bp": ["100", "300", "500", "700"], "end_bp": ["200", "400", "600", "800"]})
base = {"species": "G", "recovered": "yes", "scaffold": "s1", "n_conserved_flanking_cynogs": 2,
        "search_route": "ordered_scaffold", "inapplicable_reason": ""}
matrix = pd.DataFrame([
    {**base, "family": BETA, "block_id": "B1", "genome": "G", "status": "EMPTY", "region_start_rank": 1,
     "region_end_rank": 2, "n_target_units": 0, "target_unit_ids": ""},
    {**base, "family": BETA, "block_id": "B2", "genome": "G", "status": "EMPTY", "region_start_rank": 0,
     "region_end_rank": 0, "n_target_units": 0, "target_unit_ids": ""},
    {**base, "family": BETA, "block_id": "B3", "genome": "G", "status": "N_D", "recovered": "no",
     "region_start_rank": "", "region_end_rank": "", "n_target_units": 0, "target_unit_ids": ""},
    {**base, "family": "GH5_2", "block_id": "B4", "genome": "G", "status": "EMPTY", "region_start_rank": 0,
     "region_end_rank": 3, "n_target_units": 0, "target_unit_ids": ""},
    {**base, "family": BETA, "block_id": "B5", "genome": "G", "status": "PRESENT", "region_start_rank": 0,
     "region_end_rank": 3, "n_target_units": 1, "target_unit_ids": "x"},
])
loci = [{"genome": "G", "family": BETA, "accepted_id": "in", "scaffold": "s1", "start": 300, "end": 600},
        {"genome": "G", "family": BETA, "accepted_id": "edge", "scaffold": "s1", "start": 150, "end": 350}]
out, audit = scan_modelless(matrix, chd, loci)
status = dict(zip(out.block_id, out.status))
assert status == {"B1": "PRESENT", "B2": "EMPTY", "B3": "N_D", "B4": "EMPTY", "B5": "PRESENT"}, status
assert out.set_index("block_id").loc["B1", "target_unit_ids"] == "in"
assert out.set_index("block_id").loc["B5", "target_unit_ids"] == "edge;in;x"
assert out.set_index("block_id").loc["B5", "n_target_units"] == 3
pos = {(r.block_id, r.accepted_id): r.position for r in audit.itertuples()}
assert pos[("B1", "edge")] == "boundary_overlap" and pos[("B2", "in")] == "outside"
assert ("B4", "in") not in pos  # other family never scanned
assert matrix.status.tolist()[0] == "EMPTY"  # input not mutated

identity = pd.DataFrame({"genome": ["G", "G"], "family": ["B-lectin", "GH5_2"],
                         "accepted_id": ["in", "m"], "order_id": ["", "m"]})
with tempfile.TemporaryDirectory() as tmp:
    path = Path(tmp) / "d.tsv"
    row = "G\tB-lectin\tin\tresolved\tnative\ts1\t300\t600\n"
    path.write_text("genome\tfamily\taccepted_id\tstatus\tframe\torder_scaffold\torder_start\torder_end\n" + row)
    assert [x["accepted_id"] for x in modelless_loci(path, identity, chd)] == ["in"]
    for bad in (row + "G\tGH5_2\tm\tresolved\tnative\ts1\t1\t2\n",      # copy with an order identity
                row.replace("s1", "s9"),                                 # outside gene-order frame
                row.replace("native\ts1\t300\t600", "no_order_scaffold\ts1\t300\t600")):
        path.write_text("genome\tfamily\taccepted_id\tstatus\tframe\torder_scaffold\torder_start\torder_end\n" + bad)
        try:
            modelless_loci(path, identity, chd)
        except ValueError:
            continue
        raise AssertionError(bad)
print("MODELLESS_SCAN_TEST_PASS")
