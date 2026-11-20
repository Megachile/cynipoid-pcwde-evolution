#!/usr/bin/env python3
"""Exercise the Figure S6 and Table S3 count gates, including mismatch fixtures."""
import argparse
import csv
import tempfile
from pathlib import Path

from blectin_domain_counts import counts, load_counts, verify_counts, verify_table_s3_caption

p = argparse.ArgumentParser()
p.add_argument("--tree", type=Path, required=True)
p.add_argument("--manifest", type=Path, required=True)
p.add_argument("--counts", type=Path, required=True)
p.add_argument("--table-s3-caption", type=Path, required=True)
p.add_argument("--figure-numbers", type=Path, required=True)
a = p.parse_args()
observed = counts(a.tree, a.manifest)
verify_counts(a.tree, a.manifest, load_counts(a.counts))
with a.figure_numbers.open(newline="") as handle:
    label = {r["figure_id"]: r["label"] for r in csv.DictReader(handle, delimiter="\t")}["b_lectin_domain_tree"]
caption = a.table_s3_caption.read_text(encoding="utf-8")
verify_table_s3_caption(caption, observed, label)
# Mismatch fixtures: the former fixed tip count, and the correct sentence with an altered count.
total = str(observed["total_domains"])
for bad_caption in (caption + f" giving {observed['total_domains'] + 1} Figure {label} tips.",
                    caption.replace(f"{total} of these are Figure {label} tips",
                                    f"{observed['total_domains'] + 1} of these are Figure {label} tips")):
    try:
        verify_table_s3_caption(bad_caption, observed, label)
    except ValueError as error:
        if "tip count disagrees" not in str(error):
            raise
    else:
        raise SystemExit("negative Table S3 caption fixture did not fail")
wrong = dict(observed, cynipoid_domains=observed["cynipoid_domains"] + 1)
try:
    verify_counts(a.tree, a.manifest, wrong)
except ValueError as error:
    if "caption counts disagree" not in str(error):
        raise
else:
    raise SystemExit("negative count fixture did not fail")
with tempfile.TemporaryDirectory() as tmp:
    bad = Path(tmp) / "bad.tsv"
    with a.manifest.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    rows[0]["tip_id"] += "_mismatch"
    with bad.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    try:
        counts(a.tree, bad)
    except ValueError as error:
        if "one-to-one tip coverage" not in str(error):
            raise
    else:
        raise SystemExit("negative manifest fixture did not fail")
print("BLECTIN_COUNT_GATE_PASS", {k: observed[k] for k in
      ("total_domains", "cynipoid_domains", "cynipoid_proteins", "bacterial_domains")})
