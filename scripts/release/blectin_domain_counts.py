"""Count the rendered B-lectin domain tree through its pinned domain manifest.

Method: read every terminal of the rendered rooted Newick, require a one-to-one
join to the manifest `tip_id` column, classify tips by the manifest `tip_class`
(focal or bacterial) and count unique focal `biological_unit_id` values as
proteins. The SHA256 of both inputs is recorded with the counts. Run once as a
script to write the build intermediate; the Table S3 caption, numbers.json and
the Figure S6 caption all consume that file rather than recounting.
"""
import argparse
import csv
import hashlib
import json
import re
from pathlib import Path


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def counts(tree, manifest):
    from Bio import Phylo

    tips = [tip.name for tip in Phylo.read(tree, "newick").get_terminals()]
    if not tips or any(not tip for tip in tips) or len(tips) != len(set(tips)):
        raise ValueError("empty, unnamed or duplicate rendered tree tips")
    with Path(manifest).open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    ids = [row["tip_id"] for row in rows]
    if len(ids) != len(set(ids)) or set(ids) != set(tips):
        raise ValueError("rendered tree and domain manifest lack one-to-one tip coverage")
    if {row["tip_class"] for row in rows} != {"focal", "bacterial"}:
        raise ValueError("domain manifest must classify focal and bacterial tips")
    if any(not row["biological_unit_id"] for row in rows):
        raise ValueError("domain manifest has empty biological unit IDs")
    focal = [row for row in rows if row["tip_class"] == "focal"]
    bacterial = [row for row in rows if row["tip_class"] == "bacterial"]
    if len({row["biological_unit_id"] for row in bacterial}) != 1:
        raise ValueError("bacterial tips do not represent one outgroup protein")
    return {
        "total_domains": len(tips),
        "cynipoid_domains": len(focal),
        "cynipoid_proteins": len({row["biological_unit_id"] for row in focal}),
        "bacterial_domains": len(bacterial),
        "tree_sha256": sha256(tree),
        "manifest_sha256": sha256(manifest),
    }


def verify_counts(tree, manifest, recorded):
    observed = counts(tree, manifest)
    if observed != recorded:
        raise ValueError("Figure S6 caption counts disagree with rendered tree or manifest")
    return observed


def load_counts(path):
    return json.loads(Path(path).read_text())


def table_s3_tree_clause(recorded, figure_label):
    """Table S3 sentence naming the rendered-tree subset of the 69 qualifying domains."""
    return f"{recorded['total_domains']} of these are Figure {figure_label} tips"


def verify_table_s3_caption(caption, recorded, figure_label):
    claims = re.findall(rf"(\d+)(?: of these are)? Figure {re.escape(figure_label)} tips", caption)
    if claims != [str(recorded["total_domains"])]:
        raise ValueError(f"Table S3 Figure {figure_label} tip count disagrees with rendered tree: {claims}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--tree", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(counts(args.tree, args.manifest), indent=2, sort_keys=True) + "\n")
