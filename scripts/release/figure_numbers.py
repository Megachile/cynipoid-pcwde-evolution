#!/usr/bin/env python3
"""The one place figure numbers live: config/figure_numbers.tsv.

Figures are produced and referenced by stable id (tree_counts, origin_trees, ...).
A number (label) is attached only at assembly: the numbered export, the
manuscript and the supplement read it from here. `source` is `build` (path under
build/figures) or `static` (path under the data root; pinned in config/build_inputs.tsv).

Run as a script to write build/numbered_figures/Figure_<label>.pdf (and .png
where the renderer wrote one) for every mapped figure.

Environment: project Python 3 (standard library only).
"""
import argparse
import csv
import shutil
from pathlib import Path

MAP = Path(__file__).resolve().parents[2] / "config/figure_numbers.tsv"


def load():
    with MAP.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    labels = [r["label"] for r in rows]
    if len(set(labels)) != len(labels) or len({r["figure_id"] for r in rows}) != len(rows):
        raise SystemExit("config/figure_numbers.tsv repeats a label or figure id")
    return rows


def source_path(row, build, data_root):
    if row["source"] == "reserved":
        raise ValueError(f"reserved figure has no source: {row['figure_id']}")
    return (build / "figures" / row["path"]) if row["source"] == "build" else (data_root / row["path"])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--build", type=Path, required=True)
    parser.add_argument("--data-root", type=Path, required=True)
    args = parser.parse_args()
    out = args.build / "numbered_figures"
    if out.exists():
        raise SystemExit(f"refusing existing output: {out}")
    out.mkdir(parents=True)
    with (out / "MAP.tsv").open("w", newline="") as handle:
        w = csv.writer(handle, delimiter="\t", lineterminator="\n")
        w.writerow(["file", "figure_id", "source"])
        for row in load():
            if row["source"] == "reserved":
                continue
            src = source_path(row, args.build, args.data_root)
            if not src.is_file():
                raise SystemExit(f"missing figure {row['figure_id']}: {src}")
            for suffix in (".pdf", ".png"):
                path = src.with_suffix(suffix)
                if path.is_file():
                    name = f"Figure_{row['label']}{suffix}"
                    shutil.copyfile(path, out / name)
                    w.writerow([name, row["figure_id"], row["path"]])
    print(f"NUMBERED_FIGURES_PASS figures={len(load())} out={out}")


if __name__ == "__main__":
    main()
