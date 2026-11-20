#!/usr/bin/env python3
"""Stage only built intron tables needed by the Figure 3 renderer."""
from __future__ import annotations

import argparse
import shutil
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--observations", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--background", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    files = [(args.observations / "ref_proteins.tsv", "ref_proteins.tsv"),
             (args.background / "background_intron_details.tsv", "baseline_intron_details.tsv")]
    for family in ("GH5_2", "PL1_8", "PL4_4"):
        files.extend((args.summary / f"{family}_{suffix}.tsv", f"{family}_{suffix}.tsv")
                     for suffix in ("groups", "context_summary", "context_rich", "intron_positions"))
    for source, name in files:
        if not source.is_file() or source.stat().st_size == 0:
            raise FileNotFoundError(source)
        shutil.copy2(source, args.out / name)
    print(f"Figure 3 inputs: {len(files)} built tables")


if __name__ == "__main__":
    main()
