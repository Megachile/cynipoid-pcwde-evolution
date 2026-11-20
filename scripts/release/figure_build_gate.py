#!/usr/bin/env python3
"""Require a clean checkout and a confined, freshly rendered figure set."""

import argparse
import hashlib
import json
import sys
import subprocess
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import figure_numbers  # noqa: E402


def git(repo, *args):
    return subprocess.check_output(["git", "-C", str(repo), *args], text=True).strip()


def inside(path, parent):
    return path.resolve().is_relative_to(parent.resolve())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--build", type=Path, required=True)
    parser.add_argument("--receipt", action="store_true")
    args = parser.parse_args()
    repo, data, build = (p.resolve() for p in (args.repo, args.data_root, args.build))
    if git(repo, "status", "--porcelain", "--untracked-files=all"):
        raise SystemExit("figure build requires a clean checkout")
    if not (inside(data, repo) or data == repo):
        # The independent, unpacked data root is a declared external input.
        if not (data / "repository_inputs" / "figures").is_dir():
            raise SystemExit("figure data root lacks declared repository inputs")
    if not inside(build, repo) and not inside(build, repo.parent):
        raise SystemExit("build output is outside the checkout or its run directory")
    commit = git(repo, "rev-parse", "HEAD")
    marker = build / "figure_render_intent.json"
    if args.receipt:
        if not marker.is_file() or json.loads(marker.read_text()).get("commit") != commit:
            raise SystemExit("missing matching fresh-render intent")
        figures = build / "figures"
        expected = {r["path"] for r in figure_numbers.load() if r["source"] == "build"}
        missing = sorted(p for p in expected if not (figures / p).is_file())
        if missing:
            raise SystemExit(f"missing fresh figures: {missing}")
        s5 = data / "build_inputs/euceroptres_genome/euceroptres_genome.pdf"
        if s5.is_symlink() or not s5.is_file():
            raise SystemExit("missing pinned Figure S5 data-root asset")
        s5_hash = hashlib.sha256(s5.read_bytes()).hexdigest()
        if s5_hash != "c968c62fb81a2da159945ffe4286170174c264ef759ae00a8328f82cd77eb78f":
            raise SystemExit("Figure S5 checksum mismatch")
        hashes = {p: hashlib.sha256((figures / p).read_bytes()).hexdigest() for p in sorted(expected)}
        hashes["data_root:build_inputs/euceroptres_genome/euceroptres_genome.pdf"] = s5_hash
        (build / "figure_build_receipt.json").write_text(json.dumps({
            "commit": commit, "data_root": str(data),
            "figures": hashes,
        }, indent=2) + "\n")
    else:
        if build.exists():
            raise SystemExit("refusing existing build output")
        build.mkdir(parents=True)
        marker.write_text(json.dumps({"commit": commit, "data_root": str(data)}) + "\n")


if __name__ == "__main__":
    main()
