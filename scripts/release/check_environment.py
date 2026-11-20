#!/usr/bin/env python3
"""Fail the build unless the Python and R environments match environment/*-lock.txt.

Each lock line is `package==version`; the first line names the interpreter
(`python==X.Y.Z`, `R==X.Y.Z`). Run with the build's Python; --rscript names the
build's Rscript. `--write` rewrites both locks from the running environments.

Environment: project Python 3 (standard library only).
"""
import argparse
import importlib.metadata
import platform
import subprocess
import sys
from pathlib import Path

ENV = Path(__file__).resolve().parents[2] / "environment"


def read_lock(path):
    return dict(line.strip().split("==", 1) for line in path.read_text().splitlines() if line.strip())


def python_versions(packages):
    found = {"python": platform.python_version()}
    for name in packages:
        try:
            found[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            found[name] = "MISSING"
    return found


def r_versions(rscript, packages):
    code = ('cat("R", as.character(getRversion()), "\\n"); for (p in commandArgs(TRUE)) '
            'cat(p, tryCatch(as.character(packageVersion(p)), error=function(e) "MISSING"), "\\n")')
    out = subprocess.run([rscript, "-e", code, *packages], check=True, capture_output=True, text=True).stdout
    return dict(line.split() for line in out.splitlines() if line.strip())


parser = argparse.ArgumentParser()
parser.add_argument("--rscript", required=True)
parser.add_argument("--write", action="store_true")
args = parser.parse_args()
py_lock, r_lock = ENV / "python-lock.txt", ENV / "r-lock.txt"
py_want, r_want = read_lock(py_lock), read_lock(r_lock)
py_have = python_versions([p for p in py_want if p != "python"])
r_have = r_versions(args.rscript, [p for p in r_want if p != "R"])
if args.write:
    py_lock.write_text("".join(f"{k}=={v}\n" for k, v in py_have.items()))
    r_lock.write_text("".join(f"{k}=={v}\n" for k, v in r_have.items()))
    print("ENVIRONMENT_LOCKS_WRITTEN")
    sys.exit(0)
drift = [f"{k}: lock {v}, found {have.get(k)}" for want, have in ((py_want, py_have), (r_want, r_have))
         for k, v in want.items() if have.get(k) != v]
for line in drift:
    print(f"FAIL\tenvironment {line}")
if drift:
    raise SystemExit(1)
print(f"ENVIRONMENT_PASS python={py_have['python']} R={r_have['R']} packages={len(py_want) + len(r_want) - 2}")
