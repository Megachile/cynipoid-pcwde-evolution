#!/usr/bin/env python3
"""Self-test of check_closure.py's deposit rules on synthetic files; run by make before the products.

A clean package passes; an unrelated extra file, an altered copy of a mapped
input, a generated file left over from an earlier build (its receipt carries
another build's token, timestamps irrelevant), and a generated file changed
after its producer recorded it each fail.

Environment: project Python 3 (standard library only).
"""
import contextlib
import io
import shutil
import time
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import check_closure as cc  # noqa: E402

LAYOUT = ("kind\tdeposit_path\tsource\treason\n"
          "generated\tpkg/README.md\tscripts/release/check_closure.py\tsynthetic generated file\n"
          "copy\tpkg/inputs/{rest}\t{rest}\tsynthetic verbatim copies\n")


def scenario(tmp, mutate):
    data, build = tmp / "data", tmp / "build"
    (data / "repository_inputs").mkdir(parents=True)
    (data / "repository_inputs/members.tsv").write_text("id\nA\n")
    out = build / "closure"
    out.mkdir(parents=True)
    (out / "reset.done").write_text("1\n")
    (out / "build_token").write_text("current-build\n")
    (out / "inputs.tsv").write_text("relative_path\tstatus\tdeclared_class\nmembers.tsv\tread\t\n")
    pkg = build / "data_deposit/pkg"
    (pkg / "inputs").mkdir(parents=True)
    shutil.copyfile(data / "repository_inputs/members.tsv", pkg / "inputs/members.tsv")
    (pkg / "README.md").write_text("readme\n")
    cc.write_receipt(out, "scripts/release/check_closure.py", build / "data_deposit")
    layout = tmp / "layout.tsv"
    layout.write_text(LAYOUT)
    mutate(pkg, out)
    log = io.StringIO()
    try:
        with contextlib.redirect_stdout(log):
            cc.check_deposit(data, build, out, layout=layout, repo=cc.REPO)
    except SystemExit:
        reasons = [line.split("\t", 1)[1] for line in log.getvalue().splitlines() if line.startswith("FAIL\t")]
        return "fail", reasons
    return "pass", []


def unrelated(pkg, out):
    (pkg / "inputs/extra.tsv").write_text("not from any input\n")
    (pkg / "notes.txt").write_text("stray\n")


def altered(pkg, out):
    (pkg / "inputs/members.tsv").write_text("id\nA\nB\n")


def stale(pkg, out):
    # Left over from an earlier build: the only receipt for it carries another build's token.
    (out / "build_token").write_text("next-build\n")


def generated_altered(pkg, out):
    # Changed after the producer recorded it (e.g. old content with refreshed metadata).
    (pkg / "README.md").write_text("old readme\n")


def main():
    expected = {"clean": ("pass", lambda pkg, out: None), "unrelated_file": ("fail", unrelated),
                "altered_copy": ("fail", altered), "stale_generated": ("fail", stale),
                "generated_altered": ("fail", generated_altered),
                "unread_declared": ("fail", lambda pkg,out: (out/"inputs.tsv").write_text("relative_path\tstatus\tdeclared_class\nmembers.tsv\tfrozen_member_layer\tfrozen_member_layer\n"))}
    results = {}
    for name, (want, mutate) in expected.items():
        with tempfile.TemporaryDirectory() as tmp:
            results[name], reasons = scenario(Path(tmp), mutate)
        for reason in reasons:
            print(f"SELFTEST_EXPECTED\t{name}\t{reason}")
    bad = {k: v for k, v in results.items() if v != expected[k][0]}
    print("CLOSURE_SELFTEST " + " ".join(f"{k}={v}" for k, v in results.items()))
    if bad:
        raise SystemExit(f"closure self-test: unexpected results {bad}")
    print("CLOSURE_SELFTEST_PASS")


if __name__ == "__main__":
    main()
