#!/usr/bin/env python3
"""Build closure: every deposit input is read by the build or declared; every reader-deposit file is a mapped copy or a declared build product.

Three modes, run by `make manuscript`:

  --reset    right after `verify` (which hashes every input). Sets each
             repository_inputs/ file's access time to its modification time,
             then proves on a canary file in the data root that reading a file
             advances its access time on this filesystem. If it does not (e.g.
             a noatime mount), the build fails: the check cannot run. Writes
             BUILD/closure/reset.done, whose time marks the start of the build.
  --inputs   after the products are built, before `data-deposit` (which
             copies every input). A file was read by the build iff its access
             time is now later than its modification time. Every unread file
             must match a row of config/closure_declared.tsv (first matching
             pattern wins); rows of class known_gap are reported whether or not
             the file was read. Reads and classes go to BUILD/closure/inputs.tsv.
  --deposit  after `data-deposit`. Every file under BUILD/data_deposit must be
             explained by config/deposit_layout.tsv, path by path:
               copy       the file's bytes equal the input at the mapped path
                          (repository_inputs/<source>), and that input was read
                          or is declared; copies of known-gap inputs are listed;
               generated  the path is declared as written by a named committed
                          producer, the producer exists in the repository, and
                          the file was created or written during this build: its
                          status-change time (ctime) is later than
                          BUILD/closure/reset.done. ctime is used because the
                          deposit builder normalizes modification times for a
                          byte-reproducible archive, and no user call can set
                          ctime back.
             A file matching no row, a mapped copy whose bytes differ from its
             source, and a generated file older than the build all fail.

scripts/release/test_check_closure.py exercises these rules on synthetic
files (an unrelated file, an altered copy and a stale generated file must
fail) and runs before the products are built.

Prints CLOSURE_INPUTS_PASS / CLOSURE_DEPOSIT_PASS with counts per class and
lists every known gap, so open gaps stay visible in every build log.

Method limits: access time records that a file was opened and read, not that
its contents reached a product; a step that reads a file and ignores it counts
as a read. Files read only by `verify` (hashing) count as unread. Reads by
later steps (data-deposit, stage-check) are not counted. Generated deposit
files are tied to their producer by declared path and write time, not by
content: an altered generated file written during the build would pass.

Environment: project Python 3 (standard library only).
"""
import argparse
import csv
import fnmatch
import hashlib
import json
import os
import secrets
import sys
from collections import Counter
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
DECLARED = REPO / "config/closure_declared.tsv"
LAYOUT = REPO / "config/deposit_layout.tsv"
CLASSES = {"frozen_member_layer", "frozen_cached_tree", "frozen_out_of_family", "static_build_input",
           "display_config", "historical_evidence", "known_gap"}


def inputs(data_root):
    root = data_root / "repository_inputs"
    return root, sorted(p for p in root.rglob("*") if p.is_file())


def read_tsv(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def declarations(path=DECLARED):
    rows = read_tsv(path)
    bad = [r["pattern"] for r in rows if r["class"] not in CLASSES or not r["reason"].strip()]
    if bad:
        raise SystemExit(f"closure declarations with an unknown class or no reason: {bad}")
    return rows


def declared_class(rel, rows):
    for r in rows:
        if fnmatch.fnmatchcase(rel, r["pattern"]):
            return r
    return None


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def reset(data_root, out):
    root, files = inputs(data_root)
    for p in files:
        st = p.stat()
        os.utime(p, ns=(st.st_mtime_ns, st.st_mtime_ns))
    canary = data_root / ".closure_canary"
    canary.write_text("canary\n")
    st = canary.stat()
    os.utime(canary, ns=(st.st_mtime_ns - 10**9, st.st_mtime_ns - 10**9))
    canary.read_text()
    tracked = canary.stat().st_atime_ns > canary.stat().st_mtime_ns
    canary.unlink()
    if not tracked:
        raise SystemExit("closure check cannot run: reading a file does not advance its access time on this filesystem")
    out.mkdir(parents=True, exist_ok=True)
    if (out / "build_token").exists():
        raise SystemExit("closure reset refuses an existing build token; start from a fresh build directory")
    (out / "build_token").write_text(secrets.token_hex(16) + "\n")
    (out / "reset.done").write_text(f"{len(files)}\n")
    print(f"CLOSURE_RESET files={len(files)} atime_tracking=yes")


def check_inputs(data_root, out, declared=DECLARED):
    if not (out / "reset.done").is_file():
        raise SystemExit("closure --inputs needs the --reset step first")
    root, files = inputs(data_root)
    rows = declarations(declared)
    counts, gaps, undeclared, unused_declared = Counter(), [], [], []
    records = []
    for p in files:
        rel = p.relative_to(root).as_posix()
        st = p.stat()
        read = st.st_atime_ns > st.st_mtime_ns
        decl = declared_class(rel, rows)
        if decl and decl["class"] == "known_gap":
            # Known gaps stay visible whether or not the build reads the file.
            counts["known_gap_read" if read else "known_gap"] += 1
            status = "known_gap"
            gaps.append(f"{rel}\t{'read' if read else 'unread'}\t{decl['reason']}")
        elif read:
            counts["read_by_build"] += 1
            status = "read"
        elif decl:
            undeclared.append(rel)
            status = "UNREAD_DECLARED"
        else:
            undeclared.append(rel)
            status = "UNDECLARED"
        records.append((rel, status, decl["class"] if decl else ""))
    for r in rows:
        if not any(fnmatch.fnmatchcase(rel, r["pattern"]) for rel, _, _ in records):
            unused_declared.append(r["pattern"])
    with (out / "inputs.tsv").open("w", newline="") as handle:
        w = csv.writer(handle, delimiter="\t", lineterminator="\n")
        w.writerow(["relative_path", "status", "declared_class"])
        w.writerows(records)
    for g in gaps:
        print(f"KNOWN_GAP\t{g}")
    for pattern in unused_declared:
        print(f"FAIL\tdeclaration matches no deposit input: {pattern}")
    for rel in undeclared:
        print(f"FAIL\tinput neither read by the build nor declared: {rel}")
    if undeclared or unused_declared:
        raise SystemExit(1)
    summary = " ".join(f"{k}={v}" for k, v in sorted(counts.items()))
    print(f"CLOSURE_INPUTS_PASS files={len(files)} {summary} known_gaps={len(gaps)}")


def write_receipt(closure_dir, producer, package_root):
    """Called by a producer when it finishes: record this build's token and the hash of every file it wrote."""
    token_file = closure_dir / "build_token"
    if not token_file.is_file():
        return None  # producer run outside a closure-checked build: no receipt, so its outputs cannot pass
    outputs = {p.relative_to(package_root).as_posix(): sha(p) for p in sorted(package_root.rglob("*")) if p.is_file()}
    receipt = closure_dir / f"receipt_{Path(producer).stem}.json"
    receipt.write_text(json.dumps({"token": token_file.read_text().strip(), "producer": producer, "outputs": outputs},
                                  indent=1, sort_keys=True) + "\n")
    return receipt


def load_receipts(out):
    """Generated files this build's producers recorded: {deposit path: (producer, sha256)}."""
    token = (out / "build_token").read_text().strip() if (out / "build_token").is_file() else None
    recorded = {}
    for path in sorted(out.glob("receipt_*.json")):
        receipt = json.loads(path.read_text())
        if token is None or receipt.get("token") != token:
            continue  # a receipt from another build proves nothing about this one
        for rel, digest in receipt["outputs"].items():
            recorded[rel] = (receipt["producer"], digest)
    return recorded


def layout_match(rel, pattern):
    """Return the {rest} part (or '' for an exact path) if rel matches the layout pattern, else None."""
    if pattern.endswith("{rest}"):
        prefix = pattern[: -len("{rest}")]
        return rel[len(prefix):] if rel.startswith(prefix) and len(rel) > len(prefix) else None
    return "" if rel == pattern else None


def check_deposit(data_root, build, out, layout=LAYOUT, repo=REPO):
    root, _ = inputs(data_root)
    status = {r["relative_path"]: r["status"] for r in read_tsv(out / "inputs.tsv")}
    rules = read_tsv(layout)
    bad = [r for r in rules if r["kind"] not in ("copy", "generated") or not r["reason"].strip()]
    if bad:
        raise SystemExit(f"deposit layout rows with an unknown kind or no reason: {bad}")
    recorded = load_receipts(out)
    deposit = build / "data_deposit"
    counts, gap_copies, failures, used = Counter(), [], [], set()
    for p in sorted(q for q in deposit.rglob("*") if q.is_file()):
        rel = p.relative_to(deposit).as_posix()
        explained, mismatched = False, []
        for i, rule in enumerate(rules):
            rest = layout_match(rel, rule["deposit_path"])
            if rest is None:
                continue
            if rule["kind"] == "generated":
                producer, digest = recorded.get(rel, (None, None))
                if not (repo / rule["source"]).is_file():
                    failures.append(f"{rel}: declared producer missing from the repository: {rule['source']}")
                elif producer is None:
                    failures.append(f"{rel}: declared as generated by {rule['source']}, but no receipt from this build lists it")
                elif producer != rule["source"]:
                    failures.append(f"{rel}: this build's receipt names {producer}, not the declared producer {rule['source']}")
                elif digest != sha(p):
                    failures.append(f"{rel}: bytes differ from what {producer} recorded writing in this build")
                else:
                    counts["generated_by_declared_producer"] += 1
                explained = True
                used.add(i)
                break
            source = rule["source"].replace("{rest}", rest)
            src = root / source
            if not src.is_file():
                continue
            if sha(src) != sha(p):
                mismatched.append(source)
                continue
            kind = status.get(source, "UNDECLARED")
            if kind == "UNDECLARED":
                failures.append(f"{rel}: copies input {source}, which was neither read nor declared")
            elif kind != "read":
                failures.append(f"{rel}: public input {source} has no committed build read ({kind})")
            else:
                counts["copy_of_read_or_declared_input"] += 1
            explained = True
            used.add(i)
            break
        if not explained:
            if mismatched:
                failures.append(f"{rel}: bytes differ from its mapped source {', '.join(mismatched)}")
            else:
                failures.append(f"{rel}: no deposit-layout rule explains this file")
    for i, rule in enumerate(rules):
        if i not in used:
            failures.append(f"deposit layout row matches no deposit file: {rule['kind']} {rule['deposit_path']}")
    for g in gap_copies:
        print(f"KNOWN_GAP_IN_READER_DEPOSIT\t{g}")
    for f in failures:
        print(f"FAIL\t{f}")
    if failures:
        raise SystemExit(1)
    summary = " ".join(f"{k}={v}" for k, v in sorted(counts.items()))
    print(f"CLOSURE_DEPOSIT_PASS files={sum(counts.values())} {summary}")


def main():
    ap = argparse.ArgumentParser()
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--reset", action="store_true")
    mode.add_argument("--inputs", action="store_true")
    mode.add_argument("--deposit", action="store_true")
    ap.add_argument("--data-root", type=Path, required=True)
    ap.add_argument("--build", type=Path, required=True)
    a = ap.parse_args()
    out = a.build / "closure"
    if a.reset:
        reset(a.data_root, out)
    elif a.inputs:
        check_inputs(a.data_root, out)
    else:
        check_deposit(a.data_root, a.build, out)


if __name__ == "__main__":
    sys.exit(main())
