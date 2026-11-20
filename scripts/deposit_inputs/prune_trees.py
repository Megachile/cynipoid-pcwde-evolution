#!/usr/bin/env python3
"""Prune exact tips from Newick/NHX trees without re-inference.

Every untouched token (labels, support, NHX comments, branch-length strings)
is copied verbatim.  When a removed tip leaves a unary internal node, that node
is collapsed and its edge length is added exactly (Decimal) to the surviving
child's edge; the child keeps its own label/support/comment, as ape::drop.tip
does.  Input must round-trip byte-identically before any pruning is applied.

Usage: prune_trees.py PLAN.tsv DATA_ROOT OUT_EVIDENCE.tsv
PLAN.tsv columns: file (relative to DATA_ROOT), tip (exact terminal label).
Files are rewritten in place; the original is copied to <evidence dir>/before/.
"""
from __future__ import annotations

import csv
import shutil
import sys
from collections import defaultdict
from decimal import Decimal
from pathlib import Path


class Node:
    __slots__ = ("children", "label", "length", "comment", "parent")

    def __init__(self):
        self.children: list[Node] = []
        self.label = ""
        self.length: str | None = None
        self.comment = ""
        self.parent: Node | None = None


def parse(text: str) -> tuple[Node, str]:
    body = text.rstrip("\n")
    trailer = text[len(body):]
    if not body.endswith(";") or "'" in body or '"' in body:
        raise SystemExit("unsupported Newick framing or quoting")
    s = body[:-1]
    i = 0

    def read_attrs(node: Node) -> None:
        nonlocal i
        start = i
        while i < len(s) and s[i] not in "():,;[":
            i += 1
        node.label = s[start:i]
        if i < len(s) and s[i] == ":":
            i += 1
            start = i
            while i < len(s) and s[i] not in "():,;[":
                i += 1
            node.length = s[start:i]
        if i < len(s) and s[i] == "[":
            end = s.index("]", i)
            node.comment = s[i:end + 1]
            i = end + 1
        # IQ-TREE/NHX writers put the comment after the length.

    def subtree() -> Node:
        nonlocal i
        node = Node()
        if s[i] == "(":
            i += 1
            while True:
                child = subtree()
                child.parent = node
                node.children.append(child)
                if s[i] == ",":
                    i += 1
                elif s[i] == ")":
                    i += 1
                    break
                else:
                    raise SystemExit(f"unexpected character at {i}: {s[i:i+20]!r}")
        read_attrs(node)
        return node

    if s[0] != "(":
        raise SystemExit("tree does not start with '('")
    tree = subtree()
    if i != len(s):
        raise SystemExit(f"trailing content at {i}: {s[i:i+20]!r}")
    return tree, trailer


def serialize(node: Node) -> str:
    out: list[str] = []

    def emit(n: Node) -> None:
        if n.children:
            out.append("(")
            for k, c in enumerate(n.children):
                if k:
                    out.append(",")
                emit(c)
            out.append(")")
        out.append(n.label)
        if n.length is not None:
            out.append(":" + n.length)
        out.append(n.comment)

    emit(node)
    return "".join(out) + ";"


def leaves(node: Node) -> list[Node]:
    if not node.children:
        return [node]
    return [x for c in node.children for x in leaves(c)]


def add_lengths(a: str | None, b: str | None) -> str | None:
    if a is None and b is None:
        return None
    total = Decimal(a or "0") + Decimal(b or "0")
    return format(total.normalize() if total == 0 else total, "f")


def prune(root: Node, tip: str, log: list[dict]) -> Node:
    hits = [x for x in leaves(root) if x.label == tip]
    if len(hits) != 1:
        raise SystemExit(f"expected exactly one exact tip {tip!r}, found {len(hits)}")
    leaf = hits[0]
    parent = leaf.parent
    parent.children.remove(leaf)
    event = {"tip": tip, "removed_tip_length": leaf.length or "",
             "collapsed_node_label": "", "collapsed_node_length": "",
             "surviving_child_old_length": "", "surviving_child_new_length": ""}
    if len(parent.children) == 1:
        child = parent.children[0]
        grand = parent.parent
        if grand is None:
            raise SystemExit(f"removal of {tip!r} would collapse the root; refusing")
        event.update(collapsed_node_label=parent.label + parent.comment,
                     collapsed_node_length=parent.length or "",
                     surviving_child_old_length=child.length or "")
        child.length = add_lengths(parent.length, child.length)
        event["surviving_child_new_length"] = child.length or ""
        child.parent = grand
        grand.children[grand.children.index(parent)] = child
    elif not parent.children:
        raise SystemExit("degenerate parent")
    log.append(event)
    return root


def main() -> None:
    sys.setrecursionlimit(100000)
    plan_path, data_root, evidence_path = map(Path, sys.argv[1:4])
    by_file: dict[str, list[str]] = defaultdict(list)
    with plan_path.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            by_file[row["file"]].append(row["tip"])
    before_dir = evidence_path.parent / "before"
    rows = []
    for rel, tips in sorted(by_file.items()):
        path = data_root / rel
        text = path.read_text(encoding="utf-8")
        root, trailer = parse(text)
        if serialize(root) + trailer != text:
            raise SystemExit(f"round-trip failure before pruning: {rel}")
        n_before = len(leaves(root))
        log: list[dict] = []
        for tip in tips:
            root = prune(root, tip, log)
        new_text = serialize(root) + trailer
        reparsed, _ = parse(new_text)
        if serialize(reparsed) + trailer != new_text:
            raise SystemExit(f"round-trip failure after pruning: {rel}")
        n_after = len(leaves(root))
        if n_after != n_before - len(tips):
            raise SystemExit(f"tip count mismatch: {rel}")
        saved = before_dir / rel
        if not saved.exists():
            saved.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(path, saved)
        path.write_text(new_text, encoding="utf-8")
        for event in log:
            rows.append({"file": rel, "tips_before": n_before, "tips_after": n_after, **event})
        print(f"{rel}\t{n_before}->{n_after}")
    with evidence_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
