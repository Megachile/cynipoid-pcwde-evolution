#!/usr/bin/env python3
"""Compare a new supplement Google Doc export with the pinned export, paragraph by paragraph.

Reports only; it never edits the template. For each section it lists paragraphs that were
added, removed or changed (word-level), and every numeral in new or changed text that the
pinned paragraph did not contain: those are the numbers to tokenise (or to list in
config/supplement_literals.tsv if they are parameters, versions, dates or citations)
when the change is ported into authority/supplement/supplement_source.md.

Usage: doc_edit_report.py --repo . --new <export.md> [--pinned <export.md>] [--out report.md]

Environment: project Python 3 (standard library only).
"""
import argparse
import csv
import difflib
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "release"))
import supplement_template as st  # noqa: E402
from author_export import supplement_export


def sections(path):
    text = path.read_text(encoding="utf-8")
    if text.startswith("<!--"):  # Drive sync header (may hold the private document URL)
        text = text.split("\n", 1)[1]
    out, name = {}, None
    for block in re.split(r"\n\s*\n", text):
        block = block.strip()
        if not block:
            continue
        if block.startswith("## "):
            name = block[3:].split("\n")[0].strip()
            out[name] = []
            block = block.split("\n", 1)[1].strip() if "\n" in block else ""
            if not block:
                continue
        out.setdefault(name, []).append(block)
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--repo", type=Path, required=True)
    p.add_argument("--new", type=Path, required=True)
    p.add_argument("--pinned", type=Path)
    p.add_argument("--out", type=Path)
    a = p.parse_args()
    pinned = a.pinned
    if pinned is None:
        pinned = supplement_export(a.repo)
    old, new = sections(pinned), sections(a.new)
    lines = [f"# Supplement export changes: {a.new.name} vs pinned {pinned.name}", ""]
    changed = 0
    for name in list(old) + [n for n in new if n not in old]:
        if name not in new:
            lines += [f"## {name}", "", "Section removed in the new export.", ""]
            changed += 1
            continue
        if name not in old:
            lines += [f"## {name}", "", "New section."]
        a_par, b_par = old.get(name, []), new[name]
        report = []
        for op, i1, i2, j1, j2 in difflib.SequenceMatcher(None, a_par, b_par, autojunk=False).get_opcodes():
            if op == "equal":
                continue
            before, after = a_par[i1:i2], b_par[j1:j2]
            for k in range(max(len(before), len(after))):
                x = before[k] if k < len(before) else ""
                y = after[k] if k < len(after) else ""
                words_x, words_y = x.split(), y.split()
                edits = []
                for wop, w1, w2, v1, v2 in difflib.SequenceMatcher(None, words_x, words_y, autojunk=False).get_opcodes():
                    if wop != "equal":
                        edits.append(f"  - `{' '.join(words_x[w1:w2]) or '∅'}` → `{' '.join(words_y[v1:v2]) or '∅'}`")
                kind = "added" if not x else "removed" if not y else "changed"
                old_numerals = [n for n, _, _ in st.numerals(x)]
                fresh = [n for n, _, _ in st.numerals(y)]
                for n in old_numerals:
                    if n in fresh:
                        fresh.remove(n)
                report.append(f"- Paragraph {kind}: “{(y or x)[:100]}…”")
                report += edits[:40] if kind == "changed" else []
                if fresh:
                    report.append("  - New numerals to tokenise or list as literals: " + ", ".join(fresh))
        if report:
            changed += 1
            if name in old:
                lines += [f"## {name}", ""]
            lines += report + [""]
    lines.append(f"{changed} section(s) with changes. This report does not apply edits.")
    text = "\n".join(lines) + "\n"
    if a.out:
        a.out.write_text(text, encoding="utf-8")
    sys.stdout.write(text)


if __name__ == "__main__":
    main()
