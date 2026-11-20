#!/usr/bin/env python3
"""Main-text supplement pointers: legacy draft label -> stable ID -> final label, plus a failing gate.

The pinned author draft (authority/manuscript/manuscript_draft_20261109T042633Z.md) cites the
supplement with the numbering in use when it was written. Each legacy label is resolved by its
subject to a stable ID here; the final label is then read from the committed number maps
(config/table_numbers.tsv, config/figure_numbers.tsv, config/supplement_sections.tsv), so a
pointer always follows its subject. All swaps happen in one pass so that a newly assigned
label is never remapped again.

`check(prose)` fails when any supplement pointer left in the main text does not resolve to the
final map (retired items, legacy section numerals, the unsplit Table S2).

Environment: project Python 3 (standard library only).
"""
from __future__ import annotations

import csv
import re
from pathlib import Path

CONFIG = Path(__file__).resolve().parents[2] / "config"

# Legacy draft table label -> stable table ID (subject). Table S11 is retired (its sentence is
# replaced by authorized wording in build_manuscript); Table S2 is split into S2a/S2b.
LEGACY_TABLE = {
    "S1": "catalog", "S3": "outgroup", "S5": "loci", "S6": "synteny", "S7": "cazy_rank",
    "S8": "chromosome_gain_summary", "S9": "gc3", "S10": "selection_catalytic",
    "S12": "intron", "S13": "expression",
}
# Legacy draft figure label -> stable figure ID. Draft Figures S7 (gain percentiles) and S8
# (gain histogram) have no figure in the final supplement; see RETIRED_FIGURE_SUBJECTS.
LEGACY_FIGURE = {
    "S1": "euceroptres_genome", "S2": "tree_counts_full", "S3": "origin_trees_full",
    "S4": "accessory_domain_trees", "S5": "complete_locus_trees", "S6": "b_lectin_domain_tree",
}
# Legacy draft supplement section numeral -> final section ID, by section subject
# (draft supplement: 1 genome, 2 panel/transcriptomes, 3 outgroups, 4 gene trees,
# 5 B-lectin domain/locus trees, 6 orthology/synteny, 7 native copy number,
# 8 chromosome gains, 9 GC3 (and introns), 10 catalytic residues/selection).
LEGACY_SECTION = {
    "1": "S-A", "2": "S-B", "3": "S-C0", "4": "S-C", "5": "S-C", "6": "S-E",
    "7": "S-CN", "8": "S-CN", "9": "S-F", "10": "S-G",
}
# The draft's two retired-figure pointers cite gain values that the final supplement reports
# only in the chromosome-gain table. Exact contexts; each must occur once.
RETIRED_FIGURE_SUBJECTS = [
    ("(5 in Cynipini; Figure S8)", "Figure S8", "chromosome_gain_summary"),
    ("relative to all eCAMI sequence clusters (Figure S7)", "Figure S7", "chromosome_gain_summary"),
]

# Pinned author export: 59 source pointers on 34 lines. Two old section pointers
# on line 196 collapse to one S4 pointer in the generated prose.
AUTHOR_CALLSITES = {36: [('Table', 'catalog'), ('Supplement', 'S-A')],
 38: [('Supplement', 'S-B'), ('Table', 'catalog')],
 46: [('Supplement', 'S-C')],
 50: [('Table', 'outgroup'), ('Supplement', 'S-C0')],
 52: [('Figure', 'tree_counts_full'), ('Table', 'outgroup')],
 56: [('Figure', 'tree_counts_full')],
 60: [('Table', 'outgroup')],
 62: [('Supplement', 'S-C'), ('Figure', 'origin_trees_full')],
 63: [('Figure', 'accessory_domain_trees')],
 66: [('Figure', 'origin_trees_full')],
 70: [('Table', 'loci')],
 72: [('Figure', 'tree_counts_full')],
 74: [('Table', 'synteny'), ('Table', 'synteny')],
 76: [('Supplement', 'S-CN'),
      ('Table', 'cazy_rank'),
      ('Table', 'chromosome_gain_summary'),
      ('Table', 'chromosome_gain_summary'),
      ('Table', 'chromosome_gain_summary')],
 80: [('Table', 'intron')],
 90: [('Table', 'gc3'), ('Supplement', 'S-F')],
 94: [('Figure', 'complete_locus_trees'), ('Figure', 'complete_locus_trees')],
 97: [('Figure', 'complete_locus_trees')],
 99: [('Figure', 'b_lectin_domain_tree')],
 102: [('Figure', 'b_lectin_domain_tree'), ('Supplement', 'S-C')],
 106: [('Table', 'expression'), ('Table', 'expression')],
 112: [('Table', 'expression')],
 114: [('Table', 'expression')],
 118: [('Table', 'selection_catalytic'), ('Table', 'selection_catalytic'), ('Table', 'selection_catalytic')],
 122: [('Table', 'synteny')],
 168: [('Table', 'catalog'), ('Figure', 'euceroptres_genome'), ('Supplement', 'S-A')],
 178: [('Supplement', 'S-E')],
 182: [('Supplement', 'S-B')],
 190: [('Supplement', 'S-C0')],
 192: [('Supplement', 'S-C0')],
 196: [('Figure', 'origin_trees_full'),
       ('Figure', 'complete_locus_trees'),
       ('Figure', 'b_lectin_domain_tree'),
       ('Supplement', 'S-C')],
 202: [('Supplement', 'S-F')],
 206: [('Supplement', 'S-G'), ('Table', 'selection_catalytic')],
 210: [('Supplement', 'S-G')]}


POINTER = re.compile(
    r"\b(?P<kind>Tables?|Figures?|Supplements?)\s+"
    r"(?P<labels>S?\d+[a-z]?(?:(?:\s*,\s*|\s+and\s+|\s*[-–]\s*)S?\d+[a-z]?)*)")
LABEL = re.compile(r"S?\d+[a-z]?")


def _tsv(name: str) -> list[dict[str, str]]:
    with (CONFIG / name).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def final_maps() -> dict[str, dict[str, str]]:
    """Stable ID -> final label for tables, figures and sections (committed maps)."""
    tables = {r["table_id"]: r["label"] for r in _tsv("table_numbers.tsv")}
    # The expression table is split for the reader: S2a sample counts, S2b expression values.
    tables["expression"] = "S2b"
    tables["sample_counts"] = "S2a"
    figures = {r["figure_id"]: r["label"] for r in _tsv("figure_numbers.tsv") if r["label"].startswith("S")}
    sections = {r["section_id"]: r["number"] for r in _tsv("supplement_sections.tsv")}
    return {"Table": tables, "Figure": figures, "Supplement": sections}


def _remap_labels(base: str, labels: str, maps) -> list[str]:
    legacy = {"Table": LEGACY_TABLE, "Figure": LEGACY_FIGURE, "Supplement": LEGACY_SECTION}[base]
    out = []
    for old in LABEL.findall(labels):
        if base == "Figure" and not old.startswith("S"):
            out.append(old)  # main figure
        elif base == "Supplement" and old.startswith("S"):
            out.append(old)  # already final-form (only the replaced GC3 sentence in the draft)
        elif old in legacy:
            out.append(maps[base][legacy[old]])
        else:
            out.append(old)  # unresolved: left for the gate to reject
    return out


def remap(prose: str) -> str:
    maps = final_maps()

    def sub(m: re.Match) -> str:
        kind, labels = m.group("kind"), m.group("labels")
        base = kind.rstrip("s")
        new = _remap_labels(base, labels, maps)
        if base == "Supplement" and len(set(new)) == 1:
            # Two legacy sections now share one final section: "Supplements 4 and 5" -> "Supplement S4".
            return f"Supplement {new[0]}"
        # Keep the author's separators; substitute labels positionally.
        it = iter(new)
        return f"{kind} " + LABEL.sub(lambda _m: next(it), labels)

    prose = POINTER.sub(sub, prose)
    for context, old, sid in RETIRED_FIGURE_SUBJECTS:
        if prose.count(context) != 1:
            raise RuntimeError(f"retired-figure pointer context changed: {context}")
        prose = prose.replace(context, context.replace(old, "Table " + maps["Table"][sid]))
    return prose


def pointers(prose: str):
    """Yield (kind, label, start) for every supplement pointer token (main figures skipped)."""
    for m in POINTER.finditer(prose):
        base = m.group("kind").rstrip("s")
        for label in LABEL.findall(m.group("labels")):
            if base == "Figure" and not label.startswith("S"):
                continue
            yield base, label, m.start()


def check(prose: str) -> list[str]:
    """Return failures: pointers that do not resolve to the final supplement map."""
    maps = final_maps()
    final = {kind: set(labels.values()) for kind, labels in maps.items()}
    failures = []
    for kind, label, start in pointers(prose):
        if label not in final[kind]:
            context = prose[max(0, start - 60):start + 40].replace("\n", " ")
            failures.append(f"{kind} {label} does not resolve to the final map: ...{context}...")
    return failures


def check_callsites(prose: str) -> list[str]:
    """Check every pinned main-text callsite against its stable subject."""
    maps = final_maps()
    lines = prose.splitlines()
    failures = []
    for line, targets in AUTHOR_CALLSITES.items():
        expected = [(kind, maps[kind][target]) for kind, target in targets]
        observed = [(kind, label) for kind, label, _ in pointers(lines[line - 1])]
        if observed != expected:
            failures.append(f"line {line}: expected {expected}, observed {observed}")
    return failures


def check_subjects(prose: str) -> list[str]:
    """Check full-manuscript callsites whose wrong topic could still be in range."""
    failures = []
    # These subjects have in-range alternatives, so label existence alone is insufficient.
    subjects = (
        (r"Full statistics are presented in Supplement (S\d+)", "S7", "GC3 statistics"),
        (r"relative to all eCAMI sequence clusters \(Table (S\d+)\)", "S7", "chromosome gains"),
        (r"PCWDE paralog proliferation is associated with expression in multiple secretory tissues \(Table (S\d+[ab]?)\)", "S2b", "expression"),
    )
    for pattern, expected, subject in subjects:
        hits = re.findall(pattern, prose)
        if len(hits) != 1 or hits[0] != expected:
            failures.append(f"{subject} pointer must be {expected}; observed {hits}")
    return failures


def require(prose: str) -> None:
    failures = check(prose)
    if failures:
        raise RuntimeError("unresolved main-text supplement pointers:\n" + "\n".join(failures))
