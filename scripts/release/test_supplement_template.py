#!/usr/bin/env python3
"""Negative fixtures for the supplement prose gate (run by the build before rendering).

Each case copies the real template, injects one defect and runs
produce_supplement_sections.py; the run must fail with the expected message.
Unit cases cover formats and producer-constant checks directly.

Environment: project Python 3 (standard library only).
"""
import argparse
import json
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import supplement_template as st  # noqa: E402

p = argparse.ArgumentParser()
p.add_argument("--repo", type=Path, required=True)
p.add_argument("--numbers", type=Path, required=True)
a = p.parse_args()
repo = a.repo.resolve()
template = (repo / "authority/supplement/supplement_source.md").read_text()
numbers = json.loads(a.numbers.read_text())
failures = []
pending_path = repo / 'authority/source_refresh_pending.json'
if pending_path.exists() and json.loads(pending_path.read_text())['exceptions']:
    # A pending source decision must block fixture and production invocations alike.
    # Numeric negative fixtures below run after the source decision is resolved.
    with tempfile.TemporaryDirectory() as tmp:
        run = subprocess.run([sys.executable, '-B', str(repo / 'scripts/release/produce_supplement_sections.py'),
                              '--repo', str(repo), '--numbers', str(a.numbers),
                              '--out', str(Path(tmp) / 'out.json')], capture_output=True, text=True)
        assert run.returncode and 'pending decisions' in run.stderr
        assert not (Path(tmp) / 'out.json').exists()
    print('SUPPLEMENT_PENDING_SOURCE_GUARD_PASS; numeric negative fixtures deferred until source resolution')
    raise SystemExit(0)


TOTAL = "A total of {{supplement.genomes_considered}} genomes were considered"
cases = {
    "stale typed result": (template.replace(TOTAL, "A total of 70 genomes were considered", 1),
                           "numeral '70' is neither a token nor a listed literal"),
    "unknown token": (template.replace(TOTAL, "A total of {{supplement.no_such_key}} genomes were considered", 1),
                      "unknown or unresolved token"),
    "CHECK note": (template.replace(TOTAL, TOTAL + " \\[CHECK \\#99: fixture\\]", 1), "[CHECK note in prose"),
    "stale word-form result": (template.replace("{{supplement.transcriptome_synergini|word}} *Synergus*", "three *Synergus*", 1),
                               "numeral 'three' is neither a token nor a listed literal"),
    "stale literal row": (template.replace("crown in at least five genomes", "crown in at least some genomes", 1), "stale literal row"),
}
with tempfile.TemporaryDirectory() as tmp:
    for name, (text, message) in cases.items():
        if text == template:
            failures.append(f"fixture {name}: injection point not found in template")
            continue
        fixture = Path(tmp) / "template.md"
        fixture.write_text(text)
        run = subprocess.run([sys.executable, "-B", str(repo / "scripts/release/produce_supplement_sections.py"),
                              "--repo", str(repo), "--numbers", str(a.numbers), "--template", str(fixture),
                              "--out", str(Path(tmp) / "out.json")], capture_output=True, text=True)
        if run.returncode == 0 or message not in run.stderr:
            failures.append(f"fixture {name}: gate did not fail with {message!r}")
    bad = dict(numbers, supplement=dict(numbers["supplement"], busco_floor_pct=45.0))
    bad_path = Path(tmp) / "numbers.json"
    bad_path.write_text(json.dumps(bad))
    run = subprocess.run([sys.executable, "-B", str(repo / "scripts/release/produce_supplement_sections.py"),
                          "--repo", str(repo), "--numbers", str(bad_path), "--out", str(Path(tmp) / "out.json")],
                         capture_output=True, text=True)
    if run.returncode != 0 or "below a 45% BUSCO" not in (Path(tmp) / "out.json").read_text():
        failures.append("fixture producer-backed BUSCO parameter did not follow changed producer")

# Exclusive result counts ("only N ..."): they fail with no source, and a descriptive
# registry row does not exempt them. Ordinary descriptive literals, valid tokens and
# the authorized author_report rows still pass.
literals = st.read_literals(repo / "config/supplement_literals.tsv")
paragraphs = template.split("\n\n")


def sentence_errors(snippet, section, extra=()):
    """Gate errors for the template sentence containing snippet, with its own literal rows."""
    para = [p for p in paragraphs if snippet in p]
    if len(para) != 1:
        failures.append(f"fixture sentence not found once: {snippet!r}")
        return None
    text = para[0]
    lo = max(text.rfind(". ", 0, text.index(snippet)) + 2, 0) if ". " in text[:text.index(snippet)] else 0
    hi = text.find(". ", text.index(snippet) + len(snippet))
    sentence = text[lo:hi + 1 if hi >= 0 else len(text)]
    rows = [r for r in literals if r["section"] == section and r["context"] in sentence] + list(extra)
    return st.check_section(section, [sentence], rows, numbers, repo)[0]


# The two formerly unsourced exhaustive claims now have exact author assertions.
for section, snippet in (("S-B", "Since only three of the genomes"),
                         ("S-C0", "only one outgroup sequence contained")):
    if sentence_errors(snippet, section) is None:
        failures.append(f"asserted claim fixture missing: {snippet}")
VALID = (("S-E", "Two methods were used to determine"), ("S-B", "Three kinds of cynipid larvae"),
         ("S-A", "of which 30 (0.0011% of reads) were discarded"), ("S-B", "Read-sorting identified two larval samples"),
         ("S-B", "A total of {{supplement.genomes_considered}} genomes"))
for section, snippet in VALID:
    errors = sentence_errors(snippet, section)
    if errors:
        failures.append(f"fixture valid sentence {snippet!r} fails: {errors[:2]}")
with tempfile.TemporaryDirectory() as tmp:
    path = Path(tmp) / "literals.tsv"
    path.write_text("\t".join(st.LITERAL_COLUMNS) + "\nS-A\t99\tof which 99 reads\tauthor_report\t\t\n")
    try:
        st.read_literals(path)
        failures.append("fixture unauthorised author_report row was accepted")
    except SystemExit:
        pass
    path.write_text("\t".join(st.LITERAL_COLUMNS) + "\nS-B\ttwo\tWe found only two pseudogenes.\tauthor_report\t\t\n")
    try:
        st.read_literals(path)
        failures.append("fixture approved numeral in wrong author claim was accepted")
    except SystemExit:
        pass
for row in [r for r in literals if r["class"] == "author_report"]:
    if not st.claim_covers(template, template.find(row["context"]) + row["context"].index(row["numeral"]),
                           template.find(row["context"]) + row["context"].index(row["numeral"]) + len(row["numeral"]), row["context"]):
        failures.append(f"approved author claim not matched: {row['numeral']}")
for phrase, section, word in (("Since only **three** of the genomes were previously annotated.", "S-B", "three"),
                              ("However, only *one* outgroup sequence contained the full intact domain architecture.", "S-C0", "one")):
    row = {"section": section, "numeral": word, "context": phrase, "class": "descriptive", "check_key": ""}
    if not any("exclusive result count" in e for e in st.check_section(section, [phrase], [row], numbers, repo)[0]):
        failures.append(f"emphasized unresolved claim accepted: {phrase}")
method = "We used only two methods: RagTag and JCVI."
row = {"section": "S-E", "numeral": "two", "context": method, "class": "descriptive", "check_key": ""}
if st.check_section("S-E", [method], [row], numbers, repo)[0]:
    failures.append("legitimate only-two-methods claim rejected")
# Approved author claims. Every author_report/author_assertion row is bound to its own
# complete claim. Fixtures use the real template paragraphs with every literal row of the
# section that occurs in them, so the positive baseline must be completely clean and any
# rejection below is attributable to the assertion itself.
BOUND = [r for r in literals if r["class"] in ("author_report", "author_assertion")]
key = lambda r: (r["section"], r["numeral"], r["context"])  # noqa: E731
claims = {}
for row in BOUND:
    claims.setdefault((row["section"], st.AUTHOR_CLAIM_BINDINGS[key(row)]), []).append(row)
claim_paras = {}
for section, claim in claims:
    found = [q for q in paragraphs if claim in q]
    if len(found) != 1:
        failures.append(f"approved claim not found once in template: {claim[:60]!r}")
        continue
    claim_paras.setdefault(section, [])
    if found[0] not in claim_paras[section]:
        claim_paras[section].append(found[0])


def fixture_rows(section, paras):
    return [r for r in literals if r["section"] == section and any(r["context"] in q for q in paras)]


def numeral_error(section, para, value, start, end):
    lo = max(0, start - 40)
    return (f"{section}: numeral {value!r} is neither a token nor a listed literal: "
            f"...{para[lo:end + 25]}...")


def stale_error(row):
    return f"{row['section']}: stale literal row no longer in the template: {row['numeral']!r} in {row['context']!r}"


def gate(section, paras, rows):
    return st.check_section(section, paras, rows, numbers, repo)[0]


real_covers = st.asserted_row_covers


def disabled(claim):
    """In-memory mutation: rows bound to claim match by section and numeral only."""
    def covers(row, paragraph, start, end):
        if st.AUTHOR_CLAIM_BINDINGS.get(key(row)) == claim:
            return True
        return real_covers(row, paragraph, start, end)
    return covers


# One biological-predicate change per claim; the numeral and every supporting literal
# context (2021, 45/two, 4, ...) stay unchanged.
PREDICATES = {
    "of which 30 (0.0011% of reads)": ("HiFiAdapterFilt", "Cutadapt"),
    "reported raw and QV scores": ("PacBio HiFi reads", "Illumina reads"),
    "The genome coverage was estimated": ("GenomeScope2", "Merqury"),
    "Read-sorting identified two larval": ("chalcidoid genomes", "bacterial genomes"),
    "EGAPx annotated 14,633": ("lncRNA", "tRNA"),
    "Since only three of the genomes": ("previously annotated", "previously sequenced"),
    "We found only one case": ("*Parnips nigripes*", "*Parnips fictus*"),
    "However, only one outgroup": ("*Actinoplanes awajinensis*", "*Actinoplanes fictus*"),
    "These include a lineage": ("one oomycete acquisition", "one fungal acquisition"),
}
claim_cases = 0
for section, paras in claim_paras.items():
    rows = fixture_rows(section, paras)
    errors = gate(section, paras, rows)
    if errors:
        failures.append(f"approved claim baseline not clean in {section}: {errors[:2]}")
        continue
    claim_cases += 1
    for (sec, claim), bound in claims.items():
        if sec != section:
            continue
        para = next(q for q in paras if claim in q)
        others = [q for q in paras if q is not para]
        # Claim deletion: exactly this claim's rows (and supporting rows inside it) go stale.
        # Supporting rows whose context crossed the deleted sentence boundary may also fail.
        cut = para.replace(claim, "", 1)
        gone = [r for r in rows if not any(r["context"] in q for q in others + [cut])]
        got = gate(section, others + [cut], rows)
        explained = {stale_error(r) for r in gone} | {e for e in got if any(
            e.startswith(f"{section}: numeral {r['numeral']!r} is neither") for r in gone)}
        if (not all(stale_error(r) in got for r in bound) or set(got) - explained
                or any(stale_error(r) in got for r in BOUND if r in rows and r not in bound + gone)):
            failures.append(f"claim deletion not attributed to its own rows: {claim[:50]!r}: {sorted(set(got) - explained)[:2]}")
        claim_cases += 1
        # Changed biological predicate: the asserted occurrences themselves are rejected.
        prefix = next((k for k in PREDICATES if k in claim), None)
        if prefix is None:
            failures.append(f"no predicate fixture for approved claim {claim[:50]!r}")
            continue
        old, new = PREDICATES[prefix]
        changed_claim = claim.replace(old, new, 1)
        if changed_claim == claim or any(old in r["context"] for r in rows if r not in bound):
            failures.append(f"predicate fixture changes a supporting literal or nothing: {prefix}")
            continue
        changed = para.replace(claim, changed_claim, 1)
        lo, hi = changed.index(changed_claim), changed.index(changed_claim) + len(changed_claim)
        values = {r["numeral"] for r in bound}
        expected = {numeral_error(section, changed, v, a_, b_) for v, a_, b_ in st.numerals(changed)
                    if v in values and lo <= a_ and b_ <= hi}
        expected |= {stale_error(r) for r in bound}
        if len(expected) < 2 * len(bound) or set(gate(section, others + [changed], rows)) != expected:
            failures.append(f"changed predicate not rejected at its asserted occurrence: {prefix}")
        # Predicate-disabled control: the same negative fixture must then pass, so the
        # rejection above comes from this claim binding and nothing else.
        st.asserted_row_covers = disabled(claim)
        try:
            if gate(section, others + [changed], rows):
                failures.append(f"negative fixture does not depend on its predicate binding: {prefix}")
        finally:
            st.asserted_row_covers = real_covers
        claim_cases += 2
    # Cross-claim reuse: a row bound to one claim is never satisfied by another approved
    # claim with the same section and numeral.
    for (sec_a, claim_a), bound_a in claims.items():
        for (sec_b, claim_b), bound_b in claims.items():
            if sec_a != section or sec_b != section or claim_a == claim_b:
                continue
            for row_b in bound_b:
                if not any(r["numeral"] == row_b["numeral"] for r in bound_a):
                    continue
                para = next(q for q in paras if claim_a in q)
                support = [r for r in fixture_rows(section, [para]) if r not in bound_a]
                got = gate(section, [para], support + [row_b])
                missing = [r for r in bound_a if r["numeral"] == row_b["numeral"]]
                if stale_error(row_b) not in got or not all(f"numeral {r['numeral']!r} is neither" in " ".join(got) for r in missing):
                    failures.append(f"cross-claim reuse accepted: {row_b['numeral']!r} from {claim_b[:40]!r} in {claim_a[:40]!r}")
                claim_cases += 1
if sorted(r["section"] for r in BOUND if r["class"] == "author_assertion") != ["S-B", "S-C0", "S-C0", "S-C0", "S-C0"]:
    failures.append("author_assertion registry rows changed")
for good in ("{{supplement.genomes_considered}}", "{{supplement.genomes_considered|word}}"):
    if st.render(good, numbers)[2]:
        failures.append(f"valid token rejected: {good}")
for malformed in ("{{supplement.genomes_considered}", "{supplement.genomes_considered}}",
                  "{{{supplement.genomes_considered}}}", "{{supplement.genomes_considered}}}",
                  "{{{supplement.genomes_considered}}", "{{supplement.genomes_considered}}}}",
                  "{{{{supplement.genomes_considered}}}}", "{{supplement.{{genomes_considered}}}}"):
    if not st.render(malformed, numbers)[2]:
        failures.append(f"malformed token accepted: {malformed}")
for value, spec, want in ((8, "word", "eight"), (4, "Word", "Four"), (920372, "comma", "920,372"), (0.0924, "f2", "0.09"), (839.0, "f0", "839"),
                         (0.013791745, "f3", "0.014"), (100.0, "f1t", "100"), (99.0826, "f1t", "99.1"), (42, "word", "forty-two")):
    if st.fmt(value, spec) != want:
        failures.append(f"format {spec}({value}) != {want}")
for value, spec in ((15.3, None), (2.5, "word"), (100, "word")):
    try:
        st.fmt(value, spec)
        failures.append(f"format {spec}({value}) should be refused")
    except ValueError:
        pass
if [n for n, _, _ in st.numerals("IQ-TREE v3.0.1 LG+G4 odb12 E139 GCF\\_010883055.1 n=10 (p < 0.05)")] != ["10", "0.05"]:
    failures.append("numeral scanner treats identifiers as numerals")
if [n for n, _, _ in st.numerals("Twenty-one loci in two-stage runs; one-third; someone; a dozen")] != \
        ["Twenty-one", "two", "one", "dozen"]:
    failures.append("word-form scanner misses or over-matches written-out counts")
# Section F GC3 tokens: notation formats, the cross-family producer, and a live binding
# (a changed GC3 statistic must change the rendered prose, never leave the old value).
import supplement_numbers  # noqa: E402
gc3_cases = 0
for value, spec, want in ((3.0517578125e-05, "sci0", "3 × 10⁻⁵"), (9.6e-3, "sci0", "1 × 10⁻²"), (2.4e-7, "sci0", "2 × 10⁻⁷"),
                          (-10, "pow10", "10⁻¹⁰"), (5, "pow10", "10⁵")):
    gc3_cases += 1
    if st.fmt(value, spec) != want:
        failures.append(f"format {spec}({value}) != {want}")
for value, spec in ((-1.0, "sci0"), (0, "sci0"), (-10.5, "pow10")):
    gc3_cases += 1
    try:
        st.fmt(value, spec)
        failures.append(f"format {spec}({value}) should be refused")
    except ValueError:
        pass
gc3 = numbers["gc3"]
facts = supplement_numbers.gc3_facts(gc3)
gc3_cases += 1
pcwde_p = [gc3["pcwde"][f]["wilcoxon_bacteria_minus_host_distance_p"] for f in ("GH5_2", "PL1_8", "PL4_4")]
if facts["pcwde_copies"] != sum(gc3["pcwde"][f]["copies"] for f in ("GH5_2", "PL1_8", "PL4_4")) or \
        not max(pcwde_p) < 10 ** facts["pcwde_wilcoxon_p_upper_exponent"] <= 10 * max(pcwde_p):
    failures.append("gc3_facts sum or P bound wrong")
gc3_cases += 1
edge = json.loads(json.dumps(gc3))
edge["pcwde"]["PL4_4"]["wilcoxon_bacteria_minus_host_distance_p"] = 1e-10
if supplement_numbers.gc3_facts(edge)["pcwde_wilcoxon_p_upper_exponent"] != -9:
    failures.append("gc3_facts P bound not strict at an exact power of ten")
gc3_cases += 1
missing = json.loads(json.dumps(gc3))
del missing["pcwde"]["PL4_4"]
try:
    supplement_numbers.gc3_facts(missing)
    failures.append("gc3_facts accepted a missing family")
except SystemExit:
    pass
gc3_para = [q for q in paragraphs if "GC3 was measured for" in q]
if len(gc3_para) != 1:
    failures.append("GC3 paragraph not found once in template")
else:
    gc3_cases += 1
    base_text, _, base_errors = st.render(gc3_para[0], numbers)
    moved = json.loads(json.dumps(numbers))
    moved["gc3"]["pcwde"]["PL4_4"]["host_distance_pp"] = 2.04
    moved["gc3"]["blectin"]["all"]["wilcoxon_distance_p_two_sided"] = 4.4e-6
    moved["supplement"]["gc3"] = supplement_numbers.gc3_facts(moved["gc3"])
    moved_text, _, moved_errors = st.render(gc3_para[0], moved)
    if base_errors or moved_errors or "within 2.0–" not in moved_text or "P \\= 4 × 10⁻⁶" not in moved_text \
            or moved_text == base_text:
        failures.append("GC3 paragraph does not follow its producer values")
if failures:
    sys.exit("SUPPLEMENT_TEMPLATE_TEST_FAIL\n" + "\n".join(failures))
print(f"SUPPLEMENT_TEMPLATE_TEST_PASS cases={len(cases) + 1 + 8 + len(VALID) + 1} claim_cases={claim_cases} gc3_cases={gc3_cases}")
