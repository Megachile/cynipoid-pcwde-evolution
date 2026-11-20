"""Render supplement prose tokens from numbers.json and check every remaining numeral.

Template syntax: {{key.path}} or {{key.path|format}}. Keys are dotted paths into
numbers.json (e.g. supplement.PL1_8.analysis_domains, synteny_modal_blocks.
GH5_2:chr10:BLOCK_0000.detected). Formats: comma (thousands separator), f0/f1/f2/f3
(fixed decimals), f1t (one decimal, a trailing ".0" dropped: 100.0 -> 100), word /
Word (integers 0-99 as English words). A float without a format is refused, so no
value is silently rounded.

Every numeral left in a paragraph after token substitution must be listed in
config/supplement_literals.tsv (section, numeral, context, class, check_key, source_note). A
listed numeral with a check_key must equal that value: a numbers.json path, or
const:<repo script>:<NAME>[.index...] for a module-level constant in a producer
script (a percentage numeral may equal a fraction constant x 100). Each literal
row must still match the template (stale rows fail). Result numbers therefore
cannot be typed into the prose: a typed count is an unlisted numeral.

Written-out counts (one ... ninety-nine, "a dozen") are gated the same way: each
must come from a word-format token or be a listed literal. Word literals may also
use class "descriptive" (a count of items named in the same sentence, or a
description of study design), which digit numerals may not. A count right after
"only" ("only three of the genomes ...", "only one outgroup sequence contained
...") is an exhaustive result over a data set, even when the sentence names the
items, so it may not be descriptive: it needs a computed token or an approved author assertion.

Environment: project Python 3 (standard library only).
"""
import ast
import csv
import math
import re

TOKEN = re.compile(r"\{\{([^{}|]+)(?:\|([^{}]+))?\}\}")
# A numeral not glued to a preceding letter, underscore, digit or point (so v3.0.1,
# odb12, E139, LG+G4 and GCF_010883055.1 are identifiers, not numerals).
NUMERAL = re.compile(r"(?<![A-Za-z_\d.])\d(?:[\d,.]*\d)?")
_UNITS = ("zero one two three four five six seven eight nine ten eleven twelve thirteen "
          "fourteen fifteen sixteen seventeen eighteen nineteen").split()
_TENS = "twenty thirty forty fifty sixty seventy eighty ninety".split()
WORDS = _UNITS + [t if u == 0 else f"{t}-{_UNITS[u]}" for t in _TENS for u in range(10)]
WORD_VALUES = {w: i for i, w in enumerate(WORDS)}
WORD_VALUES["dozen"] = 12
# A written-out count (not part of a hyphenated compound such as "two-stage" is
# still a count; "one-third" style fractions are caught by "one").
WORD_NUMERAL = re.compile(r"(?<![A-Za-z-])(?:" + "|".join(sorted(WORD_VALUES, key=len, reverse=True)) +
                          r")(?![A-Za-z])", re.I)
LITERAL_CLASSES = {"parameter", "version", "date", "citation", "protocol", "descriptive", "author_report", "author_assertion"}
# Full approved sentences, including the biological predicate and disposition.
_HIFI = 'of which 30 (0.0011% of reads) were discarded after filtering with HiFiAdapterFilt.'
_QV = 'K-mer analysis of the final assembly relative to the PacBio HiFi reads used to create the contig assemblies reported raw and QV scores of 65.867 and 63.392 respectively.'
_COVERAGE = 'The genome coverage was estimated at 44.6× coverage as reported by GenomeScope2.'
_EGAPX = 'EGAPx annotated 14,633 genes, including {{supplement.e_primus_proteins|comma}} protein-coding genes and 1,737 lncRNA.'
_READ_SORT = 'Read-sorting identified two larval samples with predominant hits to chalcidoid genomes, which were removed from analysis.'
_ANNOTATED = 'Since only three of the genomes (*B. kinseyi*, *Leptopilina boulardi* and *L. heterotoma*) were previously annotated, Helixer v0.3.6 (model invertebrate\\_v0.3\\_m\\_0100; Stiehler et al. 2021\\) was used to predict genomic elements.'
_PSEUDOGENE = 'We found only one case for which a gene family was represented in a genome only by a pseudogene. *Parnips nigripes* contains a PL1\\_8-like region whose tblastn hits fell below the 45% identity threshold, which was characterized separately by protein-to-genome alignment of two cynipid PL1\\_8 proteins (*Aulacidea tavakolii* and *Ceroptres masudai*) with miniprot v0.9 and exonerate v2.4.0 (protein2genome), followed by a dbCAN-sub hmmsearch of the translated alignment.'
_ARCHITECTURE = 'However, only one outgroup sequence contained the full intact domain architecture including PI-PLC-X, GH30\\_3, and tandem B-lectin: bacterial species *Actinoplanes awajinensis*.'
_ACQUISITIONS = 'These include a lineage of Nematodes possessing bacterial GH5\\_2, two distinct bee acquisitions (Figure S3), one midge acquisition, and one oomycete acquisition of PL1\\_8, and acquisitions of PL4\\_4 by weevils, *Aspergillus* molds, and the parasitic plant *Cuscuta*.'
# Each author_report / author_assertion registry row (section, numeral, context) is bound
# to its own approved claim. The row covers a numeral only where its own claim and its own
# context both cover that occurrence, so one row can never be satisfied (or marked used)
# by another approved claim with the same section and numeral.
AUTHOR_CLAIM_BINDINGS = {
    ('S-A', '30', 're obtained, of which 30 (0.0011% of reads) were '): _HIFI,
    ('S-A', '0.0011', 'btained, of which 30 (0.0011% of reads) were discarde'): _HIFI,
    ('S-A', '65.867', ' raw and QV scores of 65.867 and 63.392 respectively.'): _QV,
    ('S-A', '63.392', ' scores of 65.867 and 63.392 respectively. The genome'): _QV,
    ('S-A', '44.6', 'rage was estimated at 44.6× coverage as reported by'): _COVERAGE,
    ('S-A', '14,633', 'gous. EGAPx annotated 14,633 genes, including {{suppl'): _EGAPX,
    ('S-A', '1,737', 'and 1,737 lncRNA.'): _EGAPX,
    ('S-B', 'two', 'ad-sorting identified two larval samples with pred'): _READ_SORT,
    ('S-B', 'three', _ANNOTATED): _ANNOTATED,
    ('S-C0', 'one', _PSEUDOGENE): _PSEUDOGENE,
    ('S-C0', 'one', _ARCHITECTURE): _ARCHITECTURE,
    ('S-C0', 'two', _ACQUISITIONS): _ACQUISITIONS,
    ('S-C0', 'one', _ACQUISITIONS): _ACQUISITIONS,
}
AUTHOR_ASSERTED_ROWS = set(AUTHOR_CLAIM_BINDINGS)
LITERAL_COLUMNS = ["section", "numeral", "context", "class", "check_key", "source_note"]
UNRESOLVED_CLAIMS = ("only three of the genomes", "only one outgroup sequence contained")


def claim_text(text):
    """Markdown-insensitive text and a map back to original character offsets."""
    chars, offsets = [], []
    for i, char in enumerate(text):
        if char in "*_\\`":
            continue
        chars.append(char)
        offsets.append(i)
    return "".join(chars), offsets


def claim_covers(paragraph, start, end, context):
    clean, offsets = claim_text(paragraph)
    needle, _ = claim_text(context)
    for match in re.finditer(re.escape(needle), clean):
        if offsets[match.start()] <= start and end <= offsets[match.end() - 1] + 1:
            return True
    return False


def unresolved_claim(paragraph, start, end):
    clean, offsets = claim_text(paragraph)
    for phrase in UNRESOLVED_CLAIMS:
        for match in re.finditer(re.escape(phrase), clean, re.I):
            if offsets[match.start()] <= start and end <= offsets[match.end() - 1] + 1:
                return True
    return False


def asserted_row_covers(row, paragraph, start, end):
    """An author_report/author_assertion row covers this occurrence only inside its own claim."""
    claim = AUTHOR_CLAIM_BINDINGS.get((row["section"], row["numeral"], row["context"]))
    return bool(claim) and claim_covers(paragraph, start, end, claim) and \
        claim_covers(paragraph, start, end, row["context"])


def token_boundaries(paragraph):
    """Reject malformed delimiters before any token can erase their evidence."""
    errors = []
    spans = list(TOKEN.finditer(paragraph))
    for match in re.finditer(r"[{}]+", paragraph):
        if not any(span.start() <= match.start() and match.end() <= span.end()
                   for span in spans):
            errors.append("malformed or unclosed token delimiter")
            break
    for span in spans:
        before = paragraph[span.start() - 1:span.start()] if span.start() else ""
        after = paragraph[span.end():span.end() + 1]
        if (before and before in "{}") or (after and after in "{}"):
            errors.append("malformed or unclosed token delimiter")
            break
    return errors


def lookup(numbers, path):
    value = numbers
    for part in path.split("."):
        if not isinstance(value, dict) or part not in value:
            raise KeyError(path)
        value = value[part]
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise KeyError(path)
    return value


_SUPERSCRIPT = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")


def _superscript(exponent):
    return str(int(exponent)).translate(_SUPERSCRIPT)


def fmt(value, spec):
    if spec == "absf1":
        return f"{abs(value):.1f}"
    if spec in (None, ""):
        if isinstance(value, float):
            raise ValueError("float token needs a format")
        return str(value)
    if spec in ("f0", "f1", "f2", "f3"):
        return f"{value:.{spec[1]}f}"
    if spec == "f1t":
        text = f"{value:.1f}"
        return text[:-2] if text.endswith(".0") else text
    if spec == "sci0":
        # P = 3 × 10⁻⁵ style: one-significant-figure mantissa, superscript exponent.
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not value > 0:
            raise ValueError("format sci0 needs a positive number")
        exponent = math.floor(math.log10(value))
        mantissa = round(value / 10 ** exponent)
        if mantissa == 10:
            mantissa, exponent = 1, exponent + 1
        return f"{mantissa} × 10{_superscript(exponent)}"
    if not (isinstance(value, int) or (isinstance(value, float) and value.is_integer())):
        raise ValueError(f"format {spec} needs an integer")
    value = int(value)
    if spec == "comma":
        return f"{value:,}"
    if spec == "pow10":
        # Integer exponent e rendered as 10ᵉ (bounds such as P < 10⁻¹⁰).
        return f"10{_superscript(value)}"
    if spec in ("word", "Word"):
        if not 0 <= value < len(WORDS):
            raise ValueError("word format supports 0-99")
        return WORDS[value].capitalize() if spec == "Word" else WORDS[value]
    raise ValueError(f"unknown format {spec}")


def render(paragraph, numbers):
    """Return (rendered text, [(key, format, value text)], [error])."""
    used, errors = [], token_boundaries(paragraph)
    if errors:
        return paragraph, used, errors

    def sub(match):
        key, spec = match.group(1).strip(), match.group(2)
        try:
            text = fmt(lookup(numbers, key), spec)
        except KeyError:
            errors.append(f"unknown or unresolved token {{{{{match.group(0)[2:-2]}}}}}")
            return match.group(0)
        except ValueError as exc:
            errors.append(f"token {match.group(0)}: {exc}")
            return match.group(0)
        used.append((key, spec or "", text))
        return text

    rendered = TOKEN.sub(sub, paragraph)
    if "{{" in rendered or "}}" in rendered:
        errors.append("malformed or unclosed token delimiter")
    return rendered, used, errors


def masked(paragraph):
    """Template text with every token blanked to same-length filler (keeps offsets)."""
    return TOKEN.sub(lambda m: "#" * len(m.group(0)), paragraph)


def numerals(paragraph):
    """Digit numerals and written-out counts left in a template paragraph, in text order."""
    text = masked(paragraph)
    found = [(m.group(), m.start(), m.end()) for m in NUMERAL.finditer(text)]
    found += [(m.group(), m.start(), m.end()) for m in WORD_NUMERAL.finditer(text)]
    return sorted(found, key=lambda x: x[1])


def numeral_value(numeral):
    word = WORD_VALUES.get(numeral.lower())
    return float(word) if word is not None else float(numeral.replace(",", ""))


def read_literals(path):
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != LITERAL_COLUMNS:
            raise SystemExit(f"{path}: columns must be {LITERAL_COLUMNS}")
        rows = list(reader)
    keys = [(r["section"], r["numeral"], r["context"]) for r in rows]
    if len(set(keys)) != len(keys):
        raise SystemExit(f"{path}: duplicate literal rows")
    bad = [r for r in rows if r["class"] not in LITERAL_CLASSES or r["numeral"] not in r["context"]
           or (r["class"] == "descriptive" and r["numeral"].lower() not in WORD_VALUES)
           or (r["class"] in ("author_report", "author_assertion") and
               (r["section"], r["numeral"], r["context"]) not in AUTHOR_ASSERTED_ROWS)]
    if bad:
        raise SystemExit(f"{path}: invalid literal rows {bad[:3]}")
    return rows


def producer_value(check_key, numbers, repo):
    if not check_key.startswith("const:"):
        return lookup(numbers, check_key)
    _, script, name = check_key.split(":", 2)
    name, *index = name.split(".")
    tree = ast.parse((repo / script).read_text())
    values = [ast.literal_eval(node.value) for node in tree.body if isinstance(node, ast.Assign)
              and any(isinstance(t, ast.Name) and t.id == name for t in node.targets)]
    if len(values) != 1:
        raise KeyError(check_key)
    value = values[0]
    for i in index:
        value = value[int(i)]
    return value


def check_section(section_id, paragraphs, literals, numbers, repo):
    """Gate the numerals of one section's template paragraphs. Returns (errors, used rows)."""
    errors, used = [], set()
    rows = [(i, r) for i, r in enumerate(literals) if r["section"] == section_id]
    rows_by_index = dict(rows)
    for paragraph in paragraphs:
        for value, start, end in numerals(paragraph):
            covering = [i for i, r in rows if r["numeral"] == value and
                        (asserted_row_covers(r, paragraph, start, end) if r["class"] in ("author_report", "author_assertion")
                         else claim_covers(paragraph, start, end, r["context"]))]
            lo = max(0, start - 40)
            if not covering:
                errors.append(f"{section_id}: numeral {value!r} is neither a token nor a listed literal: "
                              f"...{paragraph[lo:end + 25]}...")
            elif unresolved_claim(paragraph, start, end) and all(
                    rows_by_index[i]["class"] == "descriptive" for i in covering):
                errors.append(f"{section_id}: exclusive result count {value!r} cannot be a descriptive literal: "
                              f"...{paragraph[lo:end + 25]}...")
            used.update(covering)
    for i, r in rows:
        if i not in used:
            errors.append(f"{section_id}: stale literal row no longer in the template: {r['numeral']!r} in {r['context']!r}")
        elif r["check_key"]:
            try:
                expected = producer_value(r["check_key"], numbers, repo)
            except (KeyError, OSError, ValueError, IndexError):
                errors.append(f"{section_id}: literal check_key has no producer value: {r['check_key']}")
                continue
            stated = numeral_value(r["numeral"])
            percent = r["context"][r["context"].index(r["numeral"]) + len(r["numeral"]):].startswith("%")
            if stated != float(expected) and not (percent and abs(stated / 100 - float(expected)) < 1e-9):
                errors.append(f"{section_id}: parameter {r['numeral']} disagrees with producer value "
                              f"{r['check_key']}={expected}")
    return errors, used
