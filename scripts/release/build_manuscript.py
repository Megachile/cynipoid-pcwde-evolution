#!/usr/bin/env python3
"""Typeset the manuscript text with the rendered figures placed and numbered from config/figure_numbers.tsv.

Text source: the pinned current draft (authority/manuscript/manuscript_source_pending.md,
Google Doc export of 2026-11-18T0422Z with its sync header and embedded images removed). In that
draft image anchor imageN sits above the Figure N legend, so figures are placed by figure
number. The draft carries no alt text; each figure's author-approved alt text is read from the
submitted manuscript record (manuscript_current.md, kept unchanged) by that figure's
submitted-text anchor in config/figure_numbers.tsv.
"""

from __future__ import annotations

import csv
import hashlib
import os
import re
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import figure_numbers  # noqa: E402
import main_text_pointers  # noqa: E402

ROOT = Path(__file__).resolve().parent
PAPER = ROOT.parents[1]
SUBMITTED = PAPER / "authority/manuscript/manuscript_current.md"
SUBMITTED_SHA = "5e6b453cbbde5f5cf4dc98d16057aff9286cdf7b619bb96e9fda83aa39394355"
SOURCE = PAPER / "authority/manuscript/manuscript_source_pending.md"
EXPECTED_DRAFT_SHA = "138e0fb97788587d3fafb77e8df4b6a314c2a38c39325310603b49a8777efc25"
FIG8_CAPTION = PAPER / "authority/manuscript/figure8_caption.md"
FIGURE_ROOT = Path(os.environ.get("PCWDE_MANUSCRIPT_FIGURES", PAPER / "build/figures")).resolve()
OUT = Path(os.environ.get("PCWDE_MANUSCRIPT_OUT", PAPER / "build/manuscript")).resolve()
TEMPLATE = PAPER / "config/manuscript_template.tex"
FIGURES = [(r["figure_id"], r["manuscript_anchor"], r["label"], r["path"])
           for r in figure_numbers.load() if r["manuscript_anchor"]]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def validate_frozen_authority() -> None:
    import json
    pending = PAPER / "authority/source_refresh_pending.json"
    if pending.exists() and json.loads(pending.read_text())["exceptions"]:
        raise RuntimeError("author source refresh has pending decisions; full build prohibited")
    if not SOURCE.is_file() or sha256(SOURCE) != EXPECTED_DRAFT_SHA:
        raise RuntimeError("manifested manuscript authority checksum mismatch")
    if SOURCE.read_text() != (PAPER / "authority/manuscript/manuscript_draft_20261118T0422Z.md").read_text():
        raise RuntimeError("manuscript source adoption is incomplete")
    if not SUBMITTED.is_file() or sha256(SUBMITTED) != SUBMITTED_SHA:
        raise RuntimeError("submitted manuscript record checksum mismatch")
    if sha256(SOURCE) == "22f987bcdd57e23ecea951f92773a7c92bd0497c52b5d0e9972e34d067c82df9":
        raise RuntimeError("stale manuscript authority rejected")


def load_manifest() -> list[dict[str, str]]:
    rows = []
    for content_id, anchor, number, relative in FIGURES:
        source = FIGURE_ROOT / relative
        if source.is_symlink() or not source.is_file() or source.stat().st_size == 0:
            raise RuntimeError(f"missing accepted rendered figure: {relative}")
        target = OUT / "figures" / Path(relative).name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        rows.append({
            "content_id": content_id, "anchor": f"image{number}", "submitted_anchor": anchor, "number": number,
            "asset": f"figures/{target.name}", "sha256": sha256(source),
        })
    return rows


def tex_escape(value: str, *, link_urls: bool = True) -> str:
    tokens: list[str] = []

    def token(tex: str) -> str:
        tokens.append(tex)
        return f"ZZZTOKEN{len(tokens)-1}ZZZ"

    # A backslash-escaped asterisk in the synchronized Markdown is literal
    # author text, not an emphasis delimiter.
    value = value.replace(r"\*", token("*"))
    value = re.sub(r"\\([\\`*{}_\[\]()#+.!<>~=:-])", r"\1", value)
    value = re.sub(
        r"<sup>(.*?)</sup>",
        lambda m: token(r"\textsuperscript{" + tex_escape(m.group(1)) + "}"),
        value,
    )
    value = re.sub(r"`([^`]+)`", lambda m: token(r"\texttt{" + tex_escape(m.group(1)) + "}"), value)
    value = re.sub(r"\[([^]]+)\]\((mailto:[^)]+)\)", lambda m: token(r"\mbox{\href{" + m.group(2) + "}{" + tex_escape(m.group(1), link_urls=False) + "}}"), value)
    value = re.sub(r"\[([^]]+)\]\(((?:https?://|mailto:)[^)]+)\)", lambda m: token(r"\href{" + m.group(2) + "}{" + tex_escape(m.group(1), link_urls=False) + "}"), value)
    # Markdown targets are protected first, so bare URLs cannot wrap them twice.
    # Sentence punctuation belongs to the prose; balanced URL parentheses remain.
    def bare_url(match: re.Match[str]) -> str:
        target = match.group()
        suffix = ""
        while target and (target[-1] in ".,;:!?" or
                          (target[-1] == ")" and target.count(")") > target.count("("))):
            suffix = target[-1] + suffix
            target = target[:-1]
        return token(r"\url{" + target + "}") + suffix

    if link_urls:
        value = re.sub(r"https?://[^\s<>\[\]{}]+", bare_url, value)
    value = re.sub(r"\*\*([^*]+)\*\*", lambda m: token(r"\textbf{" + tex_escape(m.group(1)) + "}"), value)
    value = re.sub(r"\*([^*]+)\*", lambda m: token(r"\textit{" + tex_escape(m.group(1)) + "}"), value)
    value = re.sub(r"\b10\^(-?\d+)", lambda m: token(r"\(10^{" + m.group(1) + r"}\)"), value)
    value = value.replace("$\\pm $", token(r"\(\pm\)"))
    superscripts = str.maketrans("⁻⁺⁰¹²³⁴⁵⁶⁷⁸⁹", "-+0123456789")
    value = re.sub(
        r"[⁻⁺⁰¹²³⁴⁵⁶⁷⁸⁹]+",
        lambda m: token(r"\textsuperscript{\ensuremath{" + m.group().translate(superscripts) + r"}}"),
        value,
    )
    replacements = {
        "\\": r"\textbackslash{}", "&": r"\&", "%": r"\%", "$": r"\$",
        "#": r"\#", "_": r"\_", "{": r"\{", "}": r"\}",
        "~": r"\textasciitilde{}", "^": r"\textasciicircum{}",
    }
    value = "".join(replacements.get(ch, ch) for ch in value)
    for i, saved in enumerate(tokens):
        value = value.replace(f"ZZZTOKEN{i}ZZZ", saved)
    return value


def clean_caption(raw: str) -> str:
    raw = raw.strip()
    raw = raw.replace("** **", " ")
    raw = re.sub(r"^\*\*Figure\s+[1-8]\\?\.\s*", "", raw)
    raw = raw.replace("**", "", 1)
    return raw.strip()


def figure_tex(row: dict[str, str], caption: str, alt_text: str) -> str:
    number = row["number"]
    asset = Path(row["asset"]).name
    return "\n".join([
        r"\par\nolinenumbers",
        r"\begin{figure}[!htbp]",
        r"\centering",
        r"\setbox\figurelegendbox=\vbox{\hsize=\textwidth",
        rf"\caption{{\textbf{{Figure {number}.}} {tex_escape(caption)}}}",
        rf"\label{{fig:{row['content_id']}}}",
        r"\small\raggedright\noindent\textbf{Alt text:} " + tex_escape(alt_text) + r"\par}",
        r"\setlength{\figureavailableheight}{\dimexpr\textheight-\ht\figurelegendbox-\dp\figurelegendbox-3\baselineskip\relax}",
        r"\ifdim\figureavailableheight>0.72\textheight\setlength{\figureavailableheight}{0.72\textheight}\fi",
        r"\ifdim\figureavailableheight<0.25\textheight\PackageError{pcwde}{Figure legend leaves insufficient image space}{}\fi",
        rf"\includegraphics[width=\textwidth,height=\figureavailableheight,keepaspectratio]{{{asset}}}",
        r"\par\medskip\usebox{\figurelegendbox}",
        r"\end{figure}",
        r"\linenumbers\par",
    ])


def submitted_alt_text() -> dict[str, str]:
    """Author-approved alt text beneath each figure legend of the submitted manuscript."""
    lines = SUBMITTED.read_text(encoding="utf-8").split("\n[image1]:", 1)[0].splitlines()
    alt, anchor = {}, None
    for line in lines:
        match = re.match(r"(?:\*\*)?!\[\]\[(image[1-8])\](?:\*\*)?", line.strip())
        if match:
            anchor = match.group(1)
        elif anchor and line.startswith("Alt text: "):
            alt[anchor] = line[len("Alt text: "):]
            anchor = None
    return alt


def apply_authorized_supplement_references(prose: str, stats_root: Path) -> str:
    """Project the pinned author draft onto the approved current table subjects."""
    import json
    from author_outgroup_counts import project
    numbers = json.loads((stats_root.parents[1] / 'numbers.json').read_text())
    prose, _ = project(prose, numbers, 'manuscript')
    from author_result_numbers import project as project_results
    prose, _ = project_results(prose, numbers, 'manuscript')
    # The current export already uses final labels. Only expression's unsplit
    # Table S2 pointers and retained legacy section forms need mechanical projection.
    prose = prose.replace("Table S2)", "Table S2b)")
    for old, new in (("Supplement 6\\.", "supplementary text S5."),
                     ("Supplement 3\\)", "supplementary text S3)")):
        prose = prose.replace(old, new)
    prose = re.sub(r"Supplement(?:ary)? (S?[1-9][0-9]*)", lambda m:
                   "supplementary text " + {"9": "S7"}.get(m.group(1), m.group(1)), prose)
    import json
    expected = json.loads((PAPER / "config/manuscript_current_pointers.json").read_text())
    observed = current_pointers(prose)
    if observed != expected:
        raise RuntimeError("current manuscript supplement pointer subjects changed")
    return prose


def current_pointers(prose: str) -> list:
    pattern = r"\b(Table|Figure|supplementary text) (S[0-9]+[ab]?)"
    return [[m.group(1), m.group(2)] for m in re.finditer(pattern, prose)]



def render(rows: list[dict[str, str]]) -> tuple[str, int]:
    raw = SOURCE.read_text(encoding="utf-8")
    alt_text = submitted_alt_text()
    prose = apply_authorized_supplement_references(raw.split("\n[image1]:", 1)[0], OUT.parent / "introns/gc3")
    if "data:image" in prose:
        raise RuntimeError("embedded payload leaked into extracted prose")
    lines = prose.splitlines()
    by_anchor = {r["anchor"]: r for r in rows}
    out: list[str] = []
    pending: dict[str, str] | None = None
    references = False
    placed: list[str] = []
    reference_count = 0
    first_heading = True
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        image_match = re.match(r"(?:\*\*)?!\[\]\[(image[1-8])\](?:\*\*)?", line)
        if image_match and image_match.group(1) not in by_anchor:
            # A figure dropped from the paper: skip its image, caption and alt text.
            j = i + 1
            while j < len(lines) and not lines[j].strip().startswith("Alt text: "):
                j += 1
            if j >= len(lines):
                raise RuntimeError(f"dropped figure {image_match.group(1)} lacks alt text")
            i = j + 1
            continue
        if image_match:
            pending = by_anchor[image_match.group(1)]
            i += 1
            continue
        if pending and line.startswith("**Figure"):
            caption = line
            caption = clean_caption(caption)
            j = i + 1
            while j < len(lines) and not lines[j].strip():
                j += 1
            if j < len(lines) and lines[j].startswith("Alt text: "):
                out.append(figure_tex(pending, caption, lines[j][len("Alt text: "):]))
                i = j + 1
            else:
                if pending["submitted_anchor"] not in alt_text:
                    raise RuntimeError(f"figure {pending['number']} lacks approved alt text")
                out.append(figure_tex(pending, caption, alt_text[pending["submitted_anchor"]]))
                i += 1
            placed.append(pending["number"])
            pending = None
            continue
        if not line:
            i += 1
            continue
        if line.startswith("# "):
            heading = line[2:].strip()
            if first_heading:
                out.extend([r"\begin{center}", rf"{{\LARGE\bfseries {tex_escape(heading)}}}\\[1.2em]"])
                first_heading = False
            elif heading == "References":
                references = True
                out.extend([r"\clearpage", r"\section*{References}", r"\addcontentsline{toc}{section}{References}", r"\begin{sloppypar}", r"\small"])
            else:
                out.extend([r"\FloatBarrier", rf"\section{{{tex_escape(heading)}}}"])
            i += 1
            continue
        if line.startswith("## "):
            heading = line[3:].strip()
            if heading == "References":
                references = True
                out.extend([r"\clearpage", r"\section*{References}", r"\addcontentsline{toc}{section}{References}", r"\begin{sloppypar}", r"\small"])
            else:
                out.extend([r"\FloatBarrier", rf"\subsection{{{tex_escape(heading)}}}"])
            i += 1
            continue
        if not out or (out and out[-1].startswith(r"{\LARGE")):
            author_line = re.sub(r"\*\*([1-6](?:,[1-6])*)\*\*", r"<sup>\1</sup>", line)
            out.append(tex_escape(author_line) + r"\\")
        elif line.startswith("1. Department"):
            out.append(r"\vspace{0.8em}")
            out.append(tex_escape(line) + r"\\")
        elif re.match(r"^[2-6]\. ", line):
            out.append(tex_escape(line) + r"\\")
            if line.startswith("6. "):
                out.extend([r"\end{center}", r"\vspace{1em}", r"\linenumbers"])
        elif references:
            reference_count += 1
            out.append(rf"\par\hangindent=1.5em\hangafter=1 {tex_escape(line)}")
        else:
            out.append(tex_escape(line) + "\n")
        i += 1
    if references:
        out.append(r"\end{sloppypar}")
    if pending:
        raise RuntimeError(f"image anchor lacked caption: {pending['anchor']}")
    if sorted(placed) != sorted(r["number"] for r in rows):
        raise RuntimeError(f"figures placed {placed} differ from the figure manifest")
    return "\n".join(out), reference_count


def compile_tex(stem: str) -> None:
    from check_pdf_glyphs import check_pdf_glyphs
    command = ["xelatex", "-interaction=nonstopmode", "-halt-on-error", f"{stem}.tex"]
    for pass_number in (1, 2):
        transcript = OUT / f"{stem}.xelatex.pass{pass_number}.stdout.txt"
        try:
            completed = subprocess.run(
                command,
                cwd=OUT,
                check=False,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
            )
        except OSError as exc:
            transcript.write_text(
                f"COMMAND: {shlex.join(command)}\nOS ERROR: {exc}\n",
                encoding="utf-8",
            )
            raise RuntimeError(f"XeLaTeX invocation failed; output retained at {transcript}") from exc
        transcript.write_text(
            f"COMMAND: {shlex.join(command)}\nRETURN CODE: {completed.returncode}\n"
            f"{completed.stdout}",
            encoding="utf-8",
        )
        if completed.returncode != 0:
            raise RuntimeError(
                f"XeLaTeX pass {pass_number} exited {completed.returncode}; "
                f"output retained at {transcript}"
            )
        check_pdf_glyphs(transcript, OUT / f"{stem}.pdf")


def check_figure_pages(pdf: Path, aux: Path, rows: list[dict[str, str]]) -> list[dict]:
    import fitz
    with fitz.open(pdf) as document:
        pages = [page.get_text() for page in document]
    labels = dict(re.findall(r"\\newlabel\{fig:([^}]+)\}\{\{[^}]*\}\{([0-9]+)\}", aux.read_text()))
    report = []
    for row in rows:
        figure_page = int(labels[row['content_id']])
        number = row['number']
        mention = re.compile(r"\bFigures?\s+" + re.escape(number) + r"\b")
        first_page = next((i for i, page in enumerate(pages, 1) if mention.search(page)), None)
        if first_page is None or not 0 <= figure_page - first_page <= 2:
            raise RuntimeError(f"Figure {number} page {figure_page}, first mention {first_page}: exceeds 2 pages")
        report.append(dict(number=number, figure_page=figure_page, first_mention_page=first_page))
    if [r['figure_page'] for r in report] != sorted(r['figure_page'] for r in report):
        raise RuntimeError("main figures are out of numerical order")
    return report


def apply_template(template: str, body: str) -> str:
    if template.count("@@BODY@@") != 1:
        raise RuntimeError("manuscript template must contain exactly one @@BODY@@ marker")
    return template.replace("@@BODY@@", body)


def main() -> None:
    validate_frozen_authority()
    if OUT.exists():
        raise RuntimeError(f"refusing existing manuscript output: {OUT}")
    OUT.mkdir(parents=True)
    rows = load_manifest()
    template = TEMPLATE.read_text(encoding="utf-8")
    body, refs = render(rows)
    (OUT / "manuscript.tex").write_text(apply_template(template, body), encoding="utf-8")
    compile_tex("manuscript")
    import json
    (OUT / "figure_page_gate.json").write_text(json.dumps(check_figure_pages(OUT / "manuscript.pdf", OUT / "manuscript.aux", rows), indent=2) + "\n")
    (OUT / "build_summary.txt").write_text(
        f"Manuscript SHA-256: {EXPECTED_DRAFT_SHA}\n"
        f"Submitted record SHA-256: {SUBMITTED_SHA}\n"
        f"manuscript.pdf: {len(rows)} accepted rendered figures; {refs} references\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
