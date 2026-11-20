#!/usr/bin/env python3
"""Deterministically typeset Table S1 (built by build_table_s1_rnaseq.py) with XeLaTeX."""

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import textwrap

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "release"))
import author_captions  # noqa: E402


ROOT = Path(__file__).resolve().parents[2]
DATA = Path(os.environ.get("PCWDE_REPOSITORY_INPUTS", ROOT / "figshare_data/repository_inputs")).resolve()
INPUT_ROOT = Path(os.environ.get("PCWDE_TABLE_S1_INPUT_DIR", DATA / "tables/accepted")).resolve()
INPUT = INPUT_ROOT / "Table_S1.tsv"
CAPTION = INPUT_ROOT / "Table_S1_caption.txt"
BUILD = Path(os.environ.get("PCWDE_BUILD_ROOT", ROOT / "build"))
DEFAULT_OUTPUT = BUILD / "tables/Table_S1.pdf"
WORK = BUILD / "table_s1_typeset/latex_build"

# These capacities are deliberately conservative approximations for TeX Gyre Termes
# at 9 pt in the fixed-width columns below.  The page budgets include headroom
# for rules and inter-row padding.  Pagination is decided here, before TeX is
# written; every generated tabular is therefore a complete, indivisible page.
WRAP_CAPACITIES = (18, 25, 29, 25, 40, 11, 12, 12)
FIRST_PAGE_LINE_BUDGET = 22
CONTINUATION_PAGE_LINE_BUDGET = 24


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def tex_escape(value: str) -> str:
    replacements = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "{": r"\{",
        "}": r"\}",
        "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
    }
    return "".join(replacements.get(char, char) for char in value)


def load_content():
    # The deposit Table S1 is pinned by config/logical_inputs.tsv; the built table it
    # becomes is checked cell by cell against it by scripts/release/check_table_s1.py.
    with INPUT.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        fields = reader.fieldnames
    expected_fields = [
        "resource_class", "species_display", "source_accession",
        "resource_preparation", "panel_status_exclusion_reason", "assembly_n50_bp",
        "assembly_busco_pct", "proteome_busco_pct",
    ]
    if fields != expected_fields or len(rows) != 99:
        raise SystemExit(f"Unexpected frozen table shape: fields={fields!r}, rows={len(rows)}")
    # Author-submitted caption (stable ID table:catalog) supersedes the deposit caption text.
    caption = author_captions.load().get("table:catalog") or CAPTION.read_text(encoding="utf-8").strip()
    if "*" in caption:
        raise SystemExit("Table S1 caption italic markup is not supported by this typesetter")
    return rows, caption


def wrapped_lines(value: str, capacity: int) -> int:
    """Conservatively estimate lines in a fixed-width TeX paragraph cell."""
    lines = textwrap.wrap(
        value,
        width=capacity,
        break_long_words=True,
        break_on_hyphens=True,
        replace_whitespace=True,
        drop_whitespace=True,
    )
    return max(1, len(lines))


def row_line_units(row) -> int:
    values = (
        row["resource_class"],
        row["species_display"],
        row["source_accession"],
        row["resource_preparation"],
        row["panel_status_exclusion_reason"],
        row["assembly_n50_bp"],
        row["assembly_busco_pct"],
        row["proteome_busco_pct"],
    )
    return max(wrapped_lines(value, capacity) for value, capacity in zip(values, WRAP_CAPACITIES))


def paginate(rows):
    pages = []
    current = []
    used = 0
    for row in rows:
        units = row_line_units(row)
        budget = FIRST_PAGE_LINE_BUDGET if not pages else CONTINUATION_PAGE_LINE_BUDGET
        if current and used + units > budget:
            pages.append(current)
            current = []
            used = 0
            budget = CONTINUATION_PAGE_LINE_BUDGET
        if units > budget:
            raise SystemExit(f"A frozen row cannot fit on one page: {row!r}")
        current.append(row)
        used += units
    if current:
        pages.append(current)
    return pages


FONT_HASHES = {
    "texgyretermes-regular.otf": "cc3fe7c707b81428d23d54df3eadd9228a2bf6a4d43125d94df56f5f63134659",
    "texgyretermes-bold.otf": "2fb3e952065fa153c7e4e64e04b98b9d79225739b6025aa3f0f0782d299ff61e",
    "texgyretermes-italic.otf": "6dd103a1672e50568cd2f8a706ccd48443d44d7d073a59d2286f4e6f746575d6",
    "texgyretermes-bolditalic.otf": "1bf6af99cb0e26c12951317032d79b96ae009551e59ccf02a5b24f325ecfec87",
}

def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def resolve_font_dir():
    override = os.environ.get("PCWDE_TEX_GYRE_DIR")
    if override:
        font_dir = Path(override)
    else:
        proc = subprocess.run(["kpsewhich", "texgyretermes-regular.otf"], text=True,
                              stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        if proc.returncode or not proc.stdout.strip():
            raise SystemExit("TeX Gyre Termes font resolution failed closed")
        font_dir = Path(proc.stdout.strip()).parent
    for filename, expected in FONT_HASHES.items():
        path = font_dir / filename
        if not path.is_file() or sha256(path) != expected:
            raise SystemExit(f"TeX Gyre Termes font identity failed closed: {filename}")
    return font_dir.resolve()

FONT_DIR = resolve_font_dir()
FONT_DECL = (rf"\setmainfont[Path={FONT_DIR.as_posix()}/,Extension=.otf,"
             r"UprightFont=*-regular,BoldFont=*-bold,ItalicFont=*-italic,"
             r"BoldItalicFont=*-bolditalic]{texgyretermes}")



def make_tex(rows, caption: str) -> str:
    pages = paginate(rows)

    def rows_tex(page_rows):
        body = []
        for row in page_rows:
            cells = [
                tex_escape(row["resource_class"]).replace("RNA-seq", r"\mbox{RNA-seq}"),  # never break inside RNA-seq
                r"\textit{" + tex_escape(row["species_display"]) + "}",
                tex_escape(row["source_accession"]),
                tex_escape(row["resource_preparation"]),
                tex_escape(row["panel_status_exclusion_reason"]),
                tex_escape(row["assembly_n50_bp"]),
                tex_escape(row["assembly_busco_pct"]),
                tex_escape(row["proteome_busco_pct"]),
            ]
            body.append(" & ".join(cells) + r" \\")
        return "\n".join(body)

    columns = r"@{}L{30mm}L{34mm}L{38mm}L{35mm}L{55mm}C{17mm}C{17mm}C{17mm}@{}"
    header = r"""\toprule
\textbf{Resource class} & \textbf{Species} & \textbf{Source/accession} &
\textbf{Resource preparation} & \textbf{Panel status / exclusion reason} & \textbf{Assembly N50 (bp)} &
\textbf{Assembly BUSCO complete (\%)} & \textbf{Proteome BUSCO complete (\%)} \\
\midrule"""
    blocks = []
    # The accepted fixed-page artifact supplies the table number as presentation
    # context without changing the caption authority.  Keeping only the first
    # caption sentence on continuation pages is also what gives every declared
    # page block the accepted vertical extent.
    title_label = "Table S1."
    title_context = f"{title_label} {caption.split('. ', 1)[0]}."
    for index, page_rows in enumerate(pages):
        page_title = f"{title_label} {caption}" if index == 0 else title_context
        continued = "" if index == 0 else r" \textit{(continued)}"
        heading = (rf"\section*{{Table S1}}\small {tex_escape(caption)}" if index == 0 else rf"\section*{{Table S1 (continued)}}\small {tex_escape(caption.split('. ', 1)[0] + '.')}")
        # Anchor every page in its own indivisible text-height box.  This makes
        # the page origin explicit instead of relying on TeX's vertical-page
        # builder to carry the next block across ``\clearpage``.  The latter
        # produced a fourth page that some native PDF renderers displayed with
        # its title and upper header clipped at the media edge.
        blocks.append(rf"""\noindent\begin{{minipage}}[t]{{\textwidth}}
\vspace{{0pt}}
\footnotesize
{heading}
\vspace{{1.2ex}}

\fontsize{{9.0}}{{11.0}}\selectfont
\noindent\begin{{tabular}}{{{columns}}}
{header}
{rows_tex(page_rows)}
\bottomrule
\end{{tabular}}\par
\end{{minipage}}""")
    tables_tex = "\n\\clearpage\n".join(blocks)
    return rf"""\documentclass[10pt]{{article}}
\usepackage[a4paper,landscape,margin=13mm,top=15mm,bottom=15mm]{{geometry}}
\usepackage{{fontspec}}
\usepackage{{booktabs,array,fancyhdr}}
\usepackage[table]{{xcolor}}
{FONT_DECL}
\pagestyle{{fancy}}
\fancyhf{{}}
\renewcommand{{\headrulewidth}}{{0pt}}
\fancyfoot[C]{{\footnotesize Table S1 --- Page \thepage}}
\setlength{{\parindent}}{{0pt}}
\setlength{{\tabcolsep}}{{2.5pt}}
\setlength{{\emergencystretch}}{{1.5em}}
\renewcommand{{\arraystretch}}{{1.15}}
\setlength{{\extrarowheight}}{{1pt}}
\definecolor{{headergray}}{{HTML}}{{E8EDF2}}
\newcommand{{\TableRaggedRight}}{{\rightskip=0pt plus 1.5em\leftskip=0pt\parfillskip=0pt plus 1fil\hyphenpenalty=50\exhyphenpenalty=50\relax}}
\newcolumntype{{L}}[1]{{>{{\TableRaggedRight\arraybackslash}}p{{#1}}}}
\newcolumntype{{C}}[1]{{>{{\centering\arraybackslash}}p{{#1}}}}
\begin{{document}}
{tables_tex}
\end{{document}}
"""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    output = args.output.resolve()
    try:
        output.relative_to(BUILD.resolve())
    except ValueError as exc:
        raise SystemExit("Output must remain beneath the declared build root") from exc

    rows, caption = load_content()
    WORK.mkdir(parents=True, exist_ok=True)
    tex_path = WORK / "Table_S1.tex"
    tex_path.write_text(make_tex(rows, caption), encoding="utf-8")
    env = os.environ.copy()
    env.update({
        "SOURCE_DATE_EPOCH": "1786579200",
        "FORCE_SOURCE_DATE": "1",
        "TZ": "UTC",
    })
    command = [
        "xelatex", "-no-shell-escape", "-interaction=nonstopmode",
        "-halt-on-error", "Table_S1.tex",
    ]
    for _ in range(2):
        # XeLaTeX is intentionally captured: stdout is the renderer's
        # single-document JSON interface, while the full TeX transcript remains
        # available in Table_S1.log for diagnostics.
        subprocess.run(
            command,
            cwd=WORK,
            env=env,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    built = WORK / "Table_S1.pdf"
    output.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(built, output)
    pages = paginate(rows)
    from pypdf import PdfReader
    actual_pages = len(PdfReader(output).pages)
    if actual_pages != len(pages):
        raise SystemExit(f"Table S1 pagination overflow: expected {len(pages)}, actual {actual_pages}")
    manifest = {
        "input_sha256": digest(INPUT),
        "caption_sha256": digest(CAPTION),
        "output_sha256": digest(output),
        "rows": len(rows),
        "pages": len(pages),
        "rows_per_page": [len(page) for page in pages],
        "estimated_line_units_per_page": [sum(row_line_units(row) for row in page) for page in pages],
    }
    print(json.dumps(manifest, sort_keys=True))


if __name__ == "__main__":
    main()
