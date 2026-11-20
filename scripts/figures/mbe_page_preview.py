#!/usr/bin/env python3
"""
mbe_page_preview.py — show a draft figure on a mock MBE journal page at the
size it would actually print, with its caption and surrounding manuscript text.

Preview only: it never edits the figure, the caption or the manuscript.

Page model (MBE, Oxford Academic): letter paper, 17.5 cm x 23.5 cm text block,
two 8.4 cm body columns, 9 pt serif body text, 8 pt caption. A figure that fits
the text block prints at 100%; a figure wider than the block but fitting
sideways is rotated 90 degrees (as a landscape float); anything else is scaled
down to fit, and the scale is reported in the margin note and in the JSON
sidecar.

Usage:
  mbe_page_preview.py --fig FIG.pdf --figure-number 4 --out OUT_STEM \
      [--manuscript draft.md] [--caption-match TEXT] [--label TEXT]

Writes OUT_STEM.pdf, OUT_STEM.png (300 dpi, anti-aliased) and OUT_STEM.json.
Run in the repository build environment with XeLaTeX and Ghostscript.
"""
import argparse
import json
import re
import subprocess
import sys
import tempfile
from pathlib import Path


CM = 72 / 2.54
TEXT_W_CM, TEXT_H_CM = 17.5, 23.5
COL_W_CM, GUTTER_CM = 8.4, 0.7
CAPTION_RESERVE_CM = 3.0   # room kept for the caption under a full-page figure
MIN_BODY_CM = 5.0          # below this, a figure page carries no body text


def pdf_size_pt(path):
    out = subprocess.run(
        ["gs", "-q", "-dNODISPLAY", "-dNOSAFER", "-c",
         f"({path}) (r) file runpdfbegin 1 pdfgetpage /MediaBox pget pop == quit"],
        check=True, text=True, capture_output=True).stdout
    x0, y0, x1, y1 = (float(v) for v in re.findall(r"[-\d.]+", out)[:4])
    return x1 - x0, y1 - y0


def md_plain(s):
    """Google Docs Markdown export -> plain text with *emphasis* kept."""
    s = s.replace("$\\pm$", "±").replace("$\\geq$", "≥").replace("$\\leq$", "≤")
    s = re.sub(r"!\[\]\[image\d+\]", "", s)
    s = re.sub(r"\\([\\`*_{}\[\]()#+\-.!<>|~])", r"\1", s)
    s = s.replace("****", "")
    return s.strip()


def tex_escape(s):
    s = s.replace("\\", r"\textbackslash{}")
    for a, b in [("&", r"\&"), ("%", r"\%"), ("$", r"\$"), ("#", r"\#"),
                 ("_", r"\_"), ("{", r"\{"), ("}", r"\}"),
                 ("~", r"\textasciitilde{}"), ("^", r"\textasciicircum{}")]:
        s = s.replace(a, b)
    s = re.sub(r"\*\*(.+?)\*\*", r"\\textbf{\1}", s)
    s = re.sub(r"\*(.+?)\*", r"\\emph{\1}", s)
    return s


def read_manuscript(path):
    lines = [md_plain(ln) for ln in Path(path).read_text().splitlines()
             if not ln.startswith("[image")]
    captions = {}
    for ln in lines:
        m = re.match(r"^\*\*Figure (\d+)\.", ln)
        if m:
            captions.setdefault(int(m.group(1)), []).append(ln)
    return lines, captions


def pick_caption(captions, fig_no, match):
    cands = captions.get(fig_no, [])
    if match:
        cands = [c for c in cands if match.lower() in c.lower()]
    if len(cands) != 1:
        listing = "\n".join("  " + c[:110] for c in captions.get(fig_no, []))
        sys.exit(f"Figure {fig_no}: {len(cands)} matching captions; pass --caption-match. Candidates:\n{listing}")
    return cands[0]


def body_paragraphs(lines, fig_no, n_words):
    """Results paragraphs that cite the figure, then the following ones."""
    cite = re.compile(rf"\bFigures? {fig_no}\b|\bFig\. ?{fig_no}\b|\bFigures? [\dS, and]*\b{fig_no}\b")
    paras = [ln for ln in lines
             if ln and not ln.startswith("#") and not ln.startswith("**Figure")
             and not ln.startswith("|") and len(ln.split()) > 25]
    start = next((i for i, p in enumerate(paras) if cite.search(p)), 0)
    out, words = [], 0
    for p in paras[start:]:
        out.append(p)
        words += len(p.split())
        if words >= n_words:
            break
    return out


def place(w_pt, h_pt):
    """Upright unless the figure is wider than the text block and fits better sideways."""
    tw, th = TEXT_W_CM * CM, TEXT_H_CM * CM
    s_up = min(1.0, tw / w_pt, (th - CAPTION_RESERVE_CM * CM) / h_pt)
    s_rot = min(1.0, th / w_pt, (tw - CAPTION_RESERVE_CM * CM) / h_pt)
    return "rotated" if (w_pt > tw and s_rot > s_up * 1.02) else "upright"


def build_tex(fig, mode, w_pt, h_pt, caption, paras, note):
    """XeLaTeX measures the caption and prints the figure at 100% unless the
    figure plus caption cannot fit; the final scale goes to the log."""
    body = "\n\n".join(tex_escape(p) for p in paras)
    # box width/height before rotation: rotated figures run along \textheight
    boxw, boxh = (r"\textwidth", r"\textheight") if mode == "upright" else (r"\textheight", r"\textwidth")
    head = rf"""\documentclass[9pt]{{extarticle}}
\usepackage[letterpaper,textwidth={TEXT_W_CM}cm,textheight={TEXT_H_CM}cm,centering]{{geometry}}
\usepackage{{fontspec}}
\setmainfont{{texgyretermes-regular.otf}}[BoldFont=texgyretermes-bold.otf,ItalicFont=texgyretermes-italic.otf,BoldItalicFont=texgyretermes-bolditalic.otf]
\usepackage{{graphicx,multicol,xcolor,fancyhdr}}
\setlength{{\columnsep}}{{{GUTTER_CM}cm}}
\setlength{{\parindent}}{{1em}}\setlength{{\parskip}}{{0pt}}
\pagestyle{{fancy}}\fancyhf{{}}\renewcommand{{\headrulewidth}}{{0pt}}
\fancyhead[C]{{\textcolor{{gray}}{{\fontsize{{7}}{{8}}\selectfont {tex_escape(note)}}}}}
\fancyfoot[C]{{\textcolor{{gray}}{{\fontsize{{7}}{{8}}\selectfont Molecular Biology and Evolution — mock page}}}}
\newsavebox\capbox\newlength\availh\newlength\figw\newlength\figh
\begin{{document}}
\sbox\capbox{{\parbox{{{boxw}}}{{\fontsize{{8}}{{9.6}}\selectfont {tex_escape(caption)}}}}}
\setlength\availh{{{boxh}}}\addtolength\availh{{-\ht\capbox}}\addtolength\availh{{-\dp\capbox}}\addtolength\availh{{-14pt}}
\setlength\figw{{{w_pt:.3f}pt}}\setlength\figh{{{h_pt:.3f}pt}}
\ifdim\figw>{boxw}\setlength\figh{{\dimexpr\figh*{boxw}/\figw\relax}}\setlength\figw{{{boxw}}}\fi
\ifdim\figh>\availh\setlength\figw{{\dimexpr\figw*\availh/\figh\relax}}\setlength\figh{{\availh}}\fi
\typeout{{PREVIEW-FIGW=\the\figw}}
\newcommand{{\figblock}}{{\begin{{minipage}}{{{boxw}}}\centering\includegraphics[width=\figw,height=\figh]{{{fig}}}\par\vspace{{6pt}}\usebox\capbox\end{{minipage}}}}
"""
    if mode == "rotated":
        # rotated block fills the page; MBE prints landscape figures as their own page
        return head \
            + "\\noindent\\vbox to\\textheight{\\vfil\\hbox to\\textwidth{\\hfil\\rotatebox{90}{\\figblock}\\hfil}\\vfil}\n\\end{document}\n"
    return head + "\\noindent\\figblock\n\n\\vspace{10pt}\n\\begin{multicols}{2}\\fontsize{9}{11}\\selectfont\n" \
        + body + "\n\\end{multicols}\n\\end{document}\n"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fig", required=True)
    ap.add_argument("--figure-number", type=int, required=True)
    ap.add_argument("--out", required=True, help="output stem, no extension")
    ap.add_argument("--manuscript", required=True)
    ap.add_argument("--caption-match", default=None,
                    help="substring selecting the caption when the draft has several for this number")
    ap.add_argument("--label", default=None, help="short name for the margin note")
    a = ap.parse_args()

    fig = Path(a.fig).resolve()
    out = Path(a.out).resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    w_pt, h_pt = pdf_size_pt(fig)
    mode = place(w_pt, h_pt)

    lines, caps = read_manuscript(a.manuscript)
    caption = pick_caption(caps, a.figure_number, a.caption_match)
    paras = body_paragraphs(lines, a.figure_number, 900)

    label = a.label or fig.stem
    def note_for(scale):
        return (f"Print preview: {label} — {w_pt/72:.2f} x {h_pt/72:.2f} in figure printed at "
                f"{scale*100:.0f}%{' (rotated page)' if mode == 'rotated' else ''}; "
                f"caption: current draft Figure {a.figure_number}")

    with tempfile.TemporaryDirectory() as td:
        tf = Path(td) / "page.tex"
        scale = None
        for _ in range(2):  # pass 1 measures the scale, pass 2 prints it in the margin note
            tf.write_text(build_tex(fig, mode, w_pt, h_pt, caption, paras,
                                    note_for(scale) if scale else note_for(1.0)))
            p = subprocess.run(["xelatex", "-interaction=nonstopmode", "-halt-on-error", "page.tex"],
                               cwd=td, text=True, capture_output=True)
            if p.returncode:
                sys.stderr.write(p.stdout[-3000:])
                sys.exit("XeLaTeX failed")
            m = re.search(r"PREVIEW-FIGW=([\d.]+)pt", p.stdout)
            scale = float(m.group(1)) / w_pt
        pdf = out.with_suffix(".pdf")
        pdf.write_bytes((Path(td) / "page.pdf").read_bytes())
    pages = int(subprocess.run(
        ["gs", "-q", "-dNODISPLAY", "-dNOSAFER", "-c", f"({pdf}) (r) file runpdfbegin pdfpagecount == quit"],
        check=True, text=True, capture_output=True).stdout.split()[0])
    subprocess.run(["gs", "-q", "-dSAFER", "-dBATCH", "-dNOPAUSE", "-sDEVICE=png16m", "-r300", "-dTextAlphaBits=4", "-dGraphicsAlphaBits=4",
                    "-dFirstPage=1", "-dLastPage=1", f"-sOutputFile={out.with_suffix('.png')}", str(pdf)],
                   check=True)
    meta = {"figure": str(fig), "figure_in": [round(w_pt / 72, 3), round(h_pt / 72, 3)],
            "mode": mode, "print_scale": round(scale, 4), "pdf_pages": pages,
            "text_block_cm": [TEXT_W_CM, TEXT_H_CM], "caption_figure": a.figure_number, "caption_match": a.caption_match,
            "manuscript": str(a.manuscript), "body_paragraphs": len(paras)}
    out.with_suffix(".json").write_text(json.dumps(meta, indent=2) + "\n")
    print(json.dumps(meta))


if __name__ == "__main__":
    main()
