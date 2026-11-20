#!/usr/bin/env python3
"""Assemble the accepted supplement from this checkout and its declared data root."""

import hashlib
import json
import csv
import os
import math
import re
import subprocess
import sys
from pathlib import Path

import fitz  # PyMuPDF (environment/python-lock.txt): page geometry and layout checks

sys.path.insert(0, str(Path(__file__).resolve().parent))
import figure_numbers  # noqa: E402
import author_captions  # noqa: E402
from build_manuscript import tex_escape  # same Markdown-to-TeX rules as manuscript
from check_pdf_glyphs import check_pdf_glyphs
from blectin_domain_counts import verify_counts as verify_blectin_domain_counts
from supplement_page_numbers import number_pages, check_numbers

repo = Path(__file__).resolve().parents[2]
build = Path(os.environ["PCWDE_BUILD_ROOT"]).resolve()
data = Path(os.environ["PCWDE_DATA_ROOT"]).resolve()
out = build / "supplement"
if out.exists():
    raise SystemExit("refusing existing supplement output")
s5 = data / "build_inputs/euceroptres_genome/euceroptres_genome.pdf"
if s5.is_symlink() or not s5.is_file() or hashlib.sha256(s5.read_bytes()).hexdigest() != "c968c62fb81a2da159945ffe4286170174c264ef759ae00a8328f82cd77eb78f":
    raise SystemExit("missing or mismatched pinned Figure S5")


def supplement_tex(template, rows, build, data, domain_counts):
    """Fill figure labels and paths from the number map; the template names figures only by stable id."""
    definitions = "\n".join(
        f"\\expandafter\\def\\csname figpath{r['figure_id']}\\endcsname{{{figure_numbers.source_path(r, build, data).as_posix()}}}\n"
        f"\\expandafter\\def\\csname figlabel{r['figure_id']}\\endcsname{{{r['label']}}}" for r in rows)
    tex = template.replace("@@FIGURE_DEFINITIONS@@", definitions)
    for r in rows:
        tex = tex.replace(f"@@FIG_{r['figure_id']}@@", figure_numbers.source_path(r, build, data).as_posix())
        tex = tex.replace(f"@@LABEL_{r['figure_id']}@@", r["label"])
    captions = author_captions.load()
    for r in rows:
        token = f"@@CAPTION_{r['figure_id']}@@"
        if token in tex:
            caption = captions.get(f"figure:{r['figure_id']}")
            if caption is None or tex.count(token) != 1:
                raise SystemExit(f"figure caption missing or duplicated for stable ID figure:{r['figure_id']}")
            tex = tex.replace(token, author_captions.to_tex(caption))
    for key in ("total_domains", "cynipoid_domains", "cynipoid_proteins"):
        # Producer counts fill the caption only where the caption text names them.
        token = f"@@BLECTIN_{key.upper()}@@"
        if tex.count(token) > 1:
            raise SystemExit(f"Figure S6 caption token duplicated: {token}")
        tex = tex.replace(token, str(domain_counts[key]))
    if "@@" in tex:
        raise SystemExit("supplement template names a figure missing from config/figure_numbers.tsv")
    return tex


# A section with no moved text whose first item is a table carries its heading
# (and author slot) on that table's first page, so no heading stands alone.
# The page is moved down, and only if needed scaled, to reserve a heading band.
HEAD_TOP, HEAD_LINE, SLOT_LINE, HEAD_GAP, BOTTOM_MARGIN, LEFT = 36.0, 18.0, 15.0, 14.0, 18.0, 39.6
FONT_FLOOR = 6.0


def plain(text):
    return text.replace("*", "")


def page_content(page):
    """Displayed-orientation bbox of all text and drawings, and the smallest font size."""
    m = page.rotation_matrix
    box = fitz.Rect()
    sizes = []
    for block in page.get_text("dict")["blocks"]:
        if block.get("lines"):
            box.include_rect(fitz.Rect(block["bbox"]) * m)
            sizes += [span["size"] for line in block["lines"] for span in line["spans"] if span["text"].strip()]
    for drawing in page.get_drawings():
        box.include_rect(drawing["rect"] * m)
    return box, min(sizes)


def heading_overlay(table_pdf, title, slot):
    doc = fitz.open(table_pdf)
    page = doc[0]
    width, height = page.rect.width, page.rect.height  # displayed size, rotation applied
    box, min_font = page_content(page)
    lines = max(1, math.ceil(len(plain(title)) * 8.2 / (width - 2 * LEFT)))
    band = HEAD_TOP + lines * HEAD_LINE + (SLOT_LINE + 4 if slot else 0) + HEAD_GAP
    room = height - BOTTOM_MARGIN - band
    scale = min(1.0, room / box.height)
    if min_font * scale < FONT_FLOOR:
        raise SystemExit(f"heading band would shrink {table_pdf.name} text below {FONT_FLOOR} pt")
    # pdfpages centres the scaled page; offset moves it so its content starts below the band.
    top_after = height / 2 + scale * (box.y0 - height / 2)
    shift = max(0.0, band - top_after)
    if height / 2 + scale * (box.y1 - height / 2) + shift > height - BOTTOM_MARGIN + 0.01:
        raise SystemExit(f"heading band does not fit on {table_pdf.name}")
    block = r"\normalfont\Large\bfseries " + tex_escape(title) + r"\par"
    if slot:
        block += r"\vspace{4bp}\normalsize\bfseries " + tex_escape(slot) + r"\par"
    # pdfpages options cannot hold \par, so the heading block is a macro defined first.
    macro = "\\supphead" + "ABCDEFGHIJ"[len(overlays)]
    picture = (r"\setlength{\unitlength}{1bp}\put(%.2f,%.2f){\parbox[t]{%.2fbp}{\raggedright %s}}"
               % (LEFT, height - HEAD_TOP - 12, width - 2 * LEFT, macro))
    tex = (f"\\newcommand{{{macro}}}{{{block}}}\n"
           f"\\includepdf[pages=1,fitpaper=true,templatesize={{{width:.2f}bp}}{{{height:.2f}bp}},"
           f"scale={scale:.4f},offset=0bp {-shift:.2f}bp,picturecommand*={{{picture}}}]{{../tables/{table_pdf.name}}}")
    if len(doc) > 1:
        tex += f"\n\\includepdf[pages=2-,fitpaper=true]{{../tables/{table_pdf.name}}}"
    record = dict(table=table_pdf.stem, heading=plain(title), slot=slot, scale=round(scale, 4), shift_bp=round(shift, 2),
                  band_bp=round(band, 2), min_font_pt=round(min_font * scale, 2))
    return tex, record


supplement = [r for r in figure_numbers.load() if r["label"].startswith("S") and r["source"] != "reserved"]
for row in supplement:
    if not figure_numbers.source_path(row, build, data).is_file():
        raise SystemExit(f"missing supplement figure {row['figure_id']}")
with (repo / "config/table_numbers.tsv").open(newline="") as handle:
    table_labels = {r["table_id"]: r["label"] for r in csv.DictReader(handle, delimiter="\t")}
old_table_ids = {1:"catalog",2:"outgroup",3:"loci",4:"synteny",5:"intron",6:"expression",7:"selection_catalytic",8:"gc3",9:"catalytic_clustering",10:"copy_expansion",11:"native_comparators",12:"locus_conflict",13:"di_sort_limits"}
old_figure_ids = {1:"origin_trees_full",2:"accessory_domain_trees",3:"complete_locus_trees",4:"b_lectin_domain_tree",5:"euceroptres_genome"}
figure_labels = {r["figure_id"]:r["label"] for r in figure_numbers.load()}
def remap_refs(s):
    s = re.sub(r"\bFigures (\d+) and S(\d+)\b",
               lambda m: f"Figures {m.group(1)} and {figure_labels[old_figure_ids[int(m.group(2))]]}", s)
    def replace(m):
        kind, old = m.group(1), int(m.group(2))
        if kind == "Table":
            stable = old_table_ids.get(old)
            return f"Table {table_labels[stable]}" + ("b" if stable == "expression" else "") if stable in table_labels else m.group(0)
        stable = old_figure_ids.get(old)
        return f"Figure {figure_labels[stable]}" if stable else m.group(0)
    return re.sub(r"\b(Table|Figure) S(\d+)\b", replace, s)
for n in sorted(int(label[1:]) for label in table_labels.values()):
    if not (build / "tables" / ("Table_S2b.pdf" if n == 2 else f"Table_S{n}.pdf")).is_file():
        raise SystemExit(f"missing rendered supplement table S{n}")
if not (build / "tables/Table_S2a.pdf").is_file():
    raise SystemExit("missing rendered sample-count panel S2a")
source = repo / "config/supplement_assembly.tex"
section_text = json.loads((build / "supplement_section_text.json").read_text())
by_section = {row["section_id"]: row["paragraphs"] for row in section_text}
with (repo / "config/supplement_sections.tsv").open(newline="") as handle:
    sections = list(csv.DictReader(handle, delimiter="\t"))
expected = {f"S{n}" for n in range(1,11)}
if set(table_labels.values()) != expected or len(table_labels) != 10:
    raise SystemExit("supplement table map must contain S1–S10 exactly once")
figure_ids = {r["figure_id"] for r in supplement}
items_seen = []
parts = []
overlays = []
for section in sections:
    # Topic sections are labelled by id; the trailing tables part carries its heading only.
    title = section["heading"] if section["section_id"] == "additional" else section["number"] + ". " + section["heading"]
    items = list(filter(None, section["items"].split(",")))
    pending_heading = None
    if not by_section.get(section["section_id"]) and items and items[0].startswith("table:"):
        pending_heading = (title, section["author_slot"])
    else:
        parts.append("\\clearpage\n\\section*{" + tex_escape(title) + "}")
        if section["author_slot"]:
            parts.append("\\noindent\\textbf{" + tex_escape(section["author_slot"]) + "}\\par")
    for paragraph in by_section.get(section["section_id"], []):
        parts.append(tex_escape(paragraph) + "\\par")
    for item in items:
        if item in items_seen: raise SystemExit(f"duplicate supplement item: {item}")
        items_seen.append(item)
        if item.startswith("figure:"):
            fid = item.split(":",1)[1]
            if fid not in figure_ids: raise SystemExit(f"unknown figure {fid}")
            parts.append(remap_refs((repo / "config" / f"supplement_{fid}.tex").read_text()))
        elif item.startswith("table:"):
            table_id = item.split(":",1)[1]
            if table_id not in table_labels:
                raise SystemExit(f"unknown or reserved table {table_id}")
            n = int(table_labels[table_id][1:])
            if table_id == "expression":
                parts.append(r"\includepdf[pages=-,fitpaper=true]{../tables/Table_S2a.pdf}")
                parts.append(r"\includepdf[pages=-,fitpaper=true]{../tables/Table_S2b.pdf}")
                continue
            if pending_heading:
                tex_part, record = heading_overlay(build / "tables" / f"Table_S{n}.pdf", *pending_heading)
                parts.append(tex_part)
                overlays.append(record)
                pending_heading = None
            else:
                parts.append(f"\\includepdf[pages=-,fitpaper=true]{{../tables/Table_S{n}.pdf}}")
        else: raise SystemExit(f"unknown supplement item {item}")
if {f"figure:{i}" for i in figure_ids} - set(items_seen):
    raise SystemExit("a supplementary figure has no section")
if {f"table:{table_id}" for table_id in table_labels} - set(items_seen):
    raise SystemExit("a supplementary table has no section")
out.mkdir(parents=True, exist_ok=True)
(out / "section_text.json").write_text(json.dumps(section_text, ensure_ascii=False, indent=2) + "\n")
tex = source.read_text()
if tex.count("@@SECTIONED_BODY@@") != 1: raise SystemExit("missing section body slot")
domain_tree = build / "figures/b_lectin_domain_tree/primary/trees/Figure_S4_B_lectin_domains.rooted.treefile"
domain_manifest = build / "figures/b_lectin_domain_tree/primary/sidecars/Figure_S4_domain_manifest.tsv"
domain_counts = json.loads((build / "numbers.json").read_text())["blectin_domain_tree"]
verify_blectin_domain_counts(domain_tree, domain_manifest, domain_counts)
tex = supplement_tex(tex.replace("@@SECTIONED_BODY@@", "\n".join(parts)), supplement, build, data, domain_counts)
(out / "supplement.tex").write_text(tex)
for pass_number in (1, 2):
    proc = subprocess.run(["xelatex", "-interaction=nonstopmode", "-halt-on-error", "supplement.tex"], cwd=out, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    (out / f"xelatex.pass{pass_number}.txt").write_text(proc.stdout)
    if proc.returncode:
        raise SystemExit(f"supplement XeLaTeX pass {pass_number} failed")
    check_pdf_glyphs(out / f"xelatex.pass{pass_number}.txt", out / "supplement.pdf")
if not (out / "supplement.pdf").is_file():
    raise SystemExit("supplement PDF missing")

# pdfpages rounds the tall Figure S1 page sizes and shifts their artwork by
# roughly half a point. Replace those three TeX placeholders with cloned source
# pages, retaining the surrounding supplement pages and their authored text.
origin = next(r for r in supplement if r["figure_id"] == "origin_trees_full")
origin_source = fitz.open(figure_numbers.source_path(origin, build, data))
if len(origin_source) != 3:
    raise SystemExit("origin_trees_full must have three pages")
assembled = fitz.open(out / "supplement.pdf")
caption_pages = [i for i, page in enumerate(assembled) if
                 f"Figure {origin['label']}. " + author_captions.plain(author_captions.load()["figure:origin_trees_full"]).split(". ")[0] + "." in " ".join(page.get_text().split())]
if len(caption_pages) != 1:
    raise SystemExit(f"origin_trees_full caption mapping ambiguous: {caption_pages}")
hits = list(range(caption_pages[0] + 1, caption_pages[0] + 4))
if hits[-1] >= len(assembled):
    raise SystemExit("origin_trees_full source pages missing")
for j, index in enumerate(hits):
    if not any(assembled.xref_stream(x[0]) == origin_source[j].read_contents()
               for x in assembled[index].get_xobjects()):
        raise SystemExit(f"origin_trees_full source page {j+1} not at expected position")
replaced = fitz.open()
for i in range(len(assembled)):
    if i in hits:
        replaced.insert_pdf(origin_source, from_page=i - hits[0], to_page=i - hits[0])
    else:
        replaced.insert_pdf(assembled, from_page=i, to_page=i)
replaced_path = out / "supplement.replaced.pdf"
replaced.save(replaced_path, garbage=4, deflate=True)
assembled.close()
replaced.close()
replaced_path.replace(out / "supplement.pdf")
check_pdf_glyphs(out / "xelatex.pass2.txt", out / "supplement.pdf")

pdf_check = fitz.open(out / "supplement.pdf")
for j, index in enumerate(hits):
    src_page, dst_page = origin_source[j], pdf_check[index]
    if src_page.rect != dst_page.rect or src_page.read_contents() != dst_page.read_contents():
        raise SystemExit(f"S1 page {j+1}: geometry or content stream changed")
    if dst_page.get_xobjects() != src_page.get_xobjects():
        # XRef numbers may be remapped by insert_pdf; compare resource streams below.
        if len(dst_page.get_xobjects()) != len(src_page.get_xobjects()):
            raise SystemExit(f"S1 page {j+1}: resource count changed")
    for scale in (1, 4):
        matrix = fitz.Matrix(scale, scale)
        if src_page.get_pixmap(matrix=matrix, alpha=False).samples != dst_page.get_pixmap(matrix=matrix, alpha=False).samples:
            raise SystemExit(f"S1 page {j+1}: {scale}x render changed")
print(f"S1_EXACT_PASS pages={','.join(str(i+1) for i in hits)}")
pdf_check.close()
origin_source.close()

# The deposited table PDFs retain their standalone footers. Only the assembled
# supplement receives its continuous numbering, after the S1 source-page clone.
numbered = fitz.open(out / "supplement.pdf")
number_pages(numbered)
check_numbers(numbered)
numbered.save(out / "supplement.numbered.pdf", garbage=4, deflate=True)
numbered.close()
(out / "supplement.numbered.pdf").replace(out / "supplement.pdf")
print(f"NUMBERS_PASS pages={len(fitz.open(out / 'supplement.pdf'))}")

# Layout checks on the assembled PDF: every text span inside its page, and each
# overlaid heading alone in its band, above all table content, at >= 6 pt.
pdf = fitz.open(out / "supplement.pdf")
for number, page in enumerate(pdf, 1):
    for block in page.get_text("dict")["blocks"]:
        for line in block.get("lines", []):
            for span in line["spans"]:
                if span["text"].strip() and not (page.rect + (-1, -1, 1, 1)).contains(fitz.Rect(span["bbox"])):
                    raise SystemExit(f"supplement page {number}: text outside page bounds: {span['text'][:40]}")
for record in overlays:
    hits = [(i, page) for i, page in enumerate(pdf, 1) if page.search_for(record["heading"])]
    if len(hits) != 1:
        raise SystemExit(f"overlaid heading must appear on exactly one page: {record['heading']}")
    number, page = hits[0]
    head = fitz.Rect()
    for text in filter(None, (record["heading"], record["slot"])):
        for rect in page.search_for(text):
            head.include_rect(rect)
    rest = fitz.Rect()
    sizes = []
    for block in page.get_text("dict")["blocks"]:
        for line in block.get("lines", []):
            for span in line["spans"]:
                r = fitz.Rect(span["bbox"])
                if span["text"].strip() and not r.intersects(head):
                    rest.include_rect(r)
                    sizes.append(span["size"])
    for drawing in page.get_drawings():
        rest.include_rect(drawing["rect"])
    if rest.y0 <= head.y1 or min(sizes) < FONT_FLOOR:
        raise SystemExit(f"supplement page {number}: heading overlaps table content or table text < {FONT_FLOOR} pt")
    record.update(page=number, heading_bottom_bp=round(head.y1, 2), content_top_bp=round(rest.y0, 2),
                  page_min_font_pt=round(min(sizes), 2))
with (out / "HEADING_OVERLAYS.tsv").open("w", newline="") as handle:
    fields = ["page", "table", "heading", "slot", "scale", "shift_bp", "band_bp", "min_font_pt",
              "heading_bottom_bp", "content_top_bp", "page_min_font_pt"]
    writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(overlays)
print(f"SUPPLEMENT_PASS pages={len(pdf)} heading_overlays={len(overlays)}")
