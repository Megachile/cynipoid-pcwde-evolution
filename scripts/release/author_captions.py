"""Author-submitted supplement captions keyed by stable ID.

Source: authority/captions/author_captions.json (verbatim author text with its
submission provenance). Renderers look captions up by stable ID, never by label,
so renumbering cannot attach a caption to the wrong item.
"""
import json
import re
from pathlib import Path

PATH = Path(__file__).resolve().parents[2] / "authority/captions/author_captions.json"
TEX = [("\\", r"\textbackslash{}"), ("&", r"\&"), ("%", r"\%"), ("$", r"\$"), ("#", r"\#"),
       ("_", r"\_"), ("{", r"\{"), ("}", r"\}"), ("~", r"\textasciitilde{}"), ("^", r"\textasciicircum{}")]


def load(path=PATH):
    return {key: value["caption"] for key, value in json.loads(Path(path).read_text())["captions"].items()}


def to_tex(text):
    """Escape TeX specials and render *italic* markup; unbalanced markup fails."""
    for old, new in TEX:
        text = text.replace(old, new)
    text = re.sub(r"\*([^*\s][^*]*?)\*", r"\\textit{\1}", text)
    if "*" in text:
        raise ValueError(f"unbalanced italic markup in caption: {text[:80]}")
    return text


def plain(text):
    """Reader text without markup (for page-text probes)."""
    return re.sub(r"\*([^*]+)\*", r"\1", text)
