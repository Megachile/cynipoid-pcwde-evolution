"""Project approved numeric tokens while preserving the pinned author words."""
import json
from pathlib import Path
from supplement_template import render, NUMERAL, WORD_NUMERAL

def project(text, numbers, document):
    audit = []
    routes = json.loads((Path(__file__).resolve().parents[2] / 'config/author_result_routes.json').read_text())
    for row in routes:
        if row['document'] != document:
            continue
        if text.count(row['author']) != 1:
            raise RuntimeError('author numeric route source drift: ' + row['author'][:70])
        value, tokens, errors = render(row['template'], numbers)
        if errors:
            raise RuntimeError(str(errors))
        mask = lambda s: NUMERAL.sub('#', WORD_NUMERAL.sub('#', s))
        if mask(row['author']) != mask(value):
            raise RuntimeError('author numeric route changed words or punctuation')
        text = text.replace(row['author'], value)
        audit.extend(dict(document=document, key=k, format=f, value=v) for k,f,v in tokens)
    return text, audit
