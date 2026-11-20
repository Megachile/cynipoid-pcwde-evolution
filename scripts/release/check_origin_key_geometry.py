"""Check actual PDF key frames, swatches and label bounds in origin artwork."""
import argparse
import json
from pathlib import Path
import fitz

TAXA = {'Cynipoid wasps', 'Midge (Diptera)', 'Bee', 'Beetle', 'Nematode',
        'Oomycete', 'Fungus', 'Plant', 'Protist'}
DOMAINS = set('GH5_2 CBM2 CBM3 CBM5 CBM6 CBM10 CBM12 CBM17 CBM28 CBM73 CBM81 PL1_8 CBM13 PL4 PL4_2 PL4_4 PL4_5 CBM22 CBM32 CBM64 GH30 GH30_C CBM4_9 fn3 Ricin_B Big_5 GH59 GH30_2 RicinB_2 CBM35_2 SLH PI-PLC-X'.split()) | {'CBM13 (weak)', 'Region A only'}

def inspect(path):
    doc = fitz.open(path)
    assert len(doc) == 1, path
    p = doc[0]
    spans = [(s['text'].strip(), fitz.Rect(s['bbox']))
             for b in p.get_text('dict')['blocks'] for l in b.get('lines', []) for s in l['spans']]
    taxa = [(t, r) for t, r in spans if t in TAXA]
    domains = [(t, r) for t, r in spans if t in DOMAINS]
    errors = []
    if not taxa or not domains:
        errors.append('missing taxon or domain labels')
    drawings = []
    clips = {}
    for drawing in p.get_drawings(extended=True):
        level = drawing.get('level', 0)
        clips = {k: v for k, v in clips.items() if k < level}
        if drawing['type'] == 'clip':
            clips[level] = drawing['scissor']
        else:
            drawing['active_clip'] = list(clips.values())
            drawings.append(drawing)
    frames = [d['rect'] for d in drawings if d['fill'] is None and d['rect'].width > 1 and d['rect'].height > 1
              and any(d['rect'].contains(r) for _, r in taxa)]
    if len(frames) != 1:
        errors.append(f'expected one taxon frame, found {len(frames)}')
    boxes = []
    for kind, labels in [('taxon', taxa), ('domain', domains)]:
        for label, r in labels:
            boxes.append((kind, label, r))
            # Associate the nearest filled rectangular swatch to the left on this row.
            choices = [d for d in drawings if d['fill'] is not None
                       and d['rect'].width > 1 and d['rect'].height > 1
                       and d['rect'].x1 <= r.x0 + .2
                       and abs((d['rect'].y0+d['rect'].y1-r.y0-r.y1)/2) < r.height/2]
            if not choices:
                errors.append(f'missing swatch: {label}')
            else:
                drawing = max(choices, key=lambda x: x['rect'].x1)
                swatch = drawing['rect']
                for clip in drawing['active_clip']:
                    if not clip.contains(r) or not clip.contains(swatch):
                        errors.append(f'clipped at panel boundary: {label}')
                boxes.append((kind, label+' swatch', swatch))
    for kind, label, r in boxes:
        if not p.rect.contains(r):
            errors.append(f'outside page: {label}')
        if kind == 'taxon' and frames and not frames[0].contains(r):
            errors.append(f'outside taxon frame: {label}')
    for frame in frames:
        for kind, label, r in boxes:
            if kind == 'domain' and frame.intersects(r):
                errors.append(f'domain intersects taxon frame: {label}')
    for i, (kind, label, r) in enumerate(boxes):
        for kind2, label2, r2 in boxes[i+1:]:
            if r.intersects(r2):
                errors.append(f'key collision: {label} / {label2}')
    return {'file': path.name, 'page': list(p.rect), 'taxon_frames': [list(r) for r in frames],
            'boxes': [{'kind': k, 'label': t, 'bbox': list(r)} for k,t,r in boxes], 'errors': errors}

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--figures', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    names = [f'origin_trees_full_{f}_ggtree.pdf' for f in ('GH5_2','PL1_8','PL4_4')] + ['origin_trees.pdf']
    reports = [inspect(args.figures/n) for n in names]
    args.out.write_text(json.dumps(reports, indent=2)+'\n')
    for report in reports:
        print(report['file'], 'FAIL '+str(report['errors']) if report['errors'] else 'PASS')
    raise SystemExit(any(r['errors'] for r in reports))
