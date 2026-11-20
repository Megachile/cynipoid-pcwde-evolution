#!/usr/bin/env python3
"""Read-only re-chaining of retained tBLASTn HSPs (job 1424546) for all five
serial-search queries. Same collinearity/<=100 kb gap/unique-global-best and
>=95% identity / >=90% query coverage bar as recover_modelless_batched.py; the
only change is that adjacent collinear HSPs may overlap by <=30 query residues
(<=90 bp genomic), which the strict nonoverlap rule rejected. Coverage is the
union of query intervals; identity is summed nident / summed alignment length.
No new search. Original BATCHED_COORDINATE_RECOVERY.tsv is left unchanged.
All ten inputs are checked against RECOVERY_PINS.tsv before use."""
import csv, sys
from collections import defaultdict
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import recovery_pins
O = Path(sys.argv[1])
QOV, GOV, GAP = 30, 90, 100000
IDS = ['I_centaureae_beta_lectin_0001', 'Qwaqwaia_scolopiae_beta_lectin_0001',
       'Neu_valhalla_beta_lectin_0004', 'Me_ostensackeni_beta_lectin_0002',
       'DK_beta_lectin_0006']
V = recovery_pins.Verifier()
V.require([(recovery_pins.RECOVERY_PREFIX + pid + ext, O / (pid + ext)) for pid in IDS for ext in ('.faa', '.tblastn.tsv')])
def qlen(p):
    return sum(len(l.strip()) for l in V.check(p).open() if not l.startswith('>'))
def union(iv):
    tot, cur = 0, None
    for a, b in sorted(iv):
        if cur is None or a > cur[1] + 1:
            if cur: tot += cur[1] - cur[0] + 1
            cur = [a, b]
        else: cur[1] = max(cur[1], b)
    return tot + (cur[1] - cur[0] + 1 if cur else 0)
out = []
for pid in IDS:
    L = qlen(O / f'{pid}.faa')
    hs = set()
    for r in csv.reader(V.check(O / f'{pid}.tblastn.tsv').open(), delimiter='\t'):
        if r[0] == 'sseqid': continue
        hs.add((r[0], int(r[1]), int(r[2]), int(r[3]), int(r[4]), int(r[5]), int(r[6]), float(r[7])))
    groups = defaultdict(list)
    for h in hs: groups[(h[0], 1 if h[4] >= h[3] else -1)].append(h)
    cands = []
    for (sc, st), g in groups.items():
        g.sort(key=lambda h: (h[1], h[2], h[3]))
        dp = []
        for i, h in enumerate(g):
            best = (h[7], [h])
            for j, p in enumerate(g[:i]):
                gap = (h[3] - p[4]) * st          # >0 forward progress
                if h[1] > p[1] and h[2] > p[2] and p[2] - h[1] < QOV and -GOV <= gap <= GAP:
                    opt = (dp[j][0] + h[7], dp[j][1] + [h])
                    if opt[0] > best[0]: best = opt
            dp.append(best)
        for sc_, ch in dp:
            cands.append(dict(scaffold=sc, start=min(min(h[3:5]) for h in ch), end=max(max(h[3:5]) for h in ch),
                              strand='+' if st == 1 else '-', identity=sum(h[5] for h in ch) / sum(h[6] for h in ch),
                              query_coverage=union([(h[1], h[2]) for h in ch]) / L, score=sc_, hsps=len(ch),
                              chain=';'.join(f'{h[1]}-{h[2]}@{h[3]}-{h[4]}' for h in ch)))
    cands.sort(key=lambda c: (-c['score'], c['scaffold'], c['start']))
    top = cands[0]
    ties = [c for c in cands if c['score'] == top['score'] and (c['scaffold'], c['start'], c['end']) != (top['scaffold'], top['start'], top['end'])]
    # best competing chain at a non-overlapping locus
    comp = next((c for c in cands if c['scaffold'] != top['scaffold'] or c['end'] < top['start'] or c['start'] > top['end']), None)
    ok = not ties and top['identity'] >= .95 and top['query_coverage'] >= .90
    out.append(dict(accepted_id=pid, query_len=L, scaffold=top['scaffold'], start=top['start'], end=top['end'], strand=top['strand'],
                    identity=round(top['identity'], 5), query_coverage=round(top['query_coverage'], 5), score=top['score'],
                    hsps=top['hsps'], chain=top['chain'], n_tied_other_loci=len(ties),
                    best_competing_locus=f"{comp['scaffold']}:{comp['start']}-{comp['end']}" if comp else '',
                    best_competing_score=comp['score'] if comp else '',
                    status='resolved' if ok else 'unresolved_tie_or_below_threshold'))
w = csv.DictWriter((O / 'OVERLAP_TOLERANT_RECHAIN.tsv').open('w'), fieldnames=list(out[0]), delimiter='\t')
w.writeheader(); w.writerows(out)
for r in out: print(r['accepted_id'], r['status'], r['scaffold'], r['start'], r['end'], r['identity'], r['query_coverage'], r['best_competing_score'])
