#!/usr/bin/env python3
"""Render the native eCAMI distribution and focal retained-gain trees."""
import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from render import COLORS


def rows(path):
    with path.open(newline='') as handle:
        yield from csv.DictReader(handle, delimiter='\t')


def main():
    p = argparse.ArgumentParser()
    for name in ('scores', 'tips', 'changes', 'tree', 'figure-numbers', 'out'):
        p.add_argument('--' + name, type=Path, required=True)
    a = p.parse_args()
    a.out.parent.mkdir(parents=True, exist_ok=True)
    labels = {r['figure_id']: r['label'] for r in rows(a.figure_numbers)}
    label = labels['chromosome_gain_distribution']
    scores = list(rows(a.scores))
    native = [float(r['deltran_rate']) for r in scores if r['kind'] == 'native']
    focal = {r['group'].removeprefix('FOCAL:'): float(r['deltran_rate'])
             for r in scores if r['kind'] == 'focal'}
    if len(native) < 10 or not all(f in focal for f in COLORS):
        raise SystemExit('incomplete figure scores')
    with PdfPages(a.out) as pdf:
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.hist(native, bins=30, color='#87979f', edgecolor='white')
        for family, color in COLORS.items():
            ax.axvline(focal[family], color=color, linewidth=2.4, label=family)
        ax.set_xlabel('Retained BK-chromosome gains per placed genome')
        ax.set_ylabel('Eligible native eCAMI groups')
        ax.set_title(f'Figure {label}. Native eCAMI chromosome gains and focal families', loc='left')
        ax.legend(frameon=False)
        fig.text(.08, .05, 'DELTRAN reconstruction; eCAMI groups are dbCAN-sub sequence clusters.', fontsize=9)
        fig.tight_layout(rect=(0, .07, 1, 1))
        pdf.savefig(fig)
        plt.close(fig)
    print(f'ECAMI_FIGURE_PASS {label} pages=1 native={len(native)}')


if __name__ == '__main__':
    main()
