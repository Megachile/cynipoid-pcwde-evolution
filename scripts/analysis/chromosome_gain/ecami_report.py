#!/usr/bin/env python3
"""Summarize the eCAMI chromosome-gain scores without selecting a null stratum."""
import argparse
import csv
import hashlib
import json
import statistics
from pathlib import Path

FAMILIES = ('GH5_2', 'PL1_8', 'PL4_4')
MODES = ('deltran', 'acctran', 'sankoff')
STRATA = ('1', '(1,2]', '(2,4]', '>4')


def rows(path):
    with path.open(newline='') as handle:
        yield from csv.DictReader(handle, delimiter='\t')


def write(path, fields, values):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(values)


def stratum(value):
    n = float(value)
    return '1' if n == 1 else '(1,2]' if n <= 2 else '(2,4]' if n <= 4 else '>4'


def midrank(value, distribution):
    return 100 * (sum(x < value for x in distribution) +
                  .5 * sum(x == value for x in distribution)) / len(distribution)


def rate(numerator, denominator):
    return numerator / denominator if denominator else None


def fmt(value, digits=1):
    return '' if value is None else f'{value:.{digits}f}'


def main():
    p = argparse.ArgumentParser()
    for name in ('scores', 'eligibility', 'scan-eligibility', 'component-sizes', 'handchecks', 'table-numbers',
                 'figure-numbers', 'out'):
        p.add_argument('--' + name, type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    labels = {r['table_id']: r['label'] for r in rows(a.table_numbers)}
    figures = {r['figure_id']: r['label'] for r in rows(a.figure_numbers)}
    table_label = labels['chromosome_gain_summary']
    scores = list(rows(a.scores))
    native = [r for r in scores if r['kind'] == 'native']
    focal = {r['group'].removeprefix('FOCAL:'): r for r in scores if r['kind'] == 'focal'}
    scan = {r['group'].removeprefix('SCAN:'): r for r in scores if r['kind'] == 'scan'}
    scan_max = {r['group'].removeprefix('SCAN_MAX_DIAGNOSTIC:'): r for r in scores
                if r['kind'] == 'scan_max_diagnostic'}
    scan_eligibility = {r['family']: r for r in rows(a.scan_eligibility)}
    eligibility = list(rows(a.eligibility))
    if len(native) < 10 or len(focal) != 4 or len(scan) != 3 or len(scan_max) != 3 or len(scan_eligibility) != 3:
        raise SystemExit('incomplete score set')
    if sum(r['reason'] == 'eligible' for r in eligibility) != len(native):
        raise SystemExit('eligibility/score mismatch')
    handchecks = list(rows(a.handchecks))
    if len(handchecks) != 10 or sum(r['group'].startswith('NATIVE:') for r in handchecks) != 5:
        raise SystemExit('missing focal/native gain checks')
    all_rows, strata_rows, scan_rows = [], [], []
    for label in STRATA:
        subset = [r for r in native if stratum(r['median_copies']) == label]
        for mode in MODES:
            rates = [float(r[mode + '_rate']) for r in subset]
            strata_rows.append(dict(copy_stratum=label, mode=mode, native_groups=len(subset),
                                    median_gains=fmt(statistics.median(int(r[mode + '_gains']) for r in subset), 2),
                                    median_rate=fmt(statistics.median(rates), 4),
                                    zero_gain_groups=sum(int(r[mode + '_gains']) == 0 for r in subset),
                                    zero_gain_share=fmt(rate(sum(int(r[mode + '_gains']) == 0 for r in subset), len(subset)), 3)))
    for family in FAMILIES + ('B-lectin',):
        r = focal[family]
        own = [x for x in native if stratum(x['median_copies']) == stratum(r['median_copies'])]
        for mode in MODES:
            value = float(r[mode + '_rate'])
            all_rates = [float(x[mode + '_rate']) for x in native]
            own_rates = [float(x[mode + '_rate']) for x in own]
            row = dict(family=family, mode=mode, root_present=r['root_present'],
                       root_ambiguous=r['root_ambiguous'], placed_tips=r['placed_tips'],
                       gains=r[mode + '_gains'], losses=r[mode + '_losses'],
                       gains_per_placed_tip=fmt(value, 4), median_copies=r['median_copies'],
                       copy_stratum=stratum(r['median_copies']), all_native_groups=len(native),
                       all_midrank_percentile=fmt(midrank(value, all_rates)),
                       own_stratum_groups=len(own),
                       own_midrank_percentile=fmt(midrank(value, own_rates)),
                       cynipini_gains=r[mode + '_cynipini_gains'],
                       cynipini_placed_tips=r['cynipini_placed_tips'],
                       cynipini_rate=fmt(rate(int(r[mode + '_cynipini_gains']), int(r['cynipini_placed_tips'])), 4),
                       elsewhere_gains=r[mode + '_elsewhere_gains'],
                       elsewhere_placed_tips=r['elsewhere_placed_tips'],
                       elsewhere_rate=fmt(rate(int(r[mode + '_elsewhere_gains']), int(r['elsewhere_placed_tips'])), 4))
            if family == 'B-lectin':
                row['all_midrank_percentile'] = row['own_midrank_percentile'] = ''
            all_rows.append(row)
        if family in scan:
            q = scan[family]
            e = scan_eligibility[family]
            scan_rows.append(dict(family=family, ecami_group=e['ecami_group'],
                                  member_genomes=e['member_genomes'], crown_eligible=e['crown_eligible'],
                                  placed_tips=q['placed_tips'],
                                  deltran_gains=q['deltran_gains'],
                                  acctran_gains=q['acctran_gains'], sankoff_gains=q['sankoff_gains'],
                                  **{mode + '_rate': q[mode + '_rate'] for mode in MODES},
                                  **{mode + '_all_midrank_percentile': fmt(midrank(float(q[mode + '_rate']),
                                      [float(x[mode + '_rate']) for x in native])) if e['crown_eligible'] == '1' else ''
                                      for mode in MODES},
                                  **{mode + '_own_midrank_percentile': fmt(midrank(float(q[mode + '_rate']),
                                      [float(x[mode + '_rate']) for x in native if stratum(x['median_copies']) == stratum(q['median_copies'])]))
                                      if e['crown_eligible'] == '1' else '' for mode in MODES},
                                  median_copies=q['median_copies'], copy_stratum=stratum(q['median_copies']),
                                  root_present=q['root_present'],
                                  root_ambiguous=q['root_ambiguous']))
    if any((float(row['all_midrank_percentile']) - 50) *
           (float(next(x for x in all_rows if x['family'] == row['family'] and x['mode'] == 'deltran')['all_midrank_percentile']) - 50) < 0
           for row in all_rows if row['family'] in FAMILIES):
        raise SystemExit('parsimony modes reverse the all-native comparison')
    reader = []
    for r in (row for row in all_rows if row['mode'] == 'deltran'):
        reader.append(dict(section='A', family_or_stratum=r['family'], mode=r['mode'],
                           gains=r['gains'], cynipini_gains=r['cynipini_gains'],
                           other_lineage_gains=r['elsewhere_gains'], losses=r['losses'],
                           placed_tips_or_native_groups=r['placed_tips'],
                           gains_per_placed_tip=r['gains_per_placed_tip'],
                           all_native_percentile=r['all_midrank_percentile'],
                           own_stratum_percentile=r['own_midrank_percentile'],
                           zero_gain_share=''))
    for r in (row for row in strata_rows if row['mode'] == 'deltran'):
        reader.append(dict(section='B', family_or_stratum=r['copy_stratum'], mode=r['mode'],
                           gains=r['median_gains'], losses='',
                           placed_tips_or_native_groups=r['native_groups'],
                           gains_per_placed_tip=r['median_rate'],
                           all_native_percentile='', own_stratum_percentile='',
                           zero_gain_share=r['zero_gain_share']))
    write(a.out / 'Table_S15.tsv', list(reader[0]), reader)
    (a.out / 'Table_S15_caption.txt').write_text(
        f'Table {table_label}. DELTRAN-inferred gains of family presence on BK reference chromosomes for the focal PCWDE families '
        f'and {len(native)} native eCAMI groups present at the cynipoid crown. Section A '
        'gives focal gains, their Cynipini versus other-lineage split, losses and normalized rates, '
        'with descriptive midrank percentiles against all native '
        'groups and each family’s own median-copy stratum. B-lectin is Cynipidae-only '
        'and has no native percentile. Section B gives native copy-stratum medians; '
        'its gains column is the median count, its denominator column is the number '
        'of native groups, and groups with no inferred gains are shown as a percentage. eCAMI groups are '
        'dbCAN-sub sequence clusters, not CAZy-curated subfamilies. Retained gains '
        'reflect event occurrence, retention and elapsed time; the strata are not '
        'adjusted null models.\n')
    write(a.out / 'focal_summary.tsv', list(all_rows[0]), [r for r in all_rows if r['mode'] == 'deltran'])
    write(a.out / 'copy_strata.tsv', list(strata_rows[0]), [r for r in strata_rows if r['mode'] == 'deltran'])
    write(a.out / 'focal_scan_sensitivity.tsv', list(scan_rows[0]), scan_rows)
    write(a.out / 'focal_scan_max_diagnostic.tsv',
          ['family', 'placed_tips', 'deltran_gains', 'acctran_gains', 'sankoff_gains'],
          [dict(family=family, placed_tips=scan_max[family]['placed_tips'],
                **{mode + '_gains': scan_max[family][mode + '_gains'] for mode in MODES})
           for family in FAMILIES])
    core = [r for r in native if r['single_copy_core'] == '1']
    if not core:
        raise SystemExit('no single-copy noise floor')
    floor = dict(single_copy_groups=len(core),
                 zero_gain_groups=sum(int(r['deltran_gains']) == 0 for r in core),
                 median_gains=statistics.median(int(r['deltran_gains']) for r in core),
                 median_rate=statistics.median(float(r['deltran_rate']) for r in core))
    (a.out / 'noise_floor.json').write_text(json.dumps(floor, indent=2) + '\n')
    native_clade = {}
    for clade in ('cynipini', 'elsewhere'):
        values = [int(r[clade + '_gains']) / int(r[clade + '_placed_tips'])
                  for r in native if int(r[clade + '_placed_tips'])]
        native_clade[clade] = dict(groups=len(values), median_rate=statistics.median(values),
                                   zero_gain_share=sum(v == 0 for v in values) / len(values),
                                   pooled_gains=sum(int(r[clade + '_gains']) for r in native),
                                   pooled_placed_tips=sum(int(r[clade + '_placed_tips']) for r in native))
        native_clade[clade]['pooled_rate'] = native_clade[clade]['pooled_gains'] / native_clade[clade]['pooled_placed_tips']
    (a.out / 'native_clade_rates.json').write_text(json.dumps(native_clade, indent=2) + '\n')
    primary_rows = [r for r in all_rows if r['family'] in FAMILIES]
    all_low, all_high = min(float(r['all_midrank_percentile']) for r in primary_rows), max(float(r['all_midrank_percentile']) for r in primary_rows)
    own_deltran = [float(r['own_midrank_percentile']) for r in primary_rows if r['mode'] == 'deltran']
    own_n = int(next(r for r in primary_rows if r['mode'] == 'deltran')['own_stratum_groups'])
    lines = ['# Retained chromosome gains against native eCAMI groups', '',
             f'Across DELTRAN, ACCTRAN and Sankoff, the three focal families rank at percentiles {all_low:.1f}–{all_high:.1f} for retained BK-chromosome gains per placed genome among {len(native)} eligible native eCAMI groups. In their own median-copy stratum (n={own_n}), their DELTRAN midranks are at percentiles {min(own_deltran):.1f}–{max(own_deltran):.1f}. These are descriptive comparisons; copy number is not adjusted away.', '',
             f'The separate historical 40-tip regression reproduces the accepted-copy scores. On the expanded 48-tip tree, PL1_8 roots at {focal["PL1_8"]["root_present"]} and PL4_4 at {focal["PL4_4"]["root_present"]}; GH5_2 has root {focal["GH5_2"]["root_present"]} and ambiguous {focal["GH5_2"]["root_ambiguous"]}. The native background has {sum(int(r["deltran_gains"]) == 0 for r in native)}/{len(native)} zero-DELTRAN-gain groups. The single-copy noise floor has {floor["zero_gain_groups"]}/{floor["single_copy_groups"]} groups with zero DELTRAN gains (median rate {floor["median_rate"]:.3f}).', '',
             f'Table {table_label} carries all modes, rates, midranks, lineage splits and copy strata; focal_summary.tsv carries the root states and clade counts. eCAMI groups are dbCAN-sub sequence clusters, not CAZy-curated subfamilies. Retained gains reflect event occurrence, retention and elapsed time; more copies provide more opportunities for retention. The PCWDE families began from one copy. No transfer direction or rearrangement mechanism is inferred.', '',
             '## Per-family results', '',
             '| Family | Mode | Gains / losses | Gains per placed tip | All native | Own stratum |',
             '| --- | --- | --- | --- | --- | --- |']
    for family in FAMILIES:
        for mode in MODES:
            r = next(x for x in all_rows if x['family'] == family and x['mode'] == mode)
            lines.append(f'| {family} | {mode.upper()} | {r["gains"]} / {r["losses"]} | {r["gains_per_placed_tip"]} | {r["all_midrank_percentile"]}% | {r["own_midrank_percentile"]}% |')
    lines += ['', '## Fixed eCAMI component sensitivity', '',
              'Each family uses one gene-collapsed component group selected across all 48 genomes before scoring. The historical genome-wise maximum trajectory is retained only as a diagnostic.', '',
              '| Family | Fixed group | Member / placed genomes | Crown eligible | DELTRAN / ACCTRAN / Sankoff gains |',
              '| --- | --- | --- | --- | --- |']
    for r in scan_rows:
        lines.append(f'| {r["family"]} | {r["ecami_group"]} | {r["member_genomes"]} / {r["placed_tips"]} | {r["crown_eligible"]} | {r["deltran_gains"]} / {r["acctran_gains"]} / {r["sankoff_gains"]} |')
    lines += ['', 'DELTRAN gains by clade: ' + '; '.join(
        f'{family} Cynipini {next(r for r in all_rows if r["family"] == family and r["mode"] == "deltran")["cynipini_gains"]}/{next(r for r in all_rows if r["family"] == family and r["mode"] == "deltran")["cynipini_placed_tips"]}, elsewhere {next(r for r in all_rows if r["family"] == family and r["mode"] == "deltran")["elsewhere_gains"]}/{next(r for r in all_rows if r["family"] == family and r["mode"] == "deltran")["elsewhere_placed_tips"]}'
        for family in FAMILIES), '']
    lines += ['', f'Native pooled clade rates: Cynipini {native_clade["cynipini"]["pooled_gains"]}/{native_clade["cynipini"]["pooled_placed_tips"]} = {native_clade["cynipini"]["pooled_rate"]:.3f} gains per placed tip; elsewhere {native_clade["elsewhere"]["pooled_gains"]}/{native_clade["elsewhere"]["pooled_placed_tips"]} = {native_clade["elsewhere"]["pooled_rate"]:.3f}. Group-median rates are {native_clade["cynipini"]["median_rate"]:.3f} and {native_clade["elsewhere"]["median_rate"]:.3f}, respectively. Clade rates are descriptive and groups differ in placed-tip coverage.', '',
              'B-lectin is Cynipidae-only and is shown descriptively in the table, without a native percentile.', '']
    (a.out / 'summary.md').write_text('\n'.join(lines))
    def digest(path):
        return hashlib.sha256(path.read_bytes()).hexdigest()
    (a.out / 'ecami_report_sha256.json').write_text(json.dumps({
        'inputs': {name: digest(getattr(a, name)) for name in
                   ('scores', 'eligibility', 'scan_eligibility', 'component_sizes',
                    'handchecks', 'table_numbers', 'figure_numbers')},
        'outputs': {name: digest(a.out / name) for name in
                    ('Table_S15.tsv', 'Table_S15_caption.txt', 'focal_summary.tsv')},
    }, indent=2, sort_keys=True) + '\n')
    print(f'ECAMI_REPORT_PASS native={len(native)} table={table_label}')


if __name__ == '__main__':
    main()
