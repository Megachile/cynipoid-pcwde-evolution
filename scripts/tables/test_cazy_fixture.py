#!/usr/bin/env python3
"""Small end-to-end fixture for the CAZy scan and rank producers."""

import csv
import hashlib
import shlex
import subprocess
import sys
import tempfile
from pathlib import Path

from cazy_subfamily_rule import assign
from build_cazy_ranks import rank
from cazy_gene_rule import map_genes
from reconcile_cazy_members import relation
from cazy_data_root import resolve_input


def write(path, columns, records):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='') as handle:
        w = csv.writer(handle, delimiter='\t', lineterminator='\n')
        w.writerow(columns)
        w.writerows(records)


def read(path):
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(script, *args):
    subprocess.run([sys.executable, '-B', str(Path(__file__).with_name(script)), *map(str, args)], check=True)


def main():
    assert assign('GH18_e1', '-') == ('GH18', 0.0)
    assert assign('GH5_e7', 'GH5_2:7;GH5:1')[0] == 'GH5_2'
    assert assign('GH5_e7', 'GH5_2:2;GH5_3:2')[0] == 'GH5'
    assert assign('GH5_e7', 'GH5:7;GH5_2:1')[0] == 'GH5'
    assert rank(2, {'GH18': 2, 'GT1_3': 2}, 2)[0:2] == (1, '2.00')
    assert len(set(map_genes('dbcan_duplicate_suffix', None,
                            {'AE_bkragtag_pseudo_chr_1_000001.1',
                             'AE_bkragtag_pseudo_chr_1_000001.1_1'}).values())) == 1
    assert relation('KTA', 'MKTAAAA') == ('contiguous_accepted_fragment', '2', '4')
    assert relation('MK.A', 'MKXA')[0] == 'normalization_only'
    assert relation('MKK', 'MKT')[0] == 'other_sequence_difference'
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        try:
            resolve_input(root, '../cazy_fixture_outside')
            raise AssertionError('outside-root path accepted')
        except ValueError:
            pass
        with tempfile.NamedTemporaryFile(dir=root.parent) as external:
            (root / 'escaping_link').symlink_to(external.name)
            try:
                resolve_input(root, 'escaping_link')
                raise AssertionError('escaping symlink accepted')
            except ValueError:
                pass
        manifest = []
        for i in range(64):
            sp = f'Species_{i:02d}'
            hits = [('H1.1', 'GH1_e72', 'GH1:68'), ('H1.2', 'GH1_e72', 'GH1:68')]
            if i == 0:
                hits += [('X.1', 'GH5_e7', 'GH5_2:7;GH5:1'),
                         ('Y.1', 'GH5_e7', 'GH5_2:7;GH5:1'),
                         ('P.1', 'PL1_e8', 'PL1_8:5'),
                         ('P.2', 'PL1_e106', 'PL1_8:6'),
                         ('N.1', 'GH18_e1', '-'),
                         ('T.1', 'GT1_e1', 'GT1_3:2'),
                         ('T.2', 'GT1_e2', 'GT1_3:2')]
                hits += [('A.1', 'GH13_e1', '-'),
                         ('A.1', 'GH13_e17', 'GH13_17:8'),
                         ('B.1', 'GH18_e1', '-'),
                         ('B.1', 'CBM14_e1', '-'),
                         ('C.1', 'GH13_e17', 'GH13_17:8'),
                         ('C.1', 'GH13_e20', 'GH13_20:8')]
            scan = root / f'{sp}.tsv'
            scored = [(p, g, c, ('1e-20' if g == 'GH13_e20' else '1e-10')) for p,g,c in hits]
            write(scan, ('Target Name', 'Subfam Name', 'Subfam Composition', 'i-Evalue'), scored)
            proteome = root / f'{sp}.faa'
            proteome.write_text(''.join(f'>{p}\nMKT' +
                                        {'X.1': 'AAAA', 'Y.1': 'BBBB'}.get(p, 'CCCC') + '\n'
                                        for p, _, _ in hits))
            manifest.append((sp, 'test', sp, proteome.name, sha(proteome), scan.name,
                             sha(scan), 'frozen_model_stem', '', '', 1))
        manifest_file = root / 'manifest.tsv'
        write(manifest_file, ('species', 'source', 'source_key', 'proteome_path',
              'proteome_sha256', 'scan_path', 'scan_sha256', 'mapping_method',
              'mapping_source', 'mapping_sha256', 'expected_GH1_e72'), manifest)
        fig1 = root / 'fig1'
        write(fig1 / 'gh1e72_gene_counts.tsv', ('species', 'GH1_e72'),
              ((f'Species_{i:02d}', 1) for i in range(64)))
        write(fig1 / 'internal_gh1e72_gene_counts.tsv', ('species', 'GH1_e72'), [])
        write(fig1 / 'gh1e72_cynipoid_counts.tsv',
              ('accession', 'gh1e72_locus_count'), [])
        counts = root / 'counts'
        run('build_cazy_scan_counts.py', '--manifest', manifest_file,
            '--source-root', root, '--fig1-counts', fig1, '--out', counts)
        units = {(r['species'], r['rank_unit']): int(r['gene_count'])
                 for r in read(counts / 'unit_counts.tsv')}
        assert units['Species_00', 'GH1'] == 1  # isoforms collapse
        assert units['Species_00', 'GH5_2'] == 2  # distinct adjacent genes remain
        assert units['Species_00', 'PL1_8'] == 1  # component groups pool by gene
        assert units['Species_00', 'GT1_3'] == 1
        assert units['Species_00', 'GH13_17'] == 1
        assert units['Species_00', 'GH13_20'] == 1
        assert ('Species_00', 'GH13') not in units
        assert units['Species_00', 'GH18'] == 2
        assert units['Species_00', 'CBM14'] == 1
        assert read(counts / 'competing_subfamilies.tsv')[0]['selected_subfamily'] == 'GH13_20'
        assert len(read(counts / 'derived_gh1e72_gene_counts.tsv')) == 64
        annotation = root / 'annotation.tsv'
        write(annotation, ('protein_id', 'gene'),
              ((p, p.split('.')[0]) for p in {h[0] for h in hits}))
        tampered = [list(row) for row in manifest]
        tampered[0][7] = 'canonical_annotation'
        tampered[0][8] = annotation.name
        tampered[0][9] = sha(annotation)
        bad_manifest = root / 'tampered_manifest.tsv'
        write(bad_manifest, ('species', 'source', 'source_key', 'proteome_path',
              'proteome_sha256', 'scan_path', 'scan_sha256', 'mapping_method',
              'mapping_source', 'mapping_sha256', 'expected_GH1_e72'), tampered)
        annotation.write_text(annotation.read_text() + 'tampered\tgene\n')
        failed = subprocess.run([sys.executable, '-B', str(Path(__file__).with_name('build_cazy_scan_counts.py')),
                                 '--manifest', str(bad_manifest), '--source-root', str(root),
                                 '--fig1-counts', str(fig1), '--out', str(root / 'bad_counts')],
                                capture_output=True, text=True)
        assert failed.returncode and 'mapping source hash mismatch' in failed.stderr
        members = root / 'members.tsv'
        write(members, ('species', 'family', 'protein_id'),
              [('Species 00', 'GH5_2', 'X.1'), ('Species 00', 'GH5_2', 'Y.1'),
               ('Species 00', 'PL1_8', 'P.1'), ('Species 00', 'PL4_4', 'Z.1'),
               ('Species 00', 'GH5_2', 'N.1')])
        processed_map = root / 'processed_map.tsv'
        write(processed_map, ('species', 'processed_path', 'processed_sha256'),
              ((r[0], r[3], r[4]) for r in manifest))
        seqdir = root / 'accepted_sequences'
        seqdir.mkdir()
        write_seq = {'GH5_2': [('X.1', 'MKTAAAA'), ('Y.1', 'MKTBBBB'),
                               ('N.1', 'MKTCCCC')],
                     'PL1_8': [('P.1', 'MKTCCCC')], 'PL4_4': []}
        for fam, entries in write_seq.items():
            (seqdir / f'complete_locus_{fam}.faa').write_text(
                ''.join(f'>{ident}\n{seq}\n' for ident, seq in entries))
        crosswalk = root / 'crosswalk'
        run('reconcile_cazy_members.py', '--manifest', manifest_file,
            '--processed-map', processed_map, '--member-table', members,
            '--sequences-dir', seqdir, '--source-root', root, '--out', crosswalk)
        cross = {r['accepted_protein_id']: r for r in read(crosswalk / 'member_crosswalk.tsv')}
        assert cross['X.1']['outcome'] == 'exact_id_same_sequence'
        assert cross['Z.1']['outcome'] == 'accepted_sequence_unavailable'
        assert cross['N.1']['dbcan_sub_hit'] == 'yes'
        ranks = root / 'ranks'
        run('build_cazy_ranks.py', '--scan-counts', counts, '--member-table', members,
            '--member-crosswalk', crosswalk / 'member_crosswalk.tsv',
            '--manifest', manifest_file, '--out', ranks)
        paired = {(r['species'], r['unit'], r['route']): r
                  for r in read(ranks / 'paired_ranks.tsv')}
        assert paired['Species_00', 'GH5_2', 'accepted']['count'] == '3'
        assert paired['Species_00', 'GH5_2', 'scan']['count'] == '2'
        assert paired['Species_00', 'PL4_4', 'accepted']['count'] == '1'
        assert paired['Species_00', 'PL4_4', 'scan']['count'] == '0'
        assert paired['Species_00', 'PL4_4', 'scan']['rank_denominator'] == paired['Species_00', 'PL4_4', 'accepted']['rank_denominator']
        component = {(r['species'], r['ecami_group']): r
                     for r in read(ranks / 'component_ranks.tsv')}
        assert int(component['Species_00', 'GH5_e7']['rank_denominator']) > int(paired['Species_00', 'GH5_2', 'scan']['rank_denominator'])
        assert len(read(ranks / 'accepted_without_scan_hit.tsv')) == 1
        assert len(read(ranks / 'accepted_different_scan_unit.tsv')) == 1
        deposit = root / 'deposit.tsv'
        run('build_cazy_deposit.py', '--manifest', manifest_file,
            '--unit-counts', counts / 'unit_counts.tsv',
            '--paired-ranks', ranks / 'paired_ranks.tsv', '--out', deposit)
        dense = read(deposit)
        assert len(dense) == 64 * len({r['rank_unit'] for r in dense})
        assert next(r for r in dense if r['species'] == 'Species_00' and r['rank_unit'] == 'PL4_4')['accepted_focal_gene_count'] == '1'
        lineages = root / 'lineages.tsv'
        write(lineages, ('species', 'lineage_group'),
              ((f'Species_{i:02d}', ('E. primus' if i == 0 else
                'Diplolepididae' if i == 1 else 'Cynipini' if i == 2 else
                'non-Cynipini Cynipidae')) for i in range(64)))
        summary_dir = root / 'summary'
        recipe = next(line for line in (Path(__file__).resolve().parents[2] / 'Makefile').read_text().splitlines()
                      if 'scripts/tables/build_cazy_lineage_summary.py ' in line)
        production_args = shlex.split(recipe)
        precision_arg = production_args[production_args.index('--primary-precision') + 1]
        assert precision_arg == '$(BUILD)/cazy_ranks/primary_ecami_precision.tsv'
        precision_input = root / precision_arg.replace('$(BUILD)/', '').replace('cazy_ranks/', 'ranks/')
        assert precision_input.is_file(), precision_input
        print('CAZY_PRODUCTION_ARG --primary-precision ' + precision_arg, flush=True)
        run('build_cazy_lineage_summary.py', '--lineages', lineages,
            '--paired-ranks', ranks / 'paired_ranks.tsv',
            '--native-counts', ranks / 'native_counts.tsv',
            '--component-ranks', ranks / 'component_ranks.tsv',
            '--primary-ranks', ranks / 'primary_ecami_ranks.tsv',
            '--native-components', ranks / 'native_component_counts.tsv',
            '--group-mapping', counts / 'group_mapping.tsv',
            '--primary-precision', precision_input, '--out', summary_dir)
        assert len(read(summary_dir / 'lineage_summary.tsv')) == 4 * 3 * 3
        paper = root / 'paper_00.faa'
        paper.write_text('>Xnew\nMKTAAAA\n>Ynew\nMKTDPPP\n')
        sources = root / 'paper_sources.tsv'
        write(sources, ('species', 'proteome_path', 'proteome_sha256', 'source_role'),
              [('Species_00', paper.name, sha(paper), 'fixture')])
        diagnostic_members = root / 'diagnostic_members.tsv'
        write(diagnostic_members, ('species', 'family', 'protein_id'),
              [('Species 00', 'GH5_2', 'Xnew'), ('Species 00', 'GH5_2', 'Ynew')])
        diag = root / 'diagnostic'
        run('diagnose_cazy_inputs.py', '--manifest', manifest_file,
            '--paper-sources', sources, '--member-table', diagnostic_members,
            '--source-root', root, '--out', diag)
        outcomes = {r['accepted_protein_id']: r['outcome']
                    for r in read(diag / 'member_reconciliation.tsv')}
        assert outcomes == {'Xnew': 'sequence_rename_one_to_one',
                            'Ynew': 'changed_model_or_missing_scan'}
        print('CAZY_FIXTURE_PASS rows=64 isoforms=1 adjacent=2 zero_hit=1 mismap=1')


if __name__ == '__main__':
    main()
