#!/usr/bin/env python3
"""Small plus/minus and phase/GC3 fixture for the frozen intron make route."""
import argparse
import tempfile
from pathlib import Path

from pyfaidx import Fasta

from build_host_background import process
from build_member_observations import build, reverse_complement


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--work', type=Path, required=True)
    a = parser.parse_args()
    a.work.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=a.work, ignore_cleanup_errors=True) as tmp:
        root = Path(tmp)
        plus = 'ATGAAAGTTTAGCCCTAA'
        minus = reverse_complement(plus)
        phase = 'AAATGAAACCCTAA'
        (root / 'tiny.fasta').write_text(f'>plus\n{plus}\n>minus\n{minus}\n>phase\n{phase}\n')
        (root / 'tiny.gff').write_text(
            'plus\tfixture\tmRNA\t1\t18\t.\t+\t.\tID=plus_tx;Parent=plus_gene\n'
            'plus\tfixture\tCDS\t1\t6\t.\t+\t0\tParent=plus_tx;protein_id=plus_pid\n'
            'plus\tfixture\tCDS\t13\t18\t.\t+\t0\tParent=plus_tx;protein_id=plus_pid\n'
            'minus\tfixture\tmRNA\t1\t18\t.\t-\t.\tID=minus_tx;Parent=minus_gene\n'
            'minus\tfixture\tCDS\t1\t6\t.\t-\t0\tParent=minus_tx;protein_id=minus_pid\n'
            'minus\tfixture\tCDS\t13\t18\t.\t-\t0\tParent=minus_tx;protein_id=minus_pid\n'
            'phase\tfixture\tmRNA\t1\t14\t.\t+\t.\tID=phase_tx;Parent=phase_gene\n'
            'phase\tfixture\tCDS\t1\t14\t.\t+\t2\tParent=phase_tx;protein_id=phase_pid\n'
            'NW_alias\tfixture\tmRNA\t1\t18\t.\t+\t.\tID=alias_tx;Parent=alias_gene\n'
            'NW_alias\tfixture\tCDS\t1\t18\t.\t+\t0\tParent=alias_tx;protein_id=alias_pid\n')
        source = dict(source_file='tiny', assembly_snapshot_path='tiny.fasta', gff_snapshot_path='tiny.gff')
        accounting = [dict(protein_id=f'{strand}_pid', source_file='tiny', family='GH5_2', genome='tiny',
                           status='single_scaffold', cds_segments='2', scaffold=strand, strand=sign)
                      for strand, sign in [('plus', '+'), ('minus', '-')]]
        segments = []
        for strand, sign in [('plus', '+'), ('minus', '-')]:
            intervals = [(1, 6), (13, 18)] if sign == '+' else [(13, 18), (1, 6)]
            for order, (start, end) in enumerate(intervals, 1):
                segments.append(dict(protein_id=f'{strand}_pid', source_file='tiny', family='GH5_2',
                                     genome='tiny', scaffold=strand, strand=sign, transcript_order=str(order),
                                     cds_start=str(start), cds_end=str(end), gff_phase='0'))
        accounting.append(dict(protein_id='phase_pid', source_file='tiny', family='GH5_2', genome='tiny',
                               status='single_scaffold', cds_segments='1', scaffold='phase', strand='+'))
        segments.append(dict(protein_id='phase_pid', source_file='tiny', family='GH5_2', genome='tiny',
                             scaffold='phase', strand='+', transcript_order='1', cds_start='1',
                             cds_end='14', gff_phase='2'))
        observations, composition, proteins = build(accounting, segments, [source], root, root / 'indexes', 3, 3)
        assert len(composition) == 3 and sum(x['intron_count'] == 1 for x in composition) == 2
        assert composition[0]['gc3'] == composition[1]['gc3']
        assert len([r for r in observations if r['intron_number']]) == 2
        assert proteins['GH5_2']['plus_pid'] == proteins['GH5_2']['minus_pid'] == 'MKP'
        assert proteins['GH5_2']['phase_pid'] == 'MKP'
        assert composition[2]['gc3'] == '' and composition[2]['gc3_phase_corrected'] != ''
        focal_ids = {'plus_pid', 'minus_pid', 'phase_pid'}
        unaliased = process(source, focal_ids, root, root / 'indexes', {})[2]
        assert unaliased['missing_sequence'] == 1 and unaliased['seqid_alias_resolved'] == 0
        alias = {'NW_alias': dict(assembly_seqid='plus', length='18')}
        background, focal, qc = process(source, focal_ids, root, root / 'indexes', alias)
        assert len(background) == 4 and len(focal) == 3
        assert qc['starting_cds_transcripts'] == qc['longest_cds_per_gene'] == 4
        assert qc['seqid_alias_resolved'] == 1 and qc['missing_sequence'] == 0
        assert qc['gc3_eligible_background'] == 3
        assert all(x['gc3'] == composition[0]['gc3'] for x in background if x['transcript_id'] != 'phase_tx')
        # The raw-CDS GC3 rule keeps a terminal stop; initial GFF phase affects
        # translation, while phase-corrected GC3 remains available as fallback.
        raw = phase
        assert raw[2:] == 'ATGAAACCCTAA'
        raw_triplets = [raw[i:i + 3] for i in range(0, len(raw) - 2, 3)]
        coding_triplets = [raw[2:][i:i + 3] for i in range(0, len(raw[2:]) - 2, 3)]
        assert raw_triplets != coding_triplets
    print('FIXTURE_PASS plus_minus=2 introns=2 phase_skip=2 gc3_fallback=1 seqid_alias=1')


if __name__ == '__main__':
    main()
