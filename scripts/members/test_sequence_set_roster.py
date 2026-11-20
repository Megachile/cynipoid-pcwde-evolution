#!/usr/bin/env python3
"""Small resolver fixture covering union links and symmetric collision handling."""
from sequence_set_roster import Resolver, unique_reciprocal_overlaps


def fixture(reverse=False):
    r = Resolver()
    current = r.add_row('GH5_2', 'AE_bkragtag_pseudo_chr_10_001131.1', link_status='current_copy')
    r.alias(current, 'cynipoid|Acraspis_erinacei_1')
    trinity = r.add_row('GH5_2', 'AF|AF_TRINITY_DN1015_c2_g1_i1.p1', link_status='locus_tip')
    r.alias(trinity, 'cynipoid|Andricus_quercusfoliatus_1')
    r.alias(trinity, 'AF_TRINITY_DN1015_c2_g1_i1.p1')
    ena = r.add_row('GH5_2', 'AHIE_ENA|CDXCXH010007311|CDXCXH010007311.1_000002.1', link_status='current_copy')
    ambiguous_a = r.add_row('GH5_2', 'candidate_a', link_status='current_copy')
    ambiguous_b = r.add_row('GH5_2', 'candidate_b', link_status='current_copy')
    r.alias(ambiguous_a, 'ambiguous_alias')
    r.alias(ambiguous_b, 'ambiguous_alias')
    collision = r.add_row('GH5_2', 'collision_target', link_status='current_copy',
                          public_species_name='Pediaspis aceris', source='genome model',
                          source_accession='GCA_fixture', bk_chromosome='chr1')
    r.alias(collision, 'collision_one')
    r.alias(collision, 'collision_two')
    records = [
        ('cynipoid|Acraspis_erinacei_1', 'locus_tree', {}),
        ('A_erinacei_pseudo_chr_12_000282', 'selection', {'link_note': 'GCA_035041075.1 model; no pinned current-model link'}),
        ('cynipoid|Andricus_quercusfoliatus_1', 'locus_tree', {}),
        ('AF_TRINITY_DN1015_c2_g1_i1.p1', 'selection', {}),
        (ena[1], 'catalytic_alignment', {}),
        (ena[1], 'origin_tree', {}),
        ('ambiguous_alias', 'selection', {}),
        ('collision_one', 'selection', {}),
        ('collision_two', 'selection', {}),
    ]
    for source, set_name, meta in (reversed(records) if reverse else records):
        r.collect('GH5_2', source, set_name, meta)
    r.assign_collected()
    assert r.source_owner['GH5_2', 'selection', 'A_erinacei_pseudo_chr_12_000282'] != current
    assert r.source_owner['GH5_2', 'selection', 'AF_TRINITY_DN1015_c2_g1_i1.p1'] == trinity
    assert r.source_owner['GH5_2', 'selection', 'ambiguous_alias'] not in (ambiguous_a, ambiguous_b)
    assert r.source_owner['GH5_2', 'selection', 'collision_one'] != collision
    assert r.source_owner['GH5_2', 'selection', 'collision_two'] != collision
    for source in ('collision_one', 'collision_two'):
        fallback = r.rows[r.source_owner['GH5_2', 'selection', source]]
        assert fallback['public_species_name'] == 'Pediaspis aceris'
        assert (fallback['source'], fallback['source_accession'], fallback['bk_chromosome']) == ('genome model', 'GCA_fixture', 'chr1')
    assert r.source_owner['GH5_2', 'catalytic_alignment', ena[1]] == ena
    assert r.source_owner['GH5_2', 'origin_tree', ena[1]] == ena
    return sorted((k, tuple(sorted(v.items()))) for k, v in r.rows.items()), sorted(r.source_owner.items())


def run():
    assert fixture() == fixture(reverse=True)
    old = dict(family='GH5_2', source_id='old', side='older', contig='chr10', strand='+', start='100', end='200')
    new = dict(old, source_id='new', side='current', start='150', end='250')
    assert unique_reciprocal_overlaps([old, new])['old'] == 'new'
    assert 'old' not in unique_reciprocal_overlaps([old, new, dict(new, source_id='other', start='170')])
    r = Resolver()
    r.collect('GH5_2', 'duplicate', 'origin_tree')
    r.collect('GH5_2', 'duplicate', 'origin_tree')
    try:
        r.assign_collected()
    except ValueError:
        pass
    else:
        raise AssertionError('duplicate source accepted')
    print('SEQUENCE_SET_FIXTURE_PASS')


if __name__ == '__main__':
    run()
