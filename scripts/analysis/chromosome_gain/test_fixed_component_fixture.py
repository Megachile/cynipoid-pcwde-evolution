#!/usr/bin/env python3
"""Bounded fixed-group selection and per-tip identity fixture."""
from collections import defaultdict

from ecami_score import fixed_component_groups

components = {
    ('a', 'g1'): {'a1', 'a2'}, ('b', 'g1'): {'b1'},
    ('a', 'g2'): {'a3'}, ('b', 'g2'): {'b2', 'b3'},
    ('a', 'p1'): {'p1'}, ('b', 'p1'): {'p2'},
    ('a', 'q1'): {'q1'}, ('b', 'q1'): {'q2'},
    ('a', 'native'): {'n1'},
}
mapping = {'g1': 'GH5_2', 'g2': 'GH5_2', 'p1': 'PL1_8',
           'q1': 'PL4_4', 'native': 'GH1'}
before = {key: set(genes) for key, genes in components.items()}
sizes, fixed = fixed_component_groups(components, mapping)
assert fixed == {'GH5_2': 'g1', 'PL1_8': 'p1', 'PL4_4': 'q1'}
assert len(sizes['GH5_2']['g1']) == len(sizes['GH5_2']['g2']) == 3
assert all(components.get((species, fixed['GH5_2'])) == expected
           for species, expected in [('a', {'a1', 'a2'}), ('b', {'b1'})])
assert components == before and components[('a', 'native')] == {'n1'}
print('FIXED_COMPONENT_FIXTURE_PASS fixed identity, lexical tie, primary/native inputs unchanged')
