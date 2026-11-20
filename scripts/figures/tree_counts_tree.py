"""Graft the accepted Figure 2 cynipoid display tree onto Figure 1."""

from copy import deepcopy
from Bio import Phylo


def clade_signature(tree, allowed=None):
    """Rooted nontrivial tip sets, optionally induced on a subset of tips."""
    allowed = set(allowed or (tip.name for tip in tree.get_terminals()))
    result = set()
    for clade in tree.find_clades():
        members = frozenset(tip.name for tip in clade.get_terminals()
                            if tip.name in allowed)
        if 1 < len(members) < len(allowed):
            result.add(members)
    return result


def verify_graft(fig1_tree, fig2_tree, combined, fig2_prune=()):
    """Check the entire Figure 2 ingroup and the Figure 1 outgroup backbone."""
    inner_tips = ({tip.name for tip in fig2_tree.get_terminals()}
                  - {'Telenomus_remus'} - set(fig2_prune))
    outer_tips = {tip.name for tip in fig1_tree.get_terminals()} - {
        tip.name for tip in fig1_tree.common_ancestor(
            'Diplolepis_rosae', 'Belonocnema_kinseyi').get_terminals()}
    expected = inner_tips | outer_tips
    observed = [tip.name for tip in combined.get_terminals()]
    if len(observed) != len(expected) or set(observed) != expected:
        raise AssertionError('Combined tip roster differs from graft inputs')
    inner = deepcopy(fig2_tree)
    inner.root_with_outgroup('Telenomus_remus')
    for name in fig2_prune:
        if name in {tip.name for tip in inner.get_terminals()}:
            inner.prune(name)
    inner.prune('Telenomus_remus')
    if clade_signature(combined, inner_tips) != clade_signature(inner):
        raise AssertionError('Cynipoid subtree differs from Figure 2')
    # Collapse the original and replacement cynipoid clades to one meta-tip.
    # This retains the Figure 1 topology of the non-cynipoid backbone.
    def collapsed_signature(tree, cynipoid_tips):
        result = set()
        universe = outer_tips | {'Cynipoidea'}
        for clade in tree.find_clades():
            members = {tip.name for tip in clade.get_terminals()}
            reduced = frozenset((members & outer_tips) |
                                ({'Cynipoidea'} if cynipoid_tips <= members else set()))
            if 1 < len(reduced) < len(universe):
                result.add(reduced)
        return result
    old_inner = {tip.name for tip in fig1_tree.common_ancestor(
        'Diplolepis_rosae', 'Belonocnema_kinseyi').get_terminals()}
    if collapsed_signature(combined, inner_tips) != collapsed_signature(fig1_tree, old_inner):
        raise AssertionError('Outgroup backbone differs from Figure 1')


def graft(fig1_tree, fig2_tree, fig2_prune=()):
    """Return a tree retaining Figure 1's outgroups and Figure 2's cynipoids."""
    outer = deepcopy(fig1_tree)
    inner = deepcopy(fig2_tree)
    inner.root_with_outgroup('Telenomus_remus')
    for name in fig2_prune:
        if name in {tip.name for tip in inner.get_terminals()}:
            inner.prune(name)
    inner.prune('Telenomus_remus')
    original = outer.common_ancestor(
        'Diplolepis_rosae', 'Belonocnema_kinseyi')
    expected = {tip.name for tip in original.get_terminals()}
    actual = {tip.name for tip in inner.get_terminals()}
    if not expected <= actual or 'Telenomus_remus' in actual:
        raise ValueError('Figure 1 cynipoid tip mapping cannot be grafted')
    for parent in outer.find_clades():
        for index, child in enumerate(parent.clades):
            if child is original:
                inner.root.branch_length = child.branch_length
                parent.clades[index] = inner.root
                names = [tip.name for tip in outer.get_terminals()]
                if len(names) != len(set(names)):
                    raise ValueError('Graft has duplicate tips')
                return outer
    raise ValueError('Figure 1 cynipoid node has no parent')
