#!/usr/bin/env python3
"""Canonical historical PCWDE apparent-fusion predicate.

The historical filter is conjunctive and uses aligned target-domain
coordinates.  A cynipoid protein is excluded only when it has exactly one
selected target-domain hit, is longer than twice its family median, and that
aligned domain covers less than half of the protein.  Multiple-domain tandems,
outgroups, and either strict boundary are retained.
"""

from __future__ import annotations

FUSION_LENGTH_MULTIPLIER = 2.0
FUSION_ALIGNED_COVERAGE_MAX = 0.50


def aligned_target_coverage(protein_length: int, aligned_target_length: int) -> float:
    """Return aligned target-domain coverage of the full protein."""
    if protein_length <= 0:
        raise ValueError("protein_length must be positive")
    if aligned_target_length < 0 or aligned_target_length > protein_length:
        raise ValueError("aligned_target_length must be within the protein")
    return aligned_target_length / protein_length


def is_apparent_gene_prediction_fusion(
    *,
    selected_target_domains: int,
    protein_length: int,
    family_median_length: float,
    aligned_target_length: int,
    is_outgroup: bool = False,
) -> bool:
    """Apply the locked historical AND rule using aligned coordinates.

    The combination of exactly one selected target domain and <50% aligned
    target coverage is the executable evidence for a single partial target
    domain in an otherwise unrelated ORF.  No gappiness, expression, fragment,
    or manual-ID criterion is part of this decision.
    """
    if selected_target_domains < 0:
        raise ValueError("selected_target_domains cannot be negative")
    if family_median_length <= 0:
        raise ValueError("family_median_length must be positive")
    coverage = aligned_target_coverage(protein_length, aligned_target_length)
    return (
        not is_outgroup
        and selected_target_domains == 1
        and protein_length > FUSION_LENGTH_MULTIPLIER * family_median_length
        and coverage < FUSION_ALIGNED_COVERAGE_MAX
    )
