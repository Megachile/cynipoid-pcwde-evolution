"""Common CAZy rank-unit assignment for competitive dbCAN-sub rows.

The scan's own composition determines the finest named subfamily. Bare-family
leaders, ties, and groups without same-family composition stay at family level.
This rule applies to focal and native groups alike.
"""

import re


def assign(group: str, composition: str) -> tuple[str, float]:
    family = group.split("_e", 1)[0]
    entries = []
    for item in composition.split(";"):
        if ":" not in item:
            continue
        label, count = item.rsplit(":", 1)
        if label == family or re.fullmatch(re.escape(family) + r"_\d+", label):
            entries.append((label, int(count)))
    if not entries:
        return family, 0.0
    largest = max(count for _, count in entries)
    winners = [label for label, count in entries if count == largest]
    unit = winners[0] if len(winners) == 1 and winners[0] != family else family
    return unit, largest / sum(count for _, count in entries)
