#!/usr/bin/env python3
"""Check that candidate recipes preserve the ordering of integration-main commands."""
import argparse
from pathlib import Path


def recipes(path):
    result = {}
    target = None
    for line in Path(path).read_text().splitlines():
        if line and not line[0].isspace() and ':' in line and not line.startswith('.'):
            target = line.split(':', 1)[0]
            if target and all(c.isalnum() or c in '_-' for c in target):
                result[target] = []
            else:
                target = None
        elif target and line.startswith('\t'):
            result[target].append(line.strip())
        elif line and not line[0].isspace():
            target = None
    return result



def check(base, candidate):
    common = set(base) & set(candidate)
    for target in sorted(common):
        actual = candidate[target]
        position = 0
        for line in base[target]:
            try:
                position = actual.index(line, position) + 1
            except ValueError:
                raise AssertionError(f'{target}: missing or reordered main recipe line: {line}')
    return len(common)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('main_makefile')
    parser.add_argument('candidate_makefile')
    args = parser.parse_args()
    count = check(recipes(args.main_makefile), recipes(args.candidate_makefile))
    print(f'RECIPE_ORDER_PASS shared_targets={count}')
