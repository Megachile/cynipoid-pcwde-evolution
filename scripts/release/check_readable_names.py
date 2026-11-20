#!/usr/bin/env python3
"""Reject content-addressed deposit basenames and committed deposit pointers."""
import argparse
from pathlib import Path
import re
import subprocess
HASH_NAME = re.compile(r'^[0-9a-fA-F]{16,}(?:\.|$)')
POINTER = re.compile(r'repository_inputs/[^\s\t\"\'<>),;]+')
def check(repo, data_root):
    failures=[]
    for p in data_root.rglob('*'):
        if p.is_file() and HASH_NAME.match(p.name): failures.append(str(p))
    files=subprocess.check_output(['git','-C',str(repo),'ls-files'],text=True).splitlines()
    for rel in files:
        p=repo/rel
        try: text=p.read_text()
        except UnicodeError: continue
        for m in POINTER.finditer(text):
            if HASH_NAME.match(Path(m.group()).name): failures.append(rel+': '+m.group())
    if failures: raise ValueError('HASH_NAME_FAIL '+repr(failures[:10]))
    print('READABLE_NAMES_PASS')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--repo',type=Path,required=True);p.add_argument('--data-root',type=Path,required=True);a=p.parse_args();check(a.repo,a.data_root)
