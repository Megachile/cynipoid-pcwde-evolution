#!/usr/bin/env python3
"""Extract the text before Google Docs' embedded image definitions from a pinned sync."""
import argparse, hashlib, re
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('source',type=Path);p.add_argument('out',type=Path);a=p.parse_args()
raw=a.source.read_bytes(); expected='18b32793af728c5342d7d0b2cb492edfb13307f12bfefca531b9452deb1fb198'
if hashlib.sha256(raw).hexdigest()!=expected: raise SystemExit('revision snapshot SHA256 mismatch')
text=raw.decode('utf-8').split('\n[image1]:',1)[0]
text=re.sub(r'^!?\[\]\[image\d+\].*\n','',text,flags=re.M)
if 'data:image' in text or re.search(r'^\[image\d+\]:',text,re.M): raise SystemExit('embedded image remained')
a.out.write_text(text)
