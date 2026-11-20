"""Fail closed if excluded public payload or legacy identifier fields reappear."""
import csv,hashlib
from pathlib import Path

def validate(inputs, repo, check_pins=True):
    inputs,repo=Path(inputs),Path(repo)
    retired=[s for s in (repo/'config/public_retired_inputs.txt').read_text().splitlines() if s and not s.startswith('#')]
    forbidden=[rel for rel in retired if (inputs/rel).exists()]
    forbidden += [str(p.relative_to(inputs)) for p in (inputs/'figures/origin_trees').glob('MC_*_domains.faa')]
    if forbidden:raise ValueError(f'retired public inputs reappeared: {forbidden}')
    for p in inputs.rglob('*.tsv'):
        with p.open(newline='') as f:
            reader=csv.reader(f,delimiter='\t');header=next(reader,[])
            if {'mc_v1','cynogv4_id'} & set(header):raise ValueError(f'legacy identifier field: {p}')
            for row in reader:
                if any('BIP_adult_OVX' in c or 'BIP_adult_VGX' in c for c in row):raise ValueError(f'excluded expression library: {p}')
                if any(i<len(row) and row[i] in {'EXPSAMPLE-019','EXPSAMPLE-020'} for i,k in enumerate(header) if k=='sample_id'):
                    raise ValueError(f'excluded expression sample ID: {p}')
    # Bind transformed scientific inputs to the committed projection receipt.
    # Top-level inventories are regenerated after the GFF comment projection.
    if not check_pins:return
    with (repo/'config/public_scope_pins.tsv').open(newline='') as f:
        pins=list(csv.DictReader(f,delimiter='\t'))
    for row in pins:
        rel=row['path'].removeprefix('repository_inputs/')
        if row['new_sha256']=='retired':
            if (inputs/rel).exists():raise ValueError(f'retired pinned input reappeared: {rel}')
        elif rel!='INPUT_MANIFEST.tsv':
            p=inputs/rel
            if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=row['new_sha256']:
                raise ValueError(f'public projection pin mismatch: {rel}')
