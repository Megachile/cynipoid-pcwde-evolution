#!/usr/bin/env python3
"""Freeze the intron producer's small, pinned outputs after genome-wide work.

This is an upstream producer, deliberately outside the Makefile. It consumes the
same source roster and corrected-model build used for the reproduction check.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
from collections import defaultdict
from pathlib import Path

from pyfaidx import Fasta
from Bio.Seq import Seq


def rows(path):
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def write(path, data, columns):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, columns, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(data)


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--data-root', type=Path, required=True)
    p.add_argument('--build', type=Path, required=True)
    p.add_argument('--sources', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=False)
    src = a.build / 'introns'
    frozen = a.out
    copied = []
    for sub, names in {
        'member_cds': ('cohort_accounting.tsv', 'member_cds_segments.tsv',
                       'blectin_accounting.tsv', 'blectin_cds_segments.tsv'),
        'host_background': ('QC.tsv', 'QC.json', 'whole_genome_comparator.json'),
    }.items():
        (frozen / sub).mkdir()
        for name in names:
            source = src / sub / name
            shutil.copyfile(source, frozen / sub / name)
            copied.append(frozen / sub / name)
    host = src / 'host_background'
    for row in rows(a.sources):
        name = row['source_file'] + '.tsv'
        shutil.copyfile(host / name, frozen / 'host_background' / name)
        copied.append(frozen / 'host_background' / name)
    details = frozen / 'host_background' / 'intron_details'
    details.mkdir()
    header = None
    handles = {}
    try:
        with (host / 'background_intron_details.tsv').open(newline='') as input_handle:
            reader = csv.reader(input_handle, delimiter='\t')
            header = next(reader)
            assert header[0] == 'source_file'
            for row in reader:
                name = row[0]
                if name not in handles:
                    output = details / (name + '.tsv')
                    h = output.open('w', newline='')
                    handles[name] = h
                    writer = csv.writer(h, delimiter='\t', lineterminator='\n')
                    writer.writerow(header)
                csv.writer(handles[name], delimiter='\t', lineterminator='\n').writerow(row)
    finally:
        for handle in handles.values():
            handle.close()
    assert set(handles) == {r['source_file'] for r in rows(a.sources)}
    copied.extend(details.glob('*.tsv'))
    source_map = {r['source_file']: r for r in rows(a.sources)}
    account = rows(frozen / 'member_cds' / 'cohort_accounting.tsv') + rows(frozen / 'member_cds' / 'blectin_accounting.tsv')
    segments = rows(frozen / 'member_cds' / 'member_cds_segments.tsv') + rows(frozen / 'member_cds' / 'blectin_cds_segments.tsv')
    by_protein = defaultdict(list)
    for row in segments:
        by_protein[row['protein_id']].append(row)
    assert len(account) == 657 and len(by_protein) == 631
    fasta = frozen / 'focal_cds.fna'
    junctions = []
    with fasta.open('w') as output:
        for source_name in sorted(source_map):
            subset = [r for r in account if r['source_file'] == source_name and r['status'] == 'single_scaffold']
            if not subset:
                continue
            source = source_map[source_name]
            assembly = a.data_root / source['assembly_snapshot_path']
            assert sha256(assembly) == source['assembly_sha256']
            genome = Fasta(str(assembly), indexname=str(src / 'indexes' / (source_name + '.fai')),
                           rebuild=False, as_raw=True, sequence_always_upper=True)
            for item in subset:
                pid = item['protein_id']
                parts = sorted(by_protein[pid], key=lambda r: int(r['transcript_order']))
                assert len(parts) == int(item['cds_segments'])
                strand = item['strand']
                seqs = [str(genome[p['scaffold']][int(p['cds_start']) - 1:int(p['cds_end'])]).upper() for p in parts]
                if strand == '-':
                    seqs = [str(Seq(s).reverse_complement()) for s in seqs]
                cds = ''.join(seqs)
                output.write(f'>{pid}\n{cds}\n')
                for number, (left, right) in enumerate(zip(parts, parts[1:]), 1):
                    start, end = ((int(left['cds_end']) + 1, int(right['cds_start']) - 1)
                                  if strand == '+' else
                                  (int(right['cds_end']) + 1, int(left['cds_start']) - 1))
                    if end < start:
                        continue
                    donor = str(genome[left['scaffold']][start - 1:min(start + 1, end)]).upper()
                    acceptor = str(genome[left['scaffold']][max(start - 1, end - 2):end]).upper()
                    if strand == '-':
                        donor, acceptor = str(Seq(acceptor).reverse_complement()), str(Seq(donor).reverse_complement())
                    junctions.append(dict(protein_id=pid, intron_number=number,
                                          donor_dinuc=donor, acceptor_dinuc=acceptor))
    assert len(junctions) > 0
    write(frozen / 'focal_junctions.tsv', junctions,
          ['protein_id', 'intron_number', 'donor_dinuc', 'acceptor_dinuc'])
    copied.extend((fasta, frozen / 'focal_junctions.tsv'))
    assert all(path.stat().st_size <= 200_000_000 for path in copied)
    manifest = [dict(path=str(path.relative_to(frozen)), bytes=path.stat().st_size, sha256=sha256(path))
                for path in sorted(copied)]
    write(frozen / 'FROZEN_MANIFEST.tsv', manifest, ['path', 'bytes', 'sha256'])
    (frozen / 'PRODUCER.json').write_text(json.dumps(dict(source_manifest_sha256=sha256(a.sources),
        sources=38, focal_copies=631, cohort=603, excluded=26, blectin=28,
        producer='upstream/introns/freeze_layer1.py'), indent=2) + '\n')
    print(json.dumps(dict(files=len(manifest), bytes=sum(r['bytes'] for r in manifest), junctions=len(junctions))))


if __name__ == '__main__':
    main()
