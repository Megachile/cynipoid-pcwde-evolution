#!/usr/bin/env python3
"""Build protein-guided codon alignments for locked current PCWDE members."""

from __future__ import annotations

import csv
import hashlib
import argparse
import os
import re
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq


FAMILIES = ("GH5_2", "PL1_8", "PL4_4")
MAX_GAP = 0.50
MIN_RECOVERY = 0.95
MAFFT_VERSION = "7.526"


def read_fasta(path: Path) -> dict[str, str]:
    return {r.id: str(r.seq).upper() for r in SeqIO.parse(path, "fasta")}


def write_fasta(path: Path, records: dict[str, str]) -> None:
    with path.open("w") as handle:
        for name, seq in records.items():
            handle.write(f">{name}\n{seq}\n")


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def clean_aa(seq: str) -> str:
    return seq.upper().rstrip("*")


def clean_cds(seq: str) -> str:
    seq = seq.upper().replace("-", "")
    if seq[-3:] in ("TAA", "TAG", "TGA"):
        seq = seq[:-3]
    return seq


def translate(seq: str) -> str:
    return str(Seq(clean_cds(seq)).translate()) if len(clean_cds(seq)) % 3 == 0 else ""


def analysis_id(protein_id: str) -> str:
    return re.sub(r"[^A-Za-z0-9_-]", "_", protein_id)


def trinity_group(protein_id: str) -> str:
    if "TRINITY_" not in protein_id:
        return ""
    value = re.sub(r"\.p\d+$", "", protein_id)
    return re.sub(r"_i\d+$", "", value)


def sha256_records(records: dict[str, str]) -> str:
    payload = "".join(f">{name}\n{records[name]}\n" for name in records)
    return hashlib.sha256(payload.encode()).hexdigest()


def require_mafft(path: Path) -> Path:
    executable = path.expanduser().resolve()
    if not executable.is_file() or not os.access(executable, os.X_OK):
        raise SystemExit(f"MAFFT executable is missing or not executable: {executable}")
    check = subprocess.run(
        [str(executable), "--version"], text=True, capture_output=True, check=False,
    )
    version_text = f"{check.stdout}\n{check.stderr}"
    if check.returncode != 0 or MAFFT_VERSION not in version_text:
        raise SystemExit(
            f"expected MAFFT {MAFFT_VERSION}, got exit {check.returncode}: "
            f"{version_text.strip()!r}"
        )
    return executable


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build protein-guided codon alignments for the accepted PCWDE membership."
    )
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument(
        "--protein-dir", required=True, type=Path,
        help="Directory containing GH5_2.faa, PL1_8.faa, and PL4_4.faa.",
    )
    parser.add_argument(
        "--cds-dir", required=True, type=Path,
        help="Directory containing <family>_locked_current_cds.fna.",
    )
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--mafft", required=True, type=Path)
    parser.add_argument(
        "--mafft-threads", type=int, choices=(8,), default=8,
        help="MAFFT threads (default: 8, matching the reported analysis).",
    )
    args = parser.parse_args()
    mafft = require_mafft(args.mafft)
    if not args.metadata.is_file():
        parser.error(f"metadata does not exist: {args.metadata}")
    for label, directory in (("protein", args.protein_dir), ("CDS", args.cds_dir)):
        if not directory.is_dir():
            parser.error(f"{label} directory does not exist: {directory}")
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=False)
    metadata = list(csv.DictReader(args.metadata.open(), delimiter="\t"))
    meta_by_family = {
        family: {r["protein_id"]: r for r in metadata if r["family"] == family}
        for family in FAMILIES
    }
    membership_rows: list[dict] = []
    qc_rows: list[dict] = []

    for family in FAMILIES:
        all_proteins = read_fasta(args.protein_dir / f"{family}.faa")
        current = {
            protein_id: clean_aa(all_proteins[protein_id])
            for protein_id in meta_by_family[family]
            if protein_id in all_proteins
        }
        cds_all = read_fasta(args.cds_dir / f"{family}_locked_current_cds.fna")
        valid: dict[str, str] = {}
        reason: dict[str, str] = {}
        for protein_id, aa in current.items():
            if protein_id not in cds_all:
                reason[protein_id] = "no_exact_cds"
                continue
            cds = clean_cds(cds_all[protein_id])
            if len(cds) % 3:
                reason[protein_id] = "cds_not_multiple_of_3"
            elif "*" in translate(cds):
                reason[protein_id] = "internal_stop"
            elif translate(cds) != aa:
                reason[protein_id] = "translation_mismatch"
            else:
                valid[protein_id] = cds
                reason[protein_id] = "included"

        recovery = len(valid) / len(current)
        if recovery < MIN_RECOVERY:
            raise SystemExit(
                f"{family}: exact CDS recovery {len(valid)}/{len(current)} ({recovery:.1%}) below {MIN_RECOVERY:.0%}"
            )

        id_map = {protein_id: analysis_id(protein_id) for protein_id in valid}
        collisions = [name for name, n in Counter(id_map.values()).items() if n > 1]
        if collisions:
            raise SystemExit(f"{family}: sanitized ID collision(s): {collisions[:5]}")

        aa_input = {id_map[name]: current[name] for name in valid}
        aa_path = out / f"{family}_current_aa.faa"
        aa_aligned_path = out / f"{family}_current_aa_aligned.faa"
        write_fasta(aa_path, aa_input)
        with aa_aligned_path.open("w") as handle:
            subprocess.run(
                [str(mafft), "--auto", "--thread", str(args.mafft_threads), str(aa_path)],
                stdout=handle, check=True,
            )
        aligned = read_fasta(aa_aligned_path)
        if set(aligned) != set(aa_input):
            raise SystemExit(f"{family}: MAFFT changed sequence IDs")

        codon: dict[str, str] = {}
        reverse_map = {value: key for key, value in id_map.items()}
        for aid, aa_aln in aligned.items():
            protein_id = reverse_map[aid]
            cds = valid[protein_id]
            chunks = []
            position = 0
            for aa in aa_aln:
                if aa in "-.":
                    chunks.append("---")
                else:
                    chunks.append(cds[position:position + 3])
                    position += 3
            if position != len(cds):
                raise SystemExit(f"{family}:{protein_id}: CDS threading consumed {position}/{len(cds)} nt")
            threaded = "".join(chunks)
            if translate(threaded.replace("---", "")) != current[protein_id]:
                raise SystemExit(f"{family}:{protein_id}: threaded translation mismatch")
            codon[aid] = threaded

        n_before = len(next(iter(codon.values()))) // 3
        keep = []
        for i in range(n_before):
            gap_count = sum(seq[3 * i:3 * i + 3] == "---" for seq in codon.values())
            if gap_count / len(codon) <= MAX_GAP:
                keep.append(i)
        filtered = {
            name: "".join(seq[3 * i:3 * i + 3] for i in keep)
            for name, seq in codon.items()
        }
        if any(len(seq) % 3 for seq in filtered.values()):
            raise SystemExit(f"{family}: non-codon filtered alignment")
        write_fasta(out / f"{family}_current_codon_aligned.fna", filtered)

        group_counts = Counter(trinity_group(name) for name in valid if trinity_group(name))
        source_counts = Counter(meta_by_family[family][name]["source"] for name in valid)
        duplicate_aa = sum(n > 1 for n in Counter(current[name] for name in valid).values())
        duplicate_cds = sum(n > 1 for n in Counter(valid.values()).values())
        qc_rows.append({
            "family": family,
            "locked_full_protein_targets": len(current),
            "exact_cds_included": len(valid),
            "excluded_no_exact_cds": len(current) - len(valid),
            "recovery_percent": f"{100 * recovery:.3f}",
            "genome_members": source_counts.get("genome", 0),
            "transcriptome_members": source_counts.get("transcriptome", 0),
            "trinity_groups_with_multiple_isoforms": sum(n > 1 for n in group_counts.values()),
            "duplicate_aa_sequence_groups": duplicate_aa,
            "duplicate_cds_sequence_groups": duplicate_cds,
            "protein_alignment_columns": len(next(iter(aligned.values()))),
            "codon_columns_before_filter": n_before,
            "codon_columns_after_filter": len(keep),
            "codon_columns_removed_gt50pct_gap": n_before - len(keep),
            "locked_protein_sha256": sha256_records(current),
            "included_cds_sha256": sha256_records(valid),
        })

        for protein_id in current:
            meta = meta_by_family[family][protein_id]
            group = trinity_group(protein_id)
            membership_rows.append({
                "family": family,
                "protein_id": protein_id,
                "analysis_id": id_map.get(protein_id, ""),
                "species": meta["species"],
                "source": meta["source"],
                "source_file": meta["source_file"],
                "trinity_gene_group": group,
                "n_isoforms_in_group": group_counts.get(group, 0) if group else 0,
                "aa_length": len(current[protein_id]),
                "cds_length_nt": len(valid.get(protein_id, "")),
                "selection_status": reason[protein_id],
            })

        print(
            f"{family}: {len(valid)}/{len(current)} exact CDS; "
            f"{len(keep)}/{n_before} codon columns retained"
        )

    write_tsv(
        out / "selection_membership_lock.tsv", membership_rows,
        ["family", "protein_id", "analysis_id", "species", "source", "source_file",
         "trinity_gene_group", "n_isoforms_in_group", "aa_length", "cds_length_nt", "selection_status"],
    )
    write_tsv(out / "alignment_qc.tsv", qc_rows, list(qc_rows[0]))


if __name__ == "__main__":
    main()
