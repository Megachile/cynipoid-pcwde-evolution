#!/usr/bin/env python3
"""Rebuild family-wide catalytic-position conservation from final alignments."""

from __future__ import annotations

import csv
import hashlib
import argparse
import shutil
import subprocess
from collections import Counter
from pathlib import Path

from Bio import SeqIO


MAFFT_THREADS_PRIMARY = 1
MAFFT_THREADS_FIXED_ALIGNMENT = 8

SITE_CONFIG = {
    "GH5_2": {
        "pdb": "1A3H.pdb",
        "chain": "A",
        "sites": [
            ("acid_base", 139, "E", "catalytic"),
            ("nucleophile", 228, "E", "catalytic"),
        ],
    },
    "PL1_8": {
        "pdb": "1QCX.pdb",
        "chain": "A",
        "sites": [
            ("supporting_arginine", 176, "R", "supporting_active_site"),
            ("catalytic_arginine", 236, "R", "catalytic"),
            ("supporting_lysine", 239, "K", "supporting_active_site"),
        ],
    },
    "PL4_4": {
        "pdb": "1NKG.pdb",
        "chain": "A",
        "sites": [
            ("proton_abstractor", 150, "K", "catalytic"),
            ("proton_donor", 210, "H", "catalytic"),
        ],
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Map reference-PDB catalytic residues onto accepted origin and locus "
            "alignments and summarize conservation in cynipoid domain sequences."
        )
    )
    parser.add_argument(
        "--paper-input-root",
        required=True,
        type=Path,
        help="Input root containing trees/origin and trees/locus.",
    )
    parser.add_argument(
        "--reference-dir",
        required=True,
        type=Path,
        help="Directory containing 1A3H.pdb, 1QCX.pdb, and 1NKG.pdb.",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="Directory for alignments and conservation tables.",
    )
    parser.add_argument(
        "--mafft",
        default=shutil.which("mafft"),
        help="MAFFT executable (default: mafft found on PATH).",
    )
    args = parser.parse_args()
    if not args.mafft:
        parser.error("MAFFT was not found on PATH; provide --mafft")
    return args

def is_cynipoid_sequence_id(sequence_id: str) -> bool:
    """Recognize legacy relabeled tips and current ENA genome tips."""
    source_prefix = sequence_id.split("|", 1)[0]
    return sequence_id.startswith("cynipoid|") or source_prefix.endswith("_ENA")


AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def pdb_chain(path: Path, chain: str) -> tuple[str, list[int]]:
    rows = []
    seen = set()
    with path.open() as handle:
        for line in handle:
            if not line.startswith("ATOM") or line[12:16].strip() != "CA":
                continue
            key = (line[21], int(line[22:26]), line[26])
            if line[21] != chain or key in seen:
                continue
            seen.add(key)
            rows.append((int(line[22:26]), line[26], AA3.get(line[17:20].strip(), "X")))
    rows.sort(key=lambda x: (x[0], x[1]))
    return "".join(x[2] for x in rows), [x[0] for x in rows]


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    paper = args.paper_input_root.resolve()
    reference_dir = args.reference_dir.resolve()
    out = args.output_dir.resolve()
    mafft = Path(args.mafft).resolve()
    out.mkdir(parents=True, exist_ok=True)
    detail = []
    summary = []
    membership_lock = []

    for family, site_cfg in SITE_CONFIG.items():
        cfg = dict(site_cfg)
        cfg["pdb"] = reference_dir / site_cfg["pdb"]
        if not cfg["pdb"].is_file():
            raise FileNotFoundError(f"Missing reference structure: {cfg['pdb']}")
        sequence, numbers = pdb_chain(cfg["pdb"], cfg["chain"])
        reference_faa = out / f"{family}_reference.faa"
        reference_faa.write_text(f">REFERENCE_{family}\n{sequence}\n")
        alignment_inputs = {
            "origin": paper / "trees" / "origin" / f"{family}_aligned.faa",
            "locus": paper / "trees" / "locus" / f"{family}_locus_aligned.faa",
        }
        for alignment_set, source_alignment in alignment_inputs.items():
            if not source_alignment.is_file():
                raise FileNotFoundError(f"Missing accepted alignment: {source_alignment}")
            augmented = out / f"{family}_{alignment_set}_with_reference_aligned.faa"
            command = [str(mafft), "--addfull", str(reference_faa), "--thread", str(MAFFT_THREADS_FIXED_ALIGNMENT), str(source_alignment)]
            with augmented.open("w") as handle:
                subprocess.run(command, check=True, stdout=handle)

            analyze_alignment(
                family, cfg, numbers, alignment_set, augmented, detail, summary
            )

        de_novo_input = out / f"{family}_de_novo_input.faa"
        combined_path = paper / "trees" / "origin" / f"{family}_combined.faa"
        if not combined_path.is_file():
            raise FileNotFoundError(f"Missing accepted domain set: {combined_path}")
        combined = combined_path.read_text()
        input_ids = [record.id for record in SeqIO.parse(combined_path, "fasta")]
        membership_lock.append({
            "family": family,
            "membership_source": str(combined_path),
            "sha256": hashlib.sha256(combined_path.read_bytes()).hexdigest(),
            "n_total_domain_sequences": len(input_ids),
            "n_cynipoid_domain_sequences": sum(is_cynipoid_sequence_id(x) for x in input_ids),
            "n_outgroup_domain_sequences": sum(not is_cynipoid_sequence_id(x) for x in input_ids),
            "membership_status": "locked_clean_20260710",
        })
        de_novo_input.write_text(combined + f">REFERENCE_{family}\n{sequence}\n")
        de_novo_alignment = out / f"{family}_de_novo_with_reference_aligned.faa"
        with de_novo_alignment.open("w") as handle:
            subprocess.run(
                [str(mafft), "--auto", "--thread", str(MAFFT_THREADS_PRIMARY), str(de_novo_input)],
                check=True,
                stdout=handle,
            )
        analyze_alignment(
            family, cfg, numbers, "de_novo", de_novo_alignment, detail, summary
        )

    write_tsv(out / "catalytic_site_conservation_by_sequence.tsv", detail, list(detail[0]))
    write_tsv(out / "catalytic_site_conservation_summary.tsv", summary, list(summary[0]))
    write_tsv(out / "membership_lock.tsv", membership_lock, list(membership_lock[0]))


def analyze_alignment(family, cfg, numbers, alignment_set, augmented, detail, summary):

    aligned = {record.id: str(record.seq) for record in SeqIO.parse(augmented, "fasta")}
    reference_id = f"REFERENCE_{family}"
    reference_aligned = aligned[reference_id]
    columns = {}
    residue_index = -1
    for column, aa in enumerate(reference_aligned):
        if aa == "-":
            continue
        residue_index += 1
        columns[numbers[residue_index]] = column

    cynipoids = sorted(name for name in aligned if is_cynipoid_sequence_id(name))
    for label, ref_resnum, expected, role in cfg["sites"]:
        if ref_resnum not in columns:
            raise RuntimeError(f"{family} reference residue {ref_resnum} was not aligned")
        column = columns[ref_resnum]
        observed = Counter()
        for sequence_id in cynipoids:
            aa = aligned[sequence_id][column].upper()
            status = "observed"
            if aa in "-.?X":
                status = "missing_or_gap"
            else:
                observed[aa] += 1
            detail.append({
                        "alignment_set": alignment_set,
                        "family": family,
                        "sequence_id": sequence_id,
                        "counting_unit": "cynipoid_domain_sequence",
                        "reference_pdb": cfg["pdb"].name,
                        "reference_chain": cfg["chain"],
                        "site_label": label,
                        "site_role": role,
                        "reference_residue": f"{expected}{ref_resnum}",
                        "alignment_column_1based": column + 1,
                        "observed_aa": aa,
                        "site_status": status,
                        "conserved": aa == expected,
            })
        n_total = len(cynipoids)
        n_observed = sum(observed.values())
        n_conserved = observed[expected]
        summary.append({
                    "alignment_set": alignment_set,
                    "family": family,
                    "counting_unit": "cynipoid_domain_sequence",
                    "reference_pdb": cfg["pdb"].name,
                    "site_label": label,
                    "site_role": role,
                    "reference_residue": f"{expected}{ref_resnum}",
                    "n_total_sequences": n_total,
                    "n_site_observed": n_observed,
                    "n_missing_or_gap": n_total - n_observed,
                    "n_expected_residue": n_conserved,
                    "pct_expected_among_observed": f"{100 * n_conserved / n_observed:.3f}" if n_observed else "NA",
                    "pct_expected_among_all": f"{100 * n_conserved / n_total:.3f}" if n_total else "NA",
                    "observed_residue_counts": ";".join(f"{aa}:{n}" for aa, n in sorted(observed.items())),
        })


if __name__ == "__main__":
    main()
