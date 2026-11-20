#!/usr/bin/env python3
"""Small baseline/candidate fixtures for the public deposit projection."""
import csv
import importlib.util
import pathlib
import os
import subprocess
import sys
import tempfile


import public_inputs as projection
SOURCE = pathlib.Path(__file__).with_name("build_data_deposit.py")


def put(path, content):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content)


PRIVATE_FIXTURE_PATH = "/" + "carc/private"


def rejects(action, expected):
    try:
        action()
    except SystemExit as exc:
        assert expected in str(exc), (expected, str(exc))
    else:
        raise AssertionError(f"expected rejection: {expected}")


def fixture(root, with_intron):
    inputs = root / "inputs"
    put(root / "config/intron_sources.tsv", "source_file\nA\n")
    for relative in projection.ROW_HASH_TSVS:
        put(inputs / relative, "id\tsha256\nx\tabc\n")
    put(inputs / "tables/authority/synteny_current/producer_build_receipt.tsv",
        "id\tsha256\nx\tDCAN\n")
    put(inputs / projection.CDS_MAPPING, f"id\tcds_path\nx\t{PRIVATE_FIXTURE_PATH}/cds\n")
    figures = inputs / "figures"
    existing = [p.relative_to(figures).as_posix() for p in figures.rglob("*") if p.is_file()
                and p != inputs / projection.FIGURE_MANIFEST]
    for i in range(8):
        put(figures / f"fixture_{i}.txt", "public\n")
        existing.append(f"fixture_{i}.txt")
    rows = [["fixture", rel, "no-DCAN" if i < 8 else "public", "0", "abc"]
            for i, rel in enumerate(existing)]
    with (inputs / projection.FIGURE_MANIFEST).open("w", newline="") as handle:
        out = csv.writer(handle, delimiter="\t", lineterminator="\n")
        out.writerow(["component", "relative_path", "role", "size_bytes", "sha256"])
        out.writerows(rows)
    if with_intron:
        put(inputs / "analysis/introns/sources/A/gff.gff",
            f"##gff-version 3\n# source {PRIVATE_FIXTURE_PATH}/model\nchr\tHelixer\tgene\t1\t2\t.\t+\t.\tID=x\n")
    return inputs


if __name__ == "__main__":
    with tempfile.TemporaryDirectory() as tmp:
        base_root = pathlib.Path(tmp) / "base"
        cand_root = pathlib.Path(tmp) / "candidate"
        base = fixture(base_root, False)
        cand = fixture(cand_root, True)
        put(cand / "taxonomy.tsv", "name\nDCAN\n")
        assert ("removed_taxon_token", "taxonomy.tsv") in projection.public_problem_set(cand)
        (cand / "taxonomy.tsv").unlink()
        put(cand / "mapping.gff", "chr\tNCBI\tgene\t1\t2\t.\t+\t.\tgene_synonym=DCAN\n")
        assert ("removed_taxon_token", "mapping.gff") not in projection.public_problem_set(cand)
        (cand / "mapping.gff").unlink()
        partial_root = pathlib.Path(tmp) / "partial"
        partial = fixture(partial_root, True)
        repo_root = partial_root
        projection.project_gff_header(partial / "analysis/introns/sources/A/gff.gff")
        rejects(lambda: projection.project_repository_inputs(partial, repo_root), "partially projected intron GFF sources")
        feature_root = pathlib.Path(tmp) / "feature"
        feature = fixture(feature_root, True)
        put(feature / "analysis/introns/sources/A/gff.gff",
            f"##gff-version 3\nchr\tHelixer\tgene\t1\t2\t.\t+\t.\tID={PRIVATE_FIXTURE_PATH}/model\n")
        repo_root = feature_root
        rejects(lambda: projection.project_repository_inputs(feature, repo_root), "private path in GFF feature line")
        unexpected_root = pathlib.Path(tmp) / "unexpected"
        unexpected = fixture(unexpected_root, True)
        put(unexpected / "unexpected.txt", f"{PRIVATE_FIXTURE_PATH}/path\n")
        repo_root = unexpected_root
        rejects(lambda: projection.project_repository_inputs(unexpected, repo_root), "extra=[('private_absolute_path', 'unexpected.txt')]")
        cli_root = pathlib.Path(tmp) / "cli"
        cli = fixture(cli_root, True)
        env = dict(os.environ, PCWDE_REPO_ROOT=str(cli_root))
        command = subprocess.run([sys.executable, str(SOURCE), "--project-inputs-only", str(cli)],
                                 env=env, capture_output=True, text=True)
        assert command.returncode == 0, command.stderr
        assert "PUBLIC_PROJECTION_PASS" in command.stdout
        assert projection.public_problem_set(cli) == set()
        repo_root = base_root
        base_spec = projection.project_repository_inputs(base, base_root)
        repo_root = cand_root
        cand_spec = projection.project_repository_inputs(cand, cand_root)
        assert len([r for r in cand_spec if r["operation"] == "neutralize_private_gff_header"]) == 1
        assert (cand / "analysis/introns/sources/A/gff.gff").read_text().splitlines()[-1] == \
               "chr\tHelixer\tgene\t1\t2\t.\t+\t.\tID=x"
        assert projection.public_problem_set(cand) == set()
        for path in base.rglob("*"):
            if path.is_file() and path.name != "INPUT_MANIFEST.tsv":
                assert path.read_bytes() == (cand / path.relative_to(base)).read_bytes()
        assert projection.validate_projected_inventory(cand) > projection.validate_projected_inventory(base)
        put(cand / "unexpected.txt", f"{PRIVATE_FIXTURE_PATH}/path\n")
        rejects(lambda: projection.project_repository_inputs(cand, cand_root), "partially projected intron GFF sources")
        (cand / "unexpected.txt").unlink()
        (cand / "analysis/introns/sources/A/gff.gff").unlink()
        rejects(lambda: projection.validate_projected_inventory(cand), "projected repository inventory closure mismatch")
    print("PUBLIC_PROJECTION_UNIT_PASS")
