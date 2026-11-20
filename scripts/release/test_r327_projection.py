#!/usr/bin/env python3
"""Declared GFF projection boundary regression."""
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import public_inputs
from test_public_projection import fixture, put, PRIVATE_FIXTURE_PATH

SOURCE = Path(__file__).with_name("build_data_deposit.py")

def run(root, inputs):
    env = dict(os.environ, PCWDE_REPO_ROOT=str(root), PYTHONDONTWRITEBYTECODE="1")
    return subprocess.run([sys.executable, "-B", str(SOURCE), "--project-inputs-only", str(inputs)],
                          env=env, text=True, capture_output=True)

with tempfile.TemporaryDirectory() as temp:
    root = Path(temp)
    declared = root / "declared"
    inputs = fixture(declared, True)
    result = run(declared, inputs)
    assert result.returncode == 0, result.stderr
    assert b"# private source path omitted\n" in (inputs / "analysis/introns/sources/A/gff.gff").read_bytes()
    assert public_inputs.public_problem_set(inputs) == set()

    undeclared = root / "undeclared"
    inputs = fixture(undeclared, False)
    put(inputs / "analysis/introns/sources/B/gff.gff",
        f"##gff-version 3\n# source {PRIVATE_FIXTURE_PATH}/model\nchr\tHelixer\tgene\t1\t2\t.\t+\t.\tID=x\n")
    result = run(undeclared, inputs)
    assert result.returncode != 0 and "private_absolute_path" in result.stderr, result.stderr

    partial = root / "partial"
    inputs = fixture(partial, True)
    public_inputs.project_gff_header(inputs / "analysis/introns/sources/A/gff.gff")
    result = run(partial, inputs)
    assert result.returncode != 0 and "partially projected intron GFF" in result.stderr, result.stderr

    namespace = {"__file__": str(SOURCE), "__name__": "deposit_r327_test"}
    source = SOURCE.read_text().split("if len(sys.argv) == 2 and sys.argv[1] == \"--validate-inputs-only\":")[0]
    exec(compile(source, str(SOURCE), "exec"), namespace)
    namespace["ROOT"] = declared
    assert namespace["declared_gff_paths"](declared / "inputs") == ["analysis/introns/sources/A/gff.gff"]
    assert namespace["declared_gff_paths"](undeclared / "inputs") == []
    assert ("private_absolute_path", "analysis/introns/sources/B/gff.gff") not in namespace["expected_preprojection_problems"](undeclared / "inputs")

    already = root / "already"
    already_inputs = fixture(already, True)
    public_inputs.project_repository_inputs(already_inputs, already)
    put(already_inputs / "analysis/introns/sources/A/gff.gff",
        "##gff-version 3\n# public source\nchr\tHelixer\tgene\t1\t2\t.\t+\t.\tID=x\n")
    namespace["ROOT"] = already
    namespace["rebuild_input_inventory"](already_inputs)
    already_rows = namespace["project_repository_inputs"](already_inputs)
    assert [(r["target"], r["rows"]) for r in already_rows if r["operation"] == "neutralize_private_gff_header"] == [("analysis/introns/sources/A/gff.gff", "0")]

    stage = Path(os.environ["PCWDE_DATA_ROOT"]).resolve(strict=True)
    prior = stage / "PUBLIC_PROJECTION.tsv"
    assert prior.is_file(), f"missing sealed projection receipt: {prior}"
    assert namespace["read_tsv"](prior)[0] == ["operation", "target", "rows", "detail"]
    assert namespace["hashlib"].sha256(prior.read_bytes()).hexdigest() == namespace["ACCEPTED_PRIOR_PROJECTION_SHA256"]
    for relative in namespace["S4_SEQUENCE_DIGEST_TSVS"]:
        header, _ = namespace["read_tsv"](stage / "repository_inputs" / relative)
        assert "sequence_sha256" in header and "sequence_sha256" not in namespace["hash_fields"](header, relative)

print("R327_PROJECTION_PASS declared undeclared partial S4_digest_prior_receipt")
