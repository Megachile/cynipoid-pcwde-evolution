#!/usr/bin/env python3
"""Boundary tests for the shared, family-blind existing-order identity contract."""
import importlib.util
import tempfile
from pathlib import Path

SOURCE = Path(__file__).with_name("resolve_existing_members.py")
spec = importlib.util.spec_from_file_location("resolver", SOURCE)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def write(path, header, records):
    path.write_text("\t".join(header) + "\n" + "".join("\t".join(r) + "\n" for r in records))
    return path


def fails(text, *args):
    try:
        mod.resolve(*args)
    except ValueError as exc:
        assert text in str(exc), exc
    else:
        raise AssertionError(f"accepted invalid input; expected {text}")


def member(genome, family, pid, unit="", scaffold=""):
    return {"genome": genome, "family": family, "protein_id": pid, "unit_id": unit,
            "curated_scaffold": scaffold}


with tempfile.TemporaryDirectory() as tmp:
    d = Path(tmp)
    order = write(d / "order.tsv", ["genome", "chrom", "position", "protein_id"],
                  [["G1", "c1", "5", "F1"], ["G2", "c1", "5", "F2"], ["G1", "c1", "9", "M1"]])
    recon = write(d / "recon.tsv", ["genome", "scaffold", "position", "protein_id", "start_bp", "accepted_segment"],
                  [["G1", "s1", "0", "R1", "10", "1"], ["G2", "s1", "0", "R2", "10", "1"],
                   ["G1", "s2", "0", "S1", "10", "1"], ["G1", "s3", "0", "S1", "50", "1"],
                   ["G1", "s4", "0", "U1", "10", "1"], ["G1", "s5", "0", "M2", "10", "1"]])
    supp = write(d / "supp.tsv", ["protein_id", "family", "genome", "scaffold", "start_bp", "end_bp", "chromosome", "reason"],
                 [["R1", "x", "G1", "s1", "10", "20", "unplaced", "r"], ["R2", "x", "G2", "s1", "10", "20", "unplaced", "r"],
                  ["S1", "x", "G1", "s2", "10", "20", "unplaced", "r"], ["S1", "x", "G1", "s3", "50", "60", "unplaced", "r"],
                  ["V1", "x", "G1", "s9", "10", "20", "unplaced", "r"], ["M2", "x", "G1", "s5", "10", "20", "unplaced", "r"]])
    empty_mo = write(d / "mo_empty.tsv", ["genome", "protein_id", "scaffold", "position"], [])
    protein = d / "p.faa"
    protein.write_text(">M1\nAAA\n>M2\nCCC\n")
    gff = d / "m.gff3"
    gff.write_text("c1\tX\tmRNA\t1\t9\t.\t+\t.\tID=M1\ns5\tX\tmRNA\t1\t9\t.\t+\t.\tID=M2\n")
    src = {"G1": (protein, mod.digest(protein), gff, mod.digest(gff))}

    # Identical evidence in different families, both routes; family labels exchanged.
    cohort = [member("G1", "GH5_2", "F1"), member("G2", "B-lectin", "F2"),
              member("G1", "B-lectin", "R1"), member("G2", "PL1_8", "R2"),
              member("G1", "B-lectin", "S1"), member("G1", "PL4_4", "U1"), member("G1", "B-lectin", "V1"),
              member("G1", "PL4_4", "NOPOS")]
    swap = {"GH5_2": "B-lectin", "B-lectin": "GH5_2", "PL1_8": "B-lectin", "PL4_4": "B-lectin"}
    base = mod.resolve(cohort, order, empty_mo, {}, src, recon, supp)
    exchanged = mod.resolve([dict(m, family=swap[m["family"]]) for m in cohort], order, empty_mo, {}, src, recon, supp)
    view = lambda rs: [(r["accepted_id"], r["order_id"], r["order_route"], r["order_segments"], r["issue"]) for r in rs]
    assert view(base) == view(exchanged), "family label changed eligibility"
    got = {r["accepted_id"]: r for r in base}
    assert got["F1"]["order_route"] == got["F2"]["order_route"] == "final_order"
    assert got["R1"]["order_route"] == got["R2"]["order_route"] == "reconstructed_order"
    assert got["S1"]["order_segments"] == 2, "multi-segment identity collapsed"
    assert not got["U1"]["order_id"], "reconstructed row without supplemental record admitted"
    assert not got["V1"]["order_id"], "raw supplemental interval admitted as position"
    assert got["NOPOS"]["issue"] == "no_existing_order_position"

    # Member-layer order is a cross-check of the same positions.
    good = write(d / "mo_good.tsv", ["genome", "protein_id", "scaffold", "position"], [["G1", "F1", "c1", "5"], ["G1", "R1", "s1", "0"]])
    mod.resolve(cohort, order, good, {}, src, recon, supp)
    bad = write(d / "mo_bad.tsv", ["genome", "protein_id", "scaffold", "position"], [["G1", "F1", "c1", "6"]])
    fails("disagrees", cohort, order, bad, {}, src, recon, supp)

    # Unique-sequence alias through either route, with frame validation.
    alias = [member("G1", "B-lectin", "A1", "U", "c1"), member("G1", "GH5_2", "A2", "W", "s5")]
    res = mod.resolve(alias, order, empty_mo, {"U": "AAA", "W": "CCC"}, src, recon, supp)
    assert [(r["order_id"], r["order_route"]) for r in res] == [("M1", "final_order"), ("M2", "reconstructed_order")]
    fails("frame mismatch", [member("G1", "B-lectin", "A1", "U", "c9")], order, empty_mo, {"U": "AAA"}, src, recon, supp)
    res = mod.resolve([member("G1", "PL1_8", "A3", "U")], order, empty_mo, {"U": "GGG"}, src, recon, supp)
    assert res[0]["issue"] == "no_exact_aa_existing_order_model"
    protein.write_text(">M1\nAAA\n>M2\nAAA\n")
    src2 = {"G1": (protein, mod.digest(protein), gff, mod.digest(gff))}
    fails("ambiguous exact-AA", [member("G1", "GH5_2", "A1", "U")], order, empty_mo, {"U": "AAA"}, src2, recon, supp)

    # Collisions fail regardless of family.
    protein.write_text(">M1\nAAA\n>M2\nCCC\n")
    fails("cross-family", [member("G1", "GH5_2", "F1"), member("G2", "B-lectin", "F1")], order, empty_mo, {}, src, recon, supp)
    fails("one ordered model", [member("G1", "GH5_2", "M1"), member("G1", "B-lectin", "A1", "U", "c1")],
          order, empty_mo, {"U": "AAA"}, src, recon, supp)
    both = write(d / "both.tsv", ["genome", "chrom", "position", "protein_id"], [["G1", "c1", "5", "R1"]])
    fails("both order routes", [member("G1", "GH5_2", "R1")], both, empty_mo, {}, src, recon, supp)
    # Final-order identities on a reconstructed scaffold take the backend's
    # effective (reconstructed) position; a superseded row is not a position.
    eff = write(d / "eff_order.tsv", ["genome", "chrom", "position", "protein_id"],
                [["G2", "s1", "99", "Q2"], ["G2", "s1", "98", "Q3"]])
    eff_recon = write(d / "eff_recon.tsv", ["genome", "scaffold", "position", "protein_id", "start_bp", "accepted_segment"],
                      [["G2", "s1", "3", "Q2", "10", "0"]])
    res = mod.resolve([member("G2", "GH5_2", "Q2"), member("G2", "B-lectin", "Q3")], eff,
                      write(d / "mo_eff.tsv", ["genome", "protein_id", "scaffold", "position"], [["G2", "Q2", "s1", "3"]]),
                      {}, src, eff_recon, write(d / "nosupp.tsv", ["protein_id", "family", "genome", "scaffold", "start_bp", "end_bp", "chromosome", "reason"], []))
    assert [(r["order_id"], r["order_segments"]) for r in res] == [("Q2", 1), ("", 0)]
print("RESOLVE_EXISTING_MEMBERS_PASS")
