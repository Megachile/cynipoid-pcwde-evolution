#!/usr/bin/env python3
"""Rebuild the revised deposit from the submitted deposit and pinned upstream outputs.

The revised `repository_inputs/` differs from the submitted Figshare deposit by
six staging steps, applied in this order:

  bp_cb_pruning    remove the CB-origin Bassettia pallida tips (prune_trees.py,
                   drop_rows.py, bp_cb_prune_plan.tsv) and their domain-preparation records
  copy_table       shared/copy_placements.tsv from the placement-conformance table
  computed_tables  drop saved Tables S2-S4; make computes them
  locus_trees      Figure 4 / S3 inputs and shared/locus_tips.tsv from the locus-tree rebuild
  reader_locus     reader-deposit complete-locus files from the same rebuild
  gh1_counts       Figure 1 GH1_e72 counts
  tblastn_presence shared/tblastn_presence.tsv: validated tBLASTn loci for genome/family pairs without a member
  member_layer     members/member_table.tsv + members/member_order.tsv: the frozen member layer make builds copies from
  derived_locus_tips members/locus_tip_crosswalk.tsv: tree-intrinsic tip crosswalk; make derives species/chromosome
  derived_copy_table  drop shared/copy_placements.tsv; make builds it from the member layer
  retired_figures  drop inputs of figures removed from the paper (synteny figure, 2026-10-29)
  stable_names     figure input folders take stable names (figures/figure3 -> figures/origin_trees, ...)
  corrected_synteny three frozen order/label inputs plus source-derived native
                   membership, block summary and recovery states for Table S5
  supplement_facts  source reports of supplement prose facts (E. primus sequencing,
                   assembly and annotation; CynOG catalog build; B-lectin outgroup
                   screen; raw-genome tblastn genome manifest), pinned in
                   config/upstream_locked_inputs.tsv rows with genome "supplement"

Upstream files are named relative to the analysis workspace and must match the
sha256 pinned in config/deposit_sources.tsv. Heavy analysis (placement,
alignment, tree inference) is not rerun here.

Usage:
  stage.py --base SUBMITTED_DATA_ROOT --workspace WORKSPACE --out NEW_DATA_ROOT [--record]
  stage.py --compare A_DATA_ROOT B_DATA_ROOT   # same form (raw or projected) required
  stage.py --pin DATA_ROOT     # point config/logical_inputs.tsv at this deposit

Environment: project Python 3 (standard library only).
"""

import sys as _source_sys
from pathlib import Path as _SourcePath
_source_repo = next(p for p in _SourcePath(__file__).resolve().parents if (p / 'config/logical_inputs.tsv').is_file())
_source_sys.path.insert(0, str(_source_repo / 'scripts/release'))
from source_aliases import source_path
import argparse
import csv
import hashlib
import re
import importlib.util
import os
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
SOURCES = REPO / "config/deposit_sources.tsv"
CHROMOSOME_GAIN_SOURCES = REPO / "config/chromosome_gain_sources.tsv"
CAZY_SOURCES = REPO / "config/cazy_stage_sources.tsv"
LOCKED_INPUTS = REPO / "config/upstream_locked_inputs.tsv"
RESCUE_SOURCE = "analysis/raw_genome_tblastn_presence_loop/outputs/fig5_validated_raw_locus_rescue.tsv"
RESCUE_SHA256 = "983a6f364d339d2066b05d4b962849575209b370e8f99278281fdb74fac8215f"
MANIFESTS = ("INPUT_MANIFEST.tsv", "figures/INPUT_MANIFEST.tsv")
EXCLUDED = {"figures/figure_s5/Figure_S5.pdf"}  # static build input, outside the public deposit

BP = "source_inputs/bp_cb_tip_removal/data_root/pcwde_supporting_data/repository_inputs/"
FLOW = "source_inputs/github_flow_integration/prs/"
REBUILD = "source_inputs/locus_tree_rebuild/"
FAMILIES = ("GH5_2", "PL1_8", "PL4_4")
CORRECTED_SYNTENY = {
    "RECONSTRUCTED_AFFECTED_ORDER.tsv": "source_inputs/synteny_corrected_placement/reconstruct/inputs/RECONSTRUCTED_AFFECTED_ORDER.tsv",
    "gca_to_species.tsv": "source_inputs/synteny_corrected_placement/reconstruct/stage/gene_order_full/gca_to_species.tsv",
    "gene_order_final.tsv": "source_inputs/synteny_corrected_placement/reconstruct/stage/final_gene_order/gene_order_final.tsv",
}
ACCEPTED_BLECTIN_SEQUENCES = "analysis/beta_uniform_transcript_membership/results/outputs/accepted_sequences.faa"

# Pruning-run outputs carried as files: the removal record, two method notes
# with revised tip counts, and the Figure S4 domain-figure inputs.
BP_CARRIED = (
    "deposit/base/sidecars/removed_cb_origin_bp_tips.tsv",
    "deposit/base/sidecars/origin_context_PL1_8_full_method.txt",
    "deposit/base/sidecars/origin_context_PL4_4_full_method.txt",
    "figures/figure_s4/Figure_S4.svg",
    "figures/figure_s4/render_qc.tsv",
)
# Domain-preparation records of the removed transcripts (FASTA ids / domtbl first column).
BP_DOMAIN_RECORDS = {
    "PL1_8": {"BP_TRINITY_DN116071_c0_g1_i2.p1"},
    "PL4_4": {"BP_TRINITY_DN7234_c0_g1_i2.p1", "BP_TRINITY_DN78954_c0_g1_i1.p1"},
}
# Rebuild Figure 4 audit files nothing reads.
LOCUS_UNUSED = ("figure6/accepted_old_vs_current_row_audit.tsv", "figure6/semantic_preservation_gate_summary.tsv")
SAVED_TABLES = tuple(f"Table_S{n}{s}" for n in (2, 3, 4) for s in (".tsv", "_caption.txt"))


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


class Upstream:
    """Resolve upstream files and check (or record) their pinned hashes."""

    def __init__(self, workspace, record):
        self.workspace, self.record, self.used = workspace, record, {}
        self.pins = {}
        if not record:
            with SOURCES.open(newline="") as handle:
                self.pins = {r["upstream_path"]: r for r in csv.DictReader(handle, delimiter="\t")}

    def __call__(self, rel, step):
        path = source_path(self.workspace, rel)
        digest = sha(path)
        if not self.record and self.pins.get(rel, {}).get("sha256") != digest:
            raise SystemExit(f"upstream file missing from or changed against config/deposit_sources.tsv: {rel}")
        self.used[rel] = (step, digest)
        return path

    def write(self):
        with SOURCES.open("w", newline="") as handle:
            w = csv.writer(handle, delimiter="\t", lineterminator="\n")
            w.writerow(["step", "upstream_path", "sha256"])
            for rel, (step, digest) in sorted(self.used.items(), key=lambda x: (x[1][0], x[0])):
                w.writerow([step, rel, digest])


def drop_records(path, ids, fasta):
    lines = path.read_text(encoding="utf-8").splitlines(keepends=True)
    out, dropped, skipping = [], set(), False
    for line in lines:
        if fasta:
            if line.startswith(">"):
                rid = line[1:].split()[0]
                skipping = rid in ids
                if skipping:
                    dropped.add(rid)
            if not skipping:
                out.append(line)
        else:
            rid = line.split(maxsplit=1)[0] if line.strip() and not line.startswith("#") else ""
            if rid in ids:
                dropped.add(rid)
            else:
                out.append(line)
    if dropped != ids:
        raise SystemExit(f"{path}: expected to drop {sorted(ids)}, dropped {sorted(dropped)}")
    path.write_text("".join(out), encoding="utf-8")


def bp_cb_pruning(data_root, inputs, up):
    work = data_root.parent / "bp_cb_evidence"
    work.mkdir()
    py = sys.executable
    subprocess.run([py, "-I", str(HERE / "drop_rows.py"), str(data_root), str(work / "row_drops.tsv")], check=True, stdout=subprocess.DEVNULL)
    subprocess.run([py, "-I", str(HERE / "prune_trees.py"), str(REPO / "config/bp_cb_prune_plan.tsv"), str(data_root),
                    str(work / "prune_events.tsv")], check=True, stdout=subprocess.DEVNULL)
    shutil.rmtree(work)
    for family, ids in BP_DOMAIN_RECORDS.items():
        prep = inputs / "analysis/domain_preparation"
        drop_records(prep / f"{family}.faa", ids, fasta=True)
        drop_records(prep / f"{family}_legacy_trimmed.faa", ids, fasta=True)
        drop_records(prep / f"{family}_domtbl.txt", ids, fasta=False)
    for rel in BP_CARRIED:
        shutil.copyfile(up(BP + rel, "bp_cb_pruning"), inputs / rel)


def copy_table(inputs, up):
    src = up(FLOW + "01_placement_conformance/supporting_data/PLACEMENT_TABLE_FINAL.tsv", "copy_table")
    with src.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    fields = [f for f in rows[0] if f != "evidence"]  # internal run paths
    out = inputs / "shared/copy_placements.tsv"
    out.parent.mkdir(exist_ok=True)
    with out.open("w", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def computed_tables(inputs):
    for name in SAVED_TABLES:
        (inputs / "tables/accepted" / name).unlink()


def locus_trees(inputs, up):
    figs = inputs / "figures"
    for d in ("figure6", "figure_s3"):
        shutil.rmtree(figs / d)
        prefix = REBUILD + f"figure_inputs/{d}/"
        members = sorted(name for name in up.pins if name.startswith(prefix))
        if not members:
            raise SystemExit(f"missing pinned locus directory: {d}")
        for name in members:
            rel = Path(name[len(prefix):])
            if rel.is_absolute() or ".." in rel.parts:
                raise SystemExit("unsafe locus directory member")
            source = up(name, "locus_trees")
            (figs / d / rel).parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source, figs / d / rel)
    for rel in LOCUS_UNUSED:
        (figs / rel).unlink()
    # roster carrying the same-gene BK venom slot (2026-10-28)
    shutil.copyfile(up("source_inputs/fig4_matched_tips_alignment/inputs/approved_roster_projection.tsv", "locus_trees"),
                    figs / "figure6/approved_roster_projection.tsv")
    with up(REBUILD + "DISPLAY_TIP_CROSSWALK.tsv", "locus_trees").open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    fields = ["family", "display_tree_tip", "raw_tip", "species", "source_kind", "source_protein_id", "corrected_chromosome"]
    with (inputs / "shared/locus_tips.tsv").open("w", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def reader_locus(inputs, up):
    base = inputs / "deposit/base"
    with up(REBUILD + "figure_inputs/figure6/raw_tip_to_normalized_display_label.tsv", "reader_locus").open(newline="") as handle:
        labels = {(r["family"], r["raw_tree_tip"]): r for r in csv.DictReader(handle, delimiter="\t")}
    with up(REBUILD + "DISPLAY_TIP_CROSSWALK.tsv", "reader_locus").open(newline="") as handle:
        cross = list(csv.DictReader(handle, delimiter="\t"))
    for fam in FAMILIES:
        aln = up(REBUILD + f"trees/{fam}/{fam}.aln.faa", "reader_locus")
        shutil.copyfile(aln, base / f"alignments/complete_locus_{fam}.aln.faa")
        shutil.copyfile(up(REBUILD + f"inputs/{fam}_tree_input.faa", "reader_locus"), base / f"sequences/complete_locus_{fam}.faa")
        shutil.copyfile(up(REBUILD + f"trees/{fam}/{fam}.bacterial_rooted.treefile", "reader_locus"),
                        base / f"trees/complete_locus_{fam}.inference_with_bacterial_root.treefile")
        shutil.copyfile(up(REBUILD + f"figure_inputs/figure_s3/complete_locus_{fam}.rooted.treefile", "reader_locus"),
                        base / f"trees/complete_locus_{fam}.display_diplolepididae_rooted_ladderized.treefile")
        raw_ids = [line[1:].split()[0] for line in aln.read_text().splitlines() if line.startswith(">")]
        raw_to_display = {r["raw_tip"]: r["display_tree_tip"] for r in cross if r["family"] == fam}
        fam_labels = {k[1]: v for k, v in labels.items() if k[0] == fam}
        by_protein = {v["tree_protein_id"]: k for k, v in fam_labels.items()}
        rows = []
        for raw in raw_ids:
            display = raw_to_display.get(raw) or by_protein.get(raw) or (raw if raw in fam_labels else None)
            if display not in fam_labels:
                raise SystemExit(f"{fam}: no display label for alignment tip {raw}")
            lab = fam_labels[display]
            rows.append({"raw_tip_id": raw, "display_tree_tip": display, "source_protein_id": lab["tree_protein_id"],
                         "source_species_label": lab["canonical_scientific_name"], "display_label": lab["normalized_display_label"]})
        with (base / f"sidecars/complete_locus_{fam}_tip_labels.tsv").open("w", newline="") as handle:
            w = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(rows)


def gh1_counts(inputs, up):
    fig1 = inputs / "figures/figure1"
    src = FLOW + "07_combined_figure1/supporting_data/"
    with up(src + "GH1E72_CYNIPOID_COUNTS.tsv", "gh1_counts").open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    with (fig1 / "gh1e72_cynipoid_counts.tsv").open("w", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=["accession", "species", "gh1e72_locus_count", "na_reason"],
                           delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    shutil.copyfile(up(src + "internal_gh1e72_gene_counts.tsv", "gh1_counts"), fig1 / "internal_gh1e72_gene_counts.tsv")


def revised_label(rel):
    if rel == "accessory_domain_trees/architecture.tsv":
        return "uniform full-protein S4 domain architecture"
    if rel.startswith(("locus_trees/", "complete_locus_trees/")):
        return "2026-10-23 locus-tree rebuild input"
    if rel.startswith(("b_lectin_tree/", "b_lectin_domain_tree/")):
        return "CB-origin B. pallida tips removed"
    if rel.startswith("presence_counts/gh1e72") or rel.startswith("presence_counts/internal_gh1e72"):
        return "GH1_e72 counts"
    raise SystemExit(f"no revision label for changed figure input {rel}")


def tblastn_presence(inputs, up):
    src = up("source_inputs/presence_count_consistency/TBLASTN_FUNNEL.tsv", "tblastn_presence")
    with src.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    fields = ["family", "species", "final_state", "best_passing_hsp", "reason"]
    with (inputs / "shared/tblastn_presence.tsv").open("w", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows({**r, "family": "B-lectin" if r["family"] == "beta_lectin" else r["family"]} for r in rows)


RUN2A = "source_inputs/member_layer_locked_inputs/outputs/"
PLACE = "source_inputs/placement_methods_conformance/"
MEMBER_COLUMNS = ["species", "genome", "family", "protein_id", "unit_id", "unit_kind", "annotation_source", "placement_route",
                  "bk_chromosome", "placement_reason", "curated_scaffold", "unit_route", "unit_record_class",
                  "unit_coordinate_subject", "unit_coordinate_start", "unit_coordinate_end", "unit_architecture_n_to_c"]
ORDER_COLUMNS = ["genome", "protein_id", "scaffold", "rank", "position", "start_bp", "end_bp", "strand", "coordinate_status"]


def member_layer(inputs, up):
    """Frozen member layer: one row per gene copy, plus each three-family copy's ordered-gene positions.

    Ordered positions come from step 2a's locked final gene order (member_order.tsv);
    coordinates from the locked annotation where present, else the placement run's
    reconstructed-scaffold inputs. Placement fields are the accepted placement run's
    per-locus result, carried on each copy. B-lectin units are curated rows.
    """
    def table(path):
        with Path(path).open(newline="") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))
    order = table(up(RUN2A + "member_layer/member_order.tsv", "member_layer"))
    attr = {(r["genome"], r["protein_id"], r["scaffold"], r["rank"]): r
            for r in table(up(RUN2A + "member_attributes.tsv", "member_layer")) if r["start_bp"] not in ("", "nan")}
    start_end = {(r["genome"], r["scaffold"], r["protein_id"]): r
                 for r in table(up(PLACE + "MISSING_ORDER_PRODUCER_INPUT.tsv", "member_layer"))}
    start_only = {(r["genome"], r["scaffold"], r["protein_id"]): r
                  for r in table(up(PLACE + "RECONSTRUCTED_AFFECTED_ORDER.tsv", "member_layer"))}
    out_order = []
    for o in order:
        row = {k: o[k] for k in ("genome", "protein_id", "scaffold", "rank", "position")}
        a = attr.get((o["genome"], o["protein_id"], o["scaffold"], o["rank"]))
        b = start_end.get((o["genome"], o["scaffold"], o["protein_id"]))
        c = start_only.get((o["genome"], o["scaffold"], o["protein_id"]))
        if a:
            row.update(start_bp=a["start_bp"], end_bp=a["end_bp"], strand=a["strand"], coordinate_status="locked_annotation")
        elif b:
            row.update(start_bp=b["start_bp"], end_bp=b["end_bp"], strand="", coordinate_status="placement_input_start_end")
        elif c:
            row.update(start_bp=c["start_bp"], end_bp="", strand="", coordinate_status="placement_input_start_only")
        else:
            row.update(start_bp="", end_bp="", strand="", coordinate_status="none")
        out_order.append(row)
    beta_rows = {}
    for r in table(inputs / "tables/authority/table_s2/beta_units.tsv"):
        beta_rows.setdefault(r["source_record_id"], []).append(r)
    # A source record can carry several curated domain units; keep each field's distinct values.
    beta = {k: {f: ";".join(dict.fromkeys(r[f] for r in rows)) for f in rows[0]} for k, rows in beta_rows.items()}
    members = []
    for r in table(inputs / "shared/copy_placements.tsv"):
        base = {"species": r["species"], "genome": r["genome"], "family": r["family"], "annotation_source": r["annotation_source"],
                "placement_route": r["placement_route"], "bk_chromosome": r["bk_chromosome"], "placement_reason": r["reason"]}
        if r["family"] == "B-lectin":
            u = beta.get(r["protein_ids"])
            if u is None:
                raise SystemExit(f"curated B-lectin unit without beta_units evidence (source_record_id): {r['protein_ids']}")
            members.append({**base, "protein_id": r["protein_ids"], "unit_id": r["locus_protein_id"], "unit_kind": "curated_unit", "curated_scaffold": r["scaffold"],
                            "unit_route": u["route"], "unit_record_class": u["record_class"],
                            "unit_coordinate_subject": u["coordinate_subject"], "unit_coordinate_start": u["coordinate_start"],
                            "unit_coordinate_end": u["coordinate_end"], "unit_architecture_n_to_c": u["architecture_n_to_c"]})
        else:
            for pid in r["protein_ids"].split(","):
                members.append({**base, "protein_id": pid, "unit_kind": "rule"})
    (inputs / "members").mkdir(exist_ok=True)
    for name, columns, rows in (("member_table.tsv", MEMBER_COLUMNS, members), ("member_order.tsv", ORDER_COLUMNS, out_order)):
        with (inputs / "members" / name).open("w", newline="") as handle:
            w = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n", restval="")
            w.writeheader()
            w.writerows(rows)


def derived_copy_table(inputs):
    (inputs / "shared/copy_placements.tsv").unlink()  # make builds it from members/ (scripts/members/build_copy_table.py)


def derived_locus_tips(inputs):
    """Keep only the tree-intrinsic tip crosswalk; make builds the full tip table (scripts/members/build_locus_tips.py)."""
    with (inputs / "shared/locus_tips.tsv").open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    fields = ["family", "display_tree_tip", "raw_tip", "source_kind", "source_protein_id", "transcript_species"]
    with (inputs / "members/locus_tip_crosswalk.tsv").open("w", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows({**{k: r[k] for k in fields[:5]},
                     "transcript_species": r["species"] if r["source_kind"] == "transcript_assembly" else ""} for r in rows)
    (inputs / "shared/locus_tips.tsv").unlink()


def retired_figures(inputs):
    shutil.rmtree(inputs / "figures/figure4")  # synteny figure dropped (Adam, 2026-10-29)


# Figure input folders by stable figure id; numbers are applied only at assembly (config/figure_numbers.tsv).
STABLE_DIRS = {
    "figures/figure1": "figures/presence_counts", "figures/figure2": "figures/chromosome_counts",
    "figures/figure3": "figures/origin_trees", "figures/figure5": "figures/introns",
    "figures/figure6": "figures/locus_trees", "figures/figure7": "figures/b_lectin_tree",
    "figures/figure8": "figures/expression", "figures/figure_s2": "figures/accessory_domain_trees",
    "figures/figure_s3": "figures/complete_locus_trees", "figures/figure_s4": "figures/b_lectin_domain_tree",
    "tables/figure5_current": "tables/introns_current",
}


def stable_path(rel):
    """Map a submitted-deposit path (relative to repository_inputs) to its stable-name path."""
    for old, new in STABLE_DIRS.items():
        if rel == old or rel.startswith(old + "/"):
            return new + rel[len(old):]
    return rel


def stable_names(inputs):
    for old, new in STABLE_DIRS.items():
        (inputs / old).rename(inputs / new)


# Reader-deposit cleanup (Adam, 2026-10-29/30): inputs no build step reads, superseded saved
# Tables S6/S7, and the stale synteny-locus reader tables (the deposit builder now ships the
# build's gene-copy family_membership.tsv and chromosome_placement_evidence.tsv).
READER_CLEANUP = (
    "figures/chromosome_counts/ae_pcwde_chr_summary.tsv",
    "figures/chromosome_counts/aqp_current_pcwde_protein_map.tsv",
    "figures/chromosome_counts/beta_units.tsv",
    "figures/chromosome_counts/dl_pcwde_chr_summary.tsv",
    "figures/chromosome_counts/new_genome_family_summary.tsv",
    "figures/chromosome_counts/new_genome_ordinary_cells.tsv",
    "figures/chromosome_counts/validated_raw_locus_rescue.tsv",
    "figures/shared/BK_chromosome_mapping.json",
    "figures/shared/gca_to_species.tsv",
    "figures/locus_trees/all_family_rna_evidence_audit.tsv",
    "tables/style/Table_S6.pdf",
    "tables/authority/synteny_current/current_chromosome_cynog_dict.tsv",
    "tables/authority/synteny_current/locus_detail_source.tsv",
    "tables/authority/synteny_current/phylogenetic_order_reference.tsv",
    "tables/authority/synteny_current/producer_public_dictionary.tsv",
    "tables/authority/tables_s5_s6_current/AUTHORITY_MANIFEST.sha256",
    "tables/authority/tables_s5_s6_current/expected/per_intron_observations.tsv",
    "tables/authority/tables_s5_s6_current/expected/recurrent_position_summary.tsv",
    "tables/authority/tables_s5_s6_current/intron_inputs/CBM13_intron_positions.tsv",
    "tables/authority/tables_s5_s6_current/intron_inputs/GH5_2_intron_positions.tsv",
    "tables/authority/tables_s5_s6_current/intron_inputs/PL1_8_intron_positions.tsv",
    "tables/authority/tables_s5_s6_current/intron_inputs/PL4_4_intron_positions.tsv",
    "tables/authority/tables_s5_s6_current/intron_inputs/extracted/CBM13_intron_metadata.tsv",
    "tables/authority/tables_s5_s6_current/intron_inputs/extracted/GH5_2_intron_metadata.tsv",
    "tables/authority/tables_s5_s6_current/intron_inputs/extracted/PL1_8_intron_metadata.tsv",
    "tables/authority/tables_s5_s6_current/intron_inputs/extracted/PL4_4_intron_metadata.tsv",
    "tables/accepted/Table_S6.tsv",
    "tables/accepted/Table_S6_caption.txt",
    "tables/accepted/Table_S7_catalytic.tsv",
    "tables/accepted/Table_S7_selection.tsv",
    "deposit/base/data/family_membership.tsv",
    "deposit/base/data/chromosome_placement_evidence.tsv",
    # Figure 2 and S1 now prune the pinned complete origin trees during make.
    "figures/origin_trees/GH5_2_pruned.treefile",
    "figures/origin_trees/GH5_2_s1.treefile",
    "figures/origin_trees/PL1_8_pruned.treefile",
    "figures/origin_trees/PL1_8_s1.treefile",
    "figures/origin_trees/PL4_4_pruned.treefile",
    "figures/origin_trees/PL4_4_s1.treefile",
)


def reader_cleanup(inputs):
    for rel in READER_CLEANUP:
        (inputs / rel).unlink()
    for rel in sorted({str(Path(r).parent) for r in READER_CLEANUP}, key=len, reverse=True):
        d = inputs / rel
        if d.is_dir() and not any(d.iterdir()):
            d.rmdir()


def source_table_s1(inputs, up):
    """Replace the saved S1 table with the frozen catalog cells and reports."""
    run = 'source_inputs/table_s1_rnaseq_carry_forward/'
    dst = inputs / 'tables/authority/table_s1_catalog'
    dst.mkdir(parents=True, exist_ok=True)
    shutil.copy2(up(run + 'S1_CATALOG_BUSCO_PIN.tsv', 's1_catalog'), dst / 'catalog_busco.tsv')
    # Raw reports contain private cluster paths. Their exact SHA-256 identities
    # are public in this manifest; the byte-identical copies stay in run evidence.
    shutil.copy2(up(run + 'S1_REPORT_PINS.tsv', 's1_catalog'), dst / 'report_pins.tsv')
    (inputs / 'tables/accepted/Table_S1.tsv').unlink()


NEW_TABLE_SOURCES = {
    "catalytic_clustering/ANALYSIS.json": "source_inputs/mbe_rev_catalytic_clustering/ANALYSIS.json",
    "locus_conflict/AU_RESULTS.tsv": "source_inputs/locus_tree_conflict_tests/AU_RESULTS.tsv",
    "locus_conflict/SCF_RESULTS.tsv": "source_inputs/locus_tree_conflict_tests/SCF_RESULTS.tsv",
    "locus_conflict/LONG_BRANCH_SENSITIVITY.tsv": "source_inputs/locus_tree_conflict_tests/LONG_BRANCH_SENSITIVITY.tsv",
    "di_sort_limits/SORT_LIMITS.md": "source_inputs/di_inquiline_sort_limits/SORT_LIMITS.md",
    "di_sort_limits/COORDINATOR_ACCEPTANCE_NOTE.md": "source_inputs/di_inquiline_sort_limits/COORDINATOR_ACCEPTANCE_NOTE.md",
}


def source_new_tables(inputs, up):
    """Carry accepted analysis summaries with SHA256 checked upstream identity."""
    root = inputs / "tables/authority/supplement_new"
    for dest, rel in NEW_TABLE_SOURCES.items():
        target = root / dest
        target.parent.mkdir(parents=True, exist_ok=True)
        source = (Path(__file__).resolve().parents[2] / "scripts/analysis/catalytic/clustering/LOCKED_ANALYSIS.json") if dest == "catalytic_clustering/ANALYSIS.json" else up(rel, "supplement_new_tables")
        # Public projection: provenance fields in two sources name absolute
        # cluster paths. The upstream SHA-256 above is checked on the exact
        # original; the deposit copy keeps workspace-relative paths only.
        payload = source.read_bytes()
        # Workspace paths become pcwde_paper/...; other project paths become
        # project-root-relative (e.g. good_data/...).
        for ws in {up.workspace, up.workspace.resolve()}:
            payload = payload.replace(str(ws).rstrip("/").encode() + b"/", b"pcwde_paper/")
            payload = payload.replace(str(ws.parents[2]).rstrip("/").encode() + b"/", b"")
        if any(t in payload for t in (b"/" + b"carc/", b"/" + b"easley/", b"/" + b"users/", b"/" + b"root/")):
            raise SystemExit(f"private absolute path remains after projection: {dest}")
        target.write_bytes(payload)
        shutil.copystat(source, target)
def retire_built_introns(inputs):
    """Remove saved intron data now generated from pinned annotation sources."""
    for rel in ("figures/introns", "tables/introns_current", "deposit/base/data/introns"):
        path = inputs / rel
        if not path.is_dir():
            raise SystemExit(f"missing saved intron input family: {rel}")
        shutil.rmtree(path)
    for name in ("Table_S5.tsv", "Table_S5_caption.txt"):
        saved = inputs / "tables/accepted" / name
        if not saved.is_file():
            raise SystemExit(f"missing saved {name}")
        saved.unlink()


INTRON_FROZEN = REPO / "config/intron_frozen_sources.tsv"


def intron_sources(inputs, up):
    """Stage only checksum-pinned small layer-one inputs from the upstream run."""
    dst = inputs / "analysis/introns/frozen"
    if dst.exists():
        raise SystemExit(f"refusing existing intron frozen inputs: {dst}")
    project = up.workspace.parents[2]
    with INTRON_FROZEN.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 94 or len({r["deposit_relative_path"] for r in rows}) != len(rows):
        raise SystemExit("intron frozen roster is incomplete or duplicated")
    for r in rows:
        source_rel, target_rel = Path(r["project_path"]), Path(r["deposit_relative_path"])
        if (source_rel.is_absolute() or target_rel.is_absolute() or
            ".." in source_rel.parts or ".." in target_rel.parts or
            target_rel.parts[:3] != ("analysis", "introns", "frozen")):
            raise SystemExit("unsafe intron frozen path")
        src, out = source_path(project, source_rel), inputs / target_rel
        if not src.is_file() or src.is_symlink() or src.stat().st_size != int(r["bytes"]) or sha(src) != r["sha256"]:
            raise SystemExit(f"intron frozen input drift: {source_rel}")
        out.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(src, out)
        if sha(out) != r["sha256"]:
            raise SystemExit(f"intron frozen copy mismatch: {target_rel}")



def write_manifests(inputs, base_inputs):
    """Deposit inventory, and the figure inventory with labels kept from the submitted deposit."""
    figs = inputs / "figures"
    with (base_inputs / "figures/INPUT_MANIFEST.tsv").open(newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    head, body = rows[0], rows[1:]
    actual = {p.relative_to(figs).as_posix() for p in figs.rglob("*") if p.is_file()} - {"INPUT_MANIFEST.tsv", "figure_s5/Figure_S5.pdf"}
    renamed = []
    for r in body:
        new = stable_path("figures/" + r[1]).removeprefix("figures/")
        if new in actual:
            if sha(figs / new) != sha(base_inputs / "figures" / r[1]):
                r[2] = revised_label(new)
            renamed.append([new.split("/")[0], new, r[2], r[3]])
    body = renamed
    for rel in sorted(actual - {r[1] for r in body}):
        body.append([rel.split("/")[0], rel, revised_label(rel), "0"])
    for r in body:
        r[3] = str((figs / r[1]).stat().st_size)
    with (figs / "INPUT_MANIFEST.tsv").open("w", newline="") as handle:
        w = csv.writer(handle, delimiter="\t", lineterminator="\n")
        w.writerow(head)
        w.writerows(body)
    write_top_manifest(inputs)


def write_top_manifest(inputs):
    """Inventory every staged file after a pinned overlay changes the deposit."""
    manifest = inputs / "INPUT_MANIFEST.tsv"
    with manifest.open("w", newline="") as handle:
        w = csv.writer(handle, delimiter="\t", lineterminator="\n")
        w.writerow(["relative_path", "size_bytes", "source_class"])
        files = sorted((p for p in inputs.rglob("*") if p.is_file() and p != manifest and
                        not p.relative_to(inputs).as_posix().startswith("cazy_inputs/")),
                       key=lambda p: (p.relative_to(inputs).as_posix().startswith("tables/authority/synteny_corrected/"),
                                      p.relative_to(inputs).as_posix()))
        with CAZY_SOURCES.open(newline="") as source_handle:
            files += [inputs.parent / row["deposit_relative_path"]
                      for row in csv.DictReader(source_handle, delimiter="\t")]
        for p in files:
            rel = p.relative_to(inputs).as_posix()
            if rel not in EXCLUDED:
                source_class = ("frozen_cazy_scan_input" if rel.startswith("cazy_inputs/") else
                                "corrected_synteny_authority" if rel.startswith("tables/authority/synteny_corrected/") else
                                "accepted_public_intermediate")
                w.writerow([rel, p.stat().st_size, source_class])


def rescue_overlay(inputs, workspace):
    """Carry the reviewed unplaced PL1_8 locus into rebuilt staging."""
    source = workspace / RESCUE_SOURCE
    if sha(source) != RESCUE_SHA256:
        raise SystemExit(f"rescue source changed: {source}")
    target = inputs / "shared/fig5_validated_raw_locus_rescue.tsv"
    shutil.copyfile(source, target)


def cazy_inputs(inputs, workspace, expected_rows=148, require_scan_pins=True):
    """Stage the frozen scan sources and the one portable provenance projection."""
    if (inputs / "cazy_inputs").exists():
        raise SystemExit("unexpected pre-existing CAZy input tree")
    project_root = workspace.resolve().parents[2]
    with CAZY_SOURCES.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != expected_rows or len({r["deposit_relative_path"] for r in rows}) != expected_rows:
        raise SystemExit(f"CAZy source map must declare {expected_rows} unique files")
    declared = set()
    for row in rows:
        src_rel = Path(row["source_relative_path"])
        dst_rel = Path(row["deposit_relative_path"])
        if (src_rel.is_absolute() or ".." in src_rel.parts or dst_rel.is_absolute() or
                ".." in dst_rel.parts or dst_rel.parts[:2] != ("repository_inputs", "cazy_inputs")):
            raise SystemExit("unsafe CAZy source path")
        source = source_path(project_root, src_rel)
        dest = inputs.parent / dst_rel
        if not source.is_file() or source.is_symlink() or dest.exists():
            raise SystemExit(f"missing, linked or duplicate CAZy source: {src_rel}")
        if row["derivation"]:
            if row["derivation"] != "portable_provenance_v1":
                raise SystemExit("unknown CAZy derivation")
            if source.stat().st_size != int(row["original_bytes"]) or sha(source) != row["original_sha256"]:
                raise SystemExit("CAZy provenance source drift")
            spec = importlib.util.spec_from_file_location("derive_cazy_provenance", REPO / "scripts/release/derive_cazy_provenance.py")
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            payload = module.derive(source)
        else:
            if row["original_sha256"] or row["original_bytes"]:
                raise SystemExit("unexpected CAZy derivation fields")
            if source.stat().st_size != int(row["bytes"]) or sha(source) != row["sha256"]:
                raise SystemExit(f"CAZy source drift: {src_rel}")
            payload = None
        dest.parent.mkdir(parents=True, exist_ok=True)
        if payload is None:
            shutil.copyfile(source, dest)
        else:
            dest.write_bytes(payload)
        if dest.stat().st_size != int(row["bytes"]) or sha(dest) != row["sha256"]:
            raise SystemExit(f"CAZy staged bytes mismatch: {dst_rel}")
        declared.add(dst_rel.as_posix())
    scan_paths = set()
    with (REPO / "config/cazy_scans.tsv").open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            scan_paths.update(row[k] for k in ("proteome_path", "scan_path", "mapping_source") if row[k])
    with (REPO / "config/cazy_processed_proteomes.tsv").open(newline="") as handle:
        scan_paths.update(row["processed_path"] for row in csv.DictReader(handle, delimiter="\t"))
    if require_scan_pins and not scan_paths <= declared:
        raise SystemExit(f"CAZy scan/processed pins absent from source map: {sorted(scan_paths - declared)}")
    print(f"CAZY_STAGE_PASS files={len(declared)} pins={len(scan_paths)}")


def corrected_synteny(inputs, up):
    target = inputs / "tables/authority/synteny_corrected"
    target.mkdir(parents=True, exist_ok=True)
    for name, rel in sorted(CORRECTED_SYNTENY.items()):
        shutil.copyfile(up(rel, "synteny_corrected_placement"), target / name)


def synteny_existing_identity(inputs, up, workspace):
    """Bind accepted copies to already ordered models using pinned identity sources."""
    project = workspace.resolve().parents[2]
    if workspace.name != "pcwde_paper" or project.name != "adam":
        raise SystemExit(f"unexpected paper/project root: {workspace}")
    accepted = up(ACCEPTED_BLECTIN_SEQUENCES, "synteny_existing_identity")
    subprocess.run([
        sys.executable, str(REPO / "scripts/analysis/synteny/resolve_existing_members.py"),
        "--project", str(project),
        "--members", str(inputs / "members/member_table.tsv"),
        "--member-order", str(inputs / "members/member_order.tsv"),
        "--order", str(inputs / "tables/authority/synteny_corrected/gene_order_final.tsv"),
        "--accepted-fasta", str(accepted),
        "--accepted-sha256", "274eb8e8839d2edd89cddbecdb6467784dfa580eecfe4a47bfd5c8632255f761",
        "--sources", str(REPO / "config/intron_sources.tsv"),
        "--genome-map", str(REPO / "config/intron_genome_sources.tsv"),
        "--order-sources", str(REPO / "config/synteny_existing_sources.tsv"),
        "--out", str(inputs / "members/synteny_existing_identity.tsv"),
    ], check=True)


def synteny_existing_producer(inputs, workspace):
    """Regenerate native loci and recovery from the frozen source layer."""
    project = workspace.resolve().parents[2]
    if workspace.name != "pcwde_paper" or project.name != "adam":
        raise SystemExit(f"unexpected paper/project root: {workspace}")
    work = inputs.parent.parent / "synteny_existing_work"
    if work.exists():
        raise SystemExit(f"refusing to reuse synteny producer work: {work}")
    work.mkdir()
    source_manifest = REPO / "config/synteny_existing_sources.tsv"
    identity = inputs / "members/synteny_existing_identity.tsv"
    backend = work / "backend"
    recovery = work / "recovery"
    common = ["--project", str(project), "--sources", str(source_manifest)]
    subprocess.run([sys.executable, str(REPO / "scripts/analysis/synteny/produce_existing_backend.py"),
                    *common, "--identity", str(identity), "--out", str(backend)], check=True)
    subprocess.run([sys.executable, str(REPO / "scripts/analysis/synteny/produce_existing_recovery.py"),
                    *common, "--identity", str(identity), "--backend", str(backend),
                    "--out", str(recovery)], check=True)
    target = inputs / "tables/authority/synteny_corrected"
    produced = {
        "blocks_long.tsv": recovery / "build/table_s06_synteny_blocks_long.tsv",
        "block_summary.tsv": recovery / "build/table_s06_synteny_block_summary.tsv",
        "current_region_status_complete.tsv": recovery / "current_region_status_complete.tsv",
    }
    observed = {name: sha(path) for name, path in produced.items()}
    pinfile = REPO / "config/synteny_derived_pins.tsv"
    if os.environ.get("PCWDE_SYNTENY_PIN_MODE") == "generate":
        with (work / "DERIVED_PINS.tsv").open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["name", "sha256"])
            writer.writerows(sorted(observed.items()))
    else:
        with pinfile.open(newline="") as handle:
            expected = {r["name"]: r["sha256"] for r in csv.DictReader(handle, delimiter="\t")}
        if observed != expected:
            raise SystemExit(f"derived synteny pins changed: observed={observed}, expected={expected}")
    for name, path in produced.items():
        shutil.copyfile(path, target / name)


def validate_projected_consumer(data_root):
    """Check the projected-input contract before consuming a staged public root."""
    release_dir = Path(__file__).resolve().parents[1] / "release"
    sys.path.insert(0, str(release_dir))
    from data_root import validate_projected_inputs
    return validate_projected_inputs(data_root)


def validate_root_form(data_root):
    """Validate a data root by its own contract: projected manifest or raw source pins."""
    release_dir = Path(__file__).resolve().parents[1] / "release"
    sys.path.insert(0, str(release_dir))
    from data_root import validate_logical_inputs
    projected = data_root / "PROJECTED_INPUTS.tsv"
    if projected.exists() or projected.is_symlink():
        validate_projected_consumer(data_root)
        return "projected"
    validate_logical_inputs(REPO, data_root)
    return "raw"


def compare(a, b):
    # Both roots must have the same form, each validated by that form's contract,
    # before bytes are compared; a raw root is never compared with a projected one.
    form_a, form_b = validate_root_form(a), validate_root_form(b)
    if form_a != form_b:
        raise SystemExit(f"STAGE_COMPARE_FORM_MISMATCH {form_a} vs {form_b}")
    print(f"STAGE_FORM_PASS form={form_a}")
    def files(root):
        return {p.relative_to(root).as_posix(): p for p in root.rglob("*") if p.is_file()}
    fa, fb = files(a / "repository_inputs"), files(b / "repository_inputs")
    diffs = sorted(set(fa) ^ set(fb)) + sorted(k for k in set(fa) & set(fb) if sha(fa[k]) != sha(fb[k]))
    for d in diffs:
        print("DIFFERS", d)
    print(f"{'STAGE_COMPARE_PASS' if not diffs else 'STAGE_COMPARE_FAIL'} files={len(fa)} differing={len(diffs)}")
    return not diffs


def pin(data_root):
    """Refresh config/logical_inputs.tsv: new hashes/sizes, rows added or removed with the deposit's files."""
    inputs = data_root / "repository_inputs"
    present = {f"repository_inputs/{p.relative_to(inputs).as_posix()}" for p in inputs.rglob("*") if p.is_file()}
    path = REPO / "config/logical_inputs.tsv"
    lines = path.read_text().splitlines()
    expected_columns = 7
    ids, paths = set(), set()
    for line_number, line in enumerate(lines[1:], 2):
        fields = line.split("\t")
        if len(fields) != expected_columns:
            raise SystemExit(f"malformed logical input row {line_number}: {len(fields)} columns, expected {expected_columns}")
        if fields[0] in ids or fields[1] in paths:
            raise SystemExit(f"duplicate logical input key or path at row {line_number}: {fields[0]}")
        ids.add(fields[0]); paths.add(fields[1])
    out, seen = [lines[0]], set()
    for line in lines[1:]:
        f = line.split("\t")
        if f[1] not in present:
            print("unpinned", f[1])
            continue
        seen.add(f[1])
        f[2], f[3] = sha(data_root / f[1]), str((data_root / f[1]).stat().st_size)
        out.append("\t".join(f))
    for rel in sorted(present - seen):
        print("pinned", rel)
        out.append("\t".join([f"repository_input::{rel.removeprefix('repository_inputs/')}", rel, sha(data_root / rel),
                              str((data_root / rel).stat().st_size), "renderer_input" if "/figures/" in rel else "scientific_input",
                              "revised_deposit_stage", "verify;manuscript;data-deposit"]))
    path.write_text("\n".join(out) + "\n")
    print(f"PIN_PASS logical_inputs={len(out) - 1}")


def chromosome_gain_inputs(inputs, workspace):
    """Carry the checksum-pinned placement evidence into the public deposit."""
    project_root = workspace.resolve().parents[2]
    target = inputs / "chromosome_gain"
    target.mkdir()
    with CHROMOSOME_GAIN_SOURCES.open(newline="") as handle:
        records = list(csv.DictReader(handle, delimiter="\t"))
    if len(records) != 59 or len({r["role"] for r in records}) != len(records):
        raise SystemExit("chromosome-gain source manifest is incomplete")
    for record in records:
        relative = Path(record["project_path"])
        role = record["role"]
        if relative.is_absolute() or ".." in relative.parts or "/" in role or role in {"", ".", ".."}:
            raise SystemExit("unsafe chromosome-gain source declaration")
        source = source_path(project_root, relative)
        if not source.is_file() or source.is_symlink() or source.stat().st_size != int(record["bytes"]) or sha(source) != record["sha256"]:
            raise SystemExit(f"chromosome-gain upstream mismatch: {role}")
        packaged = target / role
        if role in {"D_kuriphilus_helixer.gff3", "GCA_038277395.1_helixer.gff3"}:
            # Helixer stamped a local model path into one provenance comment.
            # Annotation rows remain byte-identical and in the same order.
            removed = 0
            with source.open("rb") as original, packaged.open("wb") as public:
                for line in original:
                    if line.startswith(b"# fdd5a925728cf56e2470d96beeeb797b  "):
                        removed += 1
                    else:
                        public.write(line)
            if removed != 1:
                raise SystemExit(f"unexpected Helixer provenance projection: {role} ({removed})")
        else:
            shutil.copyfile(source, packaged)
        if packaged.stat().st_size != int(record["deposit_bytes"]) or sha(packaged) != record["deposit_sha256"]:
            raise SystemExit(f"chromosome-gain deposit projection mismatch: {role}")


def s4_architecture_inputs(inputs, up):
    """Carry the uniformly annotated S4 layer and its frozen source evidence."""
    run = "source_inputs/figure_s4_s6_legibility/work/"
    files = {
        "figures/accessory_domain_trees/architecture.tsv": "s4_architecture.tsv",
        "analysis/s4_architecture/full_sequence_mapping.tsv": "s4_full_sequence_mapping_preflight.tsv",
        "analysis/s4_architecture/full_sequences.faa": "s4_unique_full_sequences.faa",
        "analysis/s4_architecture/pfam_hits.domtblout": "s4_uniform_pfam_hits.domtblout",
    }
    for target, source in files.items():
        dest = inputs / target
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(up(run + source, "s4_architecture"), dest)


SUPPLEMENT_FACT_KIND = re.compile(r"^supplement fact source ([A-Za-z0-9_.]+) \(class [CU]\)$")


def supplement_fact_inputs(inputs, workspace):
    """Copy the pinned source reports of supplement prose facts into supplement_facts/."""
    project_root = workspace.resolve().parents[2]
    with LOCKED_INPUTS.open(newline="") as handle:
        records = [r for r in csv.DictReader(handle, delimiter="\t") if r["genome"] == "supplement"]
    target = inputs / "supplement_facts"
    target.mkdir()
    roles = set()
    for record in records:
        match = SUPPLEMENT_FACT_KIND.match(record["kind"])
        relative = Path(record["project_path"])
        if not match or relative.is_absolute() or ".." in relative.parts or match.group(1) in roles:
            raise SystemExit(f"unsafe or duplicate supplement fact source declaration: {record['kind']}")
        role = match.group(1)
        roles.add(role)
        source = source_path(project_root, relative)
        if not source.is_file() or source.stat().st_size != int(record["bytes"]) or sha(source) != record["sha256"]:
            raise SystemExit(f"supplement fact source mismatch: {role}")
        payload = source.read_bytes()
        # Public projection: run logs and BUSCO parameters name absolute cluster
        # paths; they become project-root-relative. Facts are not touched.
        for root in {project_root, project_root.resolve()}:
            payload = payload.replace(str(root).rstrip("/").encode() + b"/", b"")
        if any(t in payload for t in (b"/" + b"carc/", b"/" + b"easley/", b"/" + b"users/", b"/" + b"root/")):
            raise SystemExit(f"private absolute path remains after projection: {role}")
        if role == "current_proteomes.tsv":
            # The source digest is checked above and retained in upstream_locked_inputs.
            # The public roster carries only the nine scientific membership rows.
            text = payload.decode("utf-8")
            lines = list(csv.reader(text.splitlines(), delimiter="\t", strict=True))
            expected = ["species_code", "annotation_id", "proteome_path", "input_proteins", "input_sha256"]
            if not lines or lines[0] != expected or len(lines) != 10 or any(len(row) != len(expected) for row in lines[1:]):
                raise SystemExit("supplement proteome source schema or row count differs from pinned record")
            if any(not re.fullmatch(r"[0-9a-f]{64}", row[-1]) for row in lines[1:]):
                raise SystemExit("supplement proteome source digest format differs from pinned record")
            payload = ("\n".join("\t".join(row[:-1]) for row in lines) + "\n").encode("utf-8")
        (target / role).write_bytes(payload)
    if len(roles) != 9:
        raise SystemExit(f"expected 9 supplement fact sources, found {len(roles)}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", type=Path)
    parser.add_argument("--workspace", type=Path)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--record", action="store_true", help="write config/deposit_sources.tsv from the files used")
    parser.add_argument("--compare", type=Path, nargs=2)
    parser.add_argument("--pin", type=Path)
    parser.add_argument("--refresh-top-manifest", type=Path,
                        help="regenerate the top-level inventory after adding a pinned input")
    parser.add_argument("--chromosome-gain-only", action="store_true",
                        help="add pinned chromosome-gain evidence to a current staged data root")
    parser.add_argument("--supplement-facts-only", action="store_true",
                        help="add pinned supplement fact sources to a current staged data root")
    args = parser.parse_args()
    if args.pin:
        return pin(args.pin)
    if args.refresh_top_manifest:
        return write_top_manifest(args.refresh_top_manifest / "repository_inputs")
    if args.compare:
        sys.exit(0 if compare(*args.compare) else 1)
    if args.out.exists():
        raise SystemExit(f"refusing existing output: {args.out}")
    base_inputs = args.base / "repository_inputs"
    shutil.copytree(args.base, args.out)
    inputs = args.out / "repository_inputs"
    if args.chromosome_gain_only:
        chromosome_gain_inputs(inputs, args.workspace)
        manifest = inputs / "INPUT_MANIFEST.tsv"
        with manifest.open("a", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            for path in sorted((inputs / "chromosome_gain").iterdir()):
                writer.writerow([path.relative_to(inputs).as_posix(), path.stat().st_size,
                                 "accepted_public_intermediate"])
        print("CHROMOSOME_GAIN_STAGE_PASS", 59)
        return
    if args.supplement_facts_only:
        supplement_fact_inputs(inputs, args.workspace)
        manifest = inputs / "INPUT_MANIFEST.tsv"
        with manifest.open("a", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            for path in sorted((inputs / "supplement_facts").iterdir()):
                writer.writerow([path.relative_to(inputs).as_posix(), path.stat().st_size,
                                 "accepted_public_intermediate"])
        print("SUPPLEMENT_FACTS_STAGE_PASS", len(list((inputs / "supplement_facts").iterdir())))
        return
    (inputs / "figures/figure_s5/Figure_S5.pdf").unlink(missing_ok=True)
    if (inputs / "figures/figure_s5").exists():
        (inputs / "figures/figure_s5").rmdir()
    up = Upstream(args.workspace, args.record)
    bp_cb_pruning(args.out, inputs, up)
    copy_table(inputs, up)
    computed_tables(inputs)
    locus_trees(inputs, up)
    reader_locus(inputs, up)
    gh1_counts(inputs, up)
    tblastn_presence(inputs, up)
    member_layer(inputs, up)
    derived_copy_table(inputs)
    derived_locus_tips(inputs)
    retired_figures(inputs)
    stable_names(inputs)
    reader_cleanup(inputs)
    source_table_s1(inputs, up)
    source_new_tables(inputs, up)
    corrected_synteny(inputs, up)
    synteny_existing_identity(inputs, up, args.workspace)
    synteny_existing_producer(inputs, args.workspace)
    chromosome_gain_inputs(inputs, args.workspace)
    cazy_inputs(inputs, args.workspace)
    supplement_fact_inputs(inputs, args.workspace)
    rescue_overlay(inputs, args.workspace)
    s4_architecture_inputs(inputs, up)
    retire_built_introns(inputs)
    intron_sources(inputs, up)
    write_manifests(inputs, base_inputs)
    from retire_audit_inputs import retire
    retire(inputs)
    from public_scope import project
    project(inputs)
    if args.record:
        up.write()
    print(f"STAGE_PASS out={args.out} upstream_files={len(up.used)}")


if __name__ == "__main__":
    main()
