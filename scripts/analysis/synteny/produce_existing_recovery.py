#!/usr/bin/env python3
"""Adapt eligible native loci to the accepted final-49 recovery calculation."""
import argparse
import csv
import hashlib
import importlib.util
import shutil
from collections import defaultdict
from pathlib import Path

import pandas as pd

from produce_existing_backend import pinned_sources


BETA = "beta_lectin/GH30_lineage"
RECALC_SHA = "01f7c278750e75e6e17b6efa2f80f49e8b6e8f8b445988ec771313f1c90f1673"


def read(path):
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    obj = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj


def metadata(identity, authority, loci, fields, panel):
    identity_by_id = {r.accepted_id: r for r in identity.itertuples(index=False)}
    if len(identity_by_id) != len(identity):
        raise ValueError("duplicate accepted identity")
    native = {}
    species = defaultdict(set)
    for row in loci.itertuples(index=False):
        species[row.genome].add(row.species)
        for protein in row.protein_ids.split(","):
            if protein:
                if protein in native:
                    raise ValueError(f"duplicate backend protein: {protein}")
                native[protein] = row
    if any(len(names) != 1 for names in species.values()):
        raise ValueError("backend species label is not unique per genome")
    panel_species = dict(zip(panel.accession, panel.species.str.replace(" ", "_", regex=False)))
    if not set(identity.genome).issubset(panel_species):
        raise ValueError("accepted genome absent from final panel")
    out = []
    for entry in authority.itertuples(index=False):
        source = entry.source_protein_id
        if source not in identity_by_id:
            raise ValueError(f"authority member absent from identity disposition: {source}")
        ident = identity_by_id[source]
        if ident.family != entry.family_id:
            raise ValueError(f"identity/authority family conflict: {source}")
        builder = ident.order_id
        locus = native.get(builder) if builder else None
        if bool(builder) != bool(locus):
            raise ValueError(f"eligible/backend identity mismatch: {source}")
        family = BETA if ident.family == "B-lectin" else ident.family
        if locus and (locus.family != ("CBM13" if family == BETA else family)
                      or locus.genome != ident.genome):
            raise ValueError(f"backend family/genome mismatch: {source}")
        label = next(iter(species[ident.genome])) if species[ident.genome] else panel_species[ident.genome]
        out.append({**dict.fromkeys(fields, ""), "family": family,
                    "accepted_locus_id": source, "species": label,
                    "source_type": "genome", "record_class": "genome_locus",
                    "builder_target_id": builder, "chromosome": entry.reference_block,
                    "coordinate_subject": entry.source_scaffold,
                    "placement_state": entry.placement_class,
                    "projection_class": entry.placement_route,
                    "unassigned_reason": "no_existing_order_identity" if not builder else "",
                    "accepted_aliases": source})
    if len(out) != len(identity) or len({r["accepted_locus_id"] for r in out}) != len(out):
        raise ValueError("accepted authority coverage incomplete")
    return pd.DataFrame(out, columns=fields)


def membership(recalc, blocks, loci, metadata, fields):
    locus_meta = loci.set_index("locus_id").to_dict("index")
    block_meta = blocks.set_index("locus_id").to_dict("index")
    if set(locus_meta) != set(block_meta):
        raise ValueError("native locus/block keys differ")
    bound = metadata.loc[metadata.builder_target_id.ne("")]
    if bound.builder_target_id.duplicated().any():
        raise ValueError("two accepted units use one native target")
    identities = bound.set_index("builder_target_id").to_dict("index")
    rows, seen = [], set()
    for locus_id, source in sorted(locus_meta.items()):
        block = block_meta[locus_id]
        family = BETA if source["family"] == "CBM13" else source["family"]
        proteins = [p for p in source["protein_ids"].split(",") if p]
        if len(proteins) != int(source["n_genes"]):
            raise ValueError(f"native locus multiplicity changed: {locus_id}")
        for protein in dict.fromkeys(proteins):
            if protein not in identities or protein in seen:
                raise ValueError(f"missing/repeated native accepted identity: {protein}")
            seen.add(protein)
            info = identities[protein]
            if info["family"] != family:
                raise ValueError(f"target family conflict: {protein}")
            row = {field: info.get(field, "") for field in fields}
            row.update(builder_target_id=protein, builder_locus_id=locus_id,
                       block_id=recalc.block_key(type("Block", (), {"status": block["status"],
                           "family": source["family"], "locus_id": locus_id,
                           "block_id": block["block_id"]})(), True),
                       assignment_status="assigned", block_status=block["status"])
            rows.append(row)
    if seen != set(identities):
        raise ValueError("backend/identity target sets differ")
    for info in metadata.loc[metadata.builder_target_id.eq("")].to_dict("records"):
        row = {field: info.get(field, "") for field in fields}
        row.update(builder_target_id="", builder_locus_id="", block_id="",
                   assignment_status="unassigned", block_status="unassigned")
        rows.append(row)
    return pd.DataFrame(rows, columns=fields).sort_values(
        ["family", "species", "accepted_locus_id"]).reset_index(drop=True)


def block_summary(members, fields):
    rows = []
    assigned = members.loc[members.assignment_status.eq("assigned")]
    for (family, block_id), group in assigned.groupby(["family", "block_id"], sort=True):
        species = sorted(set(group.species))
        rows.append({"family": family, "block_id": block_id,
                     "block_type": "singleton" if len(group) == 1 else "synteny_block",
                     "loci": len(group), "species_count": len(species),
                     "multi_species": "yes" if len(species) >= 2 else "no",
                     "species": ";".join(species),
                     "accepted_locus_ids": ";".join(sorted(group.accepted_locus_id)),
                     "old_packaged_relation": ""})
    return pd.DataFrame(rows, columns=fields)


MODELLESS_FIELDS = {"genome", "family", "accepted_id", "status", "frame",
                    "order_scaffold", "order_start", "order_end"}
SCANNABLE_FRAMES = {"native", "agp_liftover", "exact_match"}


def modelless_loci(path, identity, chd):
    """Accepted copies without an ordered gene model, in the gene-order frame.

    They never enter block definition. Every accepted copy without an order
    identity must have exactly one disposition row; no other row is allowed.
    """
    frame = read(path)
    if not MODELLESS_FIELDS.issubset(frame):
        raise ValueError("model-less locus schema changed")
    no_order = identity.loc[identity.order_id.eq("")]
    expected = dict(zip(zip(no_order.genome, no_order.accepted_id), no_order.family))
    keys = list(zip(frame.genome, frame.accepted_id))
    if len(set(keys)) != len(keys) or set(keys) != set(expected):
        raise ValueError("model-less dispositions differ from accepted copies without order")
    scaffolds = defaultdict(set)
    for genome, chrom in zip(chd.genome, chd.chrom):
        scaffolds[genome].add(chrom)
    loci = []
    for row in frame.itertuples(index=False):
        if expected[(row.genome, row.accepted_id)] != row.family:
            raise ValueError(f"model-less family conflict: {row.accepted_id}")
        if row.status != "resolved" or row.frame not in SCANNABLE_FRAMES:
            if row.order_scaffold:
                raise ValueError(f"unscannable model-less locus has order coordinates: {row.accepted_id}")
            continue
        start, end = int(row.order_start), int(row.order_end)
        if row.order_scaffold not in scaffolds[row.genome] or not 0 < start <= end:
            raise ValueError(f"model-less locus outside gene-order frame: {row.accepted_id}")
        loci.append({"genome": row.genome, "family": BETA if row.family == "B-lectin" else row.family,
                     "accepted_id": row.accepted_id, "scaffold": row.order_scaffold,
                     "start": start, "end": end})
    return loci


def scan_modelless(matrix, chd, loci):
    """Count accepted model-less loci inside detected block search regions.

    Applies after block detection, to every family alike: a locus whose
    coordinates lie wholly within the bp span of the genes in the detected
    region counts as a target hit (PRESENT), the same as an annotated hit.
    Cells without a detected region, and block-member cells, are not changed.
    """
    spans = defaultdict(list)
    for row in chd.itertuples(index=False):
        spans[(row.genome, row.chrom)].append((int(row.position), int(row.start_bp), int(row.end_bp)))
    by_key = defaultdict(list)
    for locus in loci:
        by_key[(locus["genome"], locus["scaffold"], locus["family"])].append(locus)
    matrix = matrix.copy()
    audit = []
    for index, row in matrix.iterrows():
        if row.recovered != "yes" or row.search_route == "accepted_block_member" or row.region_start_rank == "":
            continue
        candidates = by_key.get((row.genome, row.scaffold, row.family), [])
        if not candidates:
            continue
        lo_rank, hi_rank = int(row.region_start_rank), int(row.region_end_rank)
        genes = [(s, e) for rank, s, e in spans[(row.genome, row.scaffold)] if lo_rank <= rank <= hi_rank]
        if not genes:
            raise ValueError(f"detected region has no ordered genes: {row.block_id} {row.genome}")
        lo, hi = min(s for s, _ in genes), max(e for _, e in genes)
        inside = []
        for locus in candidates:
            if lo <= locus["start"] and locus["end"] <= hi:
                where, gap = "inside", 0
                inside.append(locus["accepted_id"])
            elif locus["end"] < lo or locus["start"] > hi:
                where, gap = "outside", min(abs(lo - locus["end"]), abs(locus["start"] - hi))
            else:
                where, gap = "boundary_overlap", 0
            audit.append({"family": row.family, "block_id": row.block_id, "genome": row.genome,
                          "accepted_id": locus["accepted_id"], "scaffold": row.scaffold,
                          "locus_start": locus["start"], "locus_end": locus["end"],
                          "region_start_bp": lo, "region_end_bp": hi, "position": where,
                          "distance_bp": gap, "prior_status": row.status})
        if inside:
            prior = [x for x in str(row.target_unit_ids).split(";") if x]
            hits = sorted(set(prior) | set(inside))
            matrix.at[index, "status"] = "PRESENT"
            matrix.at[index, "n_target_units"] = len(hits)
            matrix.at[index, "target_unit_ids"] = ";".join(hits)
    fields = ["family", "block_id", "genome", "accepted_id", "scaffold", "locus_start", "locus_end",
              "region_start_bp", "region_end_bp", "position", "distance_bp", "prior_status"]
    return matrix, pd.DataFrame(audit, columns=fields)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--project", type=Path, required=True)
    p.add_argument("--sources", type=Path, required=True)
    p.add_argument("--identity", type=Path, required=True)
    p.add_argument("--backend", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    sources = pinned_sources(args.project, args.sources)
    recalc_path = Path(__file__).with_name("accepted") / "recalculate_synteny_paragraph.py"
    if hashlib.sha256(recalc_path.read_bytes()).hexdigest() != RECALC_SHA:
        raise ValueError("accepted recovery source hash changed")
    recalc = module(recalc_path, "existing_recovery")
    recon = module(sources["recon"], "existing_summary")
    identity = read(args.identity)
    authority = read(sources["placement_authority"])
    loci = read(args.backend / "pcwde_loci.tsv")
    blocks = read(args.backend / "synteny_blocks.tsv")
    panel = read(sources["historical_genome_panel"])
    if len(panel.loc[panel.status.eq("active")]) != 49:
        raise ValueError("final comparison panel changed")
    panel = panel.loc[panel.status.eq("active")].copy()
    fields = list(read(sources["historical_fields"]))
    block_fields = list(read(sources["historical_block_summary"]))
    meta = metadata(identity, authority, loci, fields, panel)
    members = membership(recalc, blocks, loci, meta, fields)
    if set(members.loc[members.assignment_status.eq("assigned"), "accepted_locus_id"]) != set(
            identity.loc[identity.order_id.ne(""), "accepted_id"]):
        raise ValueError("membership does not equal existing-order eligibility")
    args.out.mkdir(parents=True, exist_ok=False)
    build = args.out / "build"
    (build / "backend").mkdir(parents=True)
    shutil.copyfile(args.backend / "pcwde_loci.tsv", build / "backend/pcwde_loci.tsv")
    shutil.copyfile(args.backend / "synteny_blocks.tsv", build / "backend/synteny_blocks.tsv")
    members.to_csv(build / "table_s06_synteny_blocks_long.tsv", sep="\t", index=False)
    block_summary(members, block_fields).to_csv(
        build / "table_s06_synteny_block_summary.tsv", sep="\t", index=False)
    work = args.out / "work"
    work.mkdir()
    protein_chr = work / "protein_to_chromosome.tsv"
    protein_chr.write_text("protein_id\tscaffold\tbk_chr\tgenome\n")
    sidecar = work / "ordinary_family_new_loci_synteny_sidecar.tsv"
    sidecar.write_text("species_code\tbk_chromosome\tcontig\n")
    native_genome = {protein: row.genome for row in loci.itertuples(index=False)
                     for protein in row.protein_ids.split(",") if protein}
    display = dict(zip(members.loc[members.assignment_status.eq("assigned"), "builder_target_id"],
                       members.loc[members.assignment_status.eq("assigned"), "accepted_locus_id"]))
    if set(display) != set(native_genome):
        raise ValueError("recovery target display/genome bindings differ")
    recalc.CURRENT_BUILD, recalc.WORK = build, work
    recalc.PROTEIN_CHR, recalc.NEW_SIDECAR = protein_chr, sidecar
    recalc.TARGET_UNIT_DISPLAY_MAP = display
    recalc.TARGET_UNIT_GENOME_MAP = native_genome
    recalc.DISPLAY_BETA = BETA
    chd = read(sources["final_gene_order"])
    chd = chd.loc[chd.genome.isin(set(panel.accession))].copy()
    targets = members.loc[members.assignment_status.eq("assigned"), ["builder_target_id", "family"]].copy()
    targets.columns = ["protein_id", "family"]
    targets["family"] = targets.family.replace({BETA: "CBM13"})
    surface, tandem, sensitivity = recalc.current_tandem_and_proximity()
    matrix, _ = recalc.complete_recovery(chd, targets, panel, surface)
    if len(matrix) != matrix[["family", "block_id"]].drop_duplicates().shape[0] * 49:
        raise ValueError("recovery matrix is not blocks x 49 genomes")
    if "modelless_loci" not in sources:
        raise ValueError("model-less locus disposition is not pinned")
    matrix, scan = scan_modelless(matrix, chd, modelless_loci(sources["modelless_loci"], identity, chd))
    scan.to_csv(args.out / "modelless_scan_audit.tsv", sep="\t", index=False)
    matrix.to_csv(args.out / "current_region_status_complete.tsv", sep="\t", index=False)
    stats = pd.DataFrame(recon.summarize("existing_order", build, matrix, tandem, pd.DataFrame(sensitivity)))
    stats.to_csv(args.out / "synteny_statistics.tsv", sep="\t", index=False)
    print(f"EXISTING_RECOVERY_PASS accepted={len(members)} eligible={len(targets)} "
          f"blocks={len(matrix) // 49} cells={len(matrix)}")


if __name__ == "__main__":
    main()
