#!/usr/bin/env python3
"""Resolve accepted copies to existing ordered models, without creating order rows.

One family-blind criterion: identity is an exact (genome, protein ID) match,
or a unique amino-acid match, to a model with a position in either pinned order
route (final gene order, or a supplemental-anchored reconstructed-scaffold
order row). A raw locus, placement, or overlapping annotation is never an
identity edge, and no position is created here.
"""

import sys as _source_sys
from pathlib import Path as _SourcePath
_source_repo = next(p for p in _SourcePath(__file__).resolve().parents if (p / 'config/logical_inputs.tsv').is_file())
_source_sys.path.insert(0, str(_source_repo / 'scripts/release'))
from source_aliases import source_path
import argparse
import csv
import hashlib
from collections import defaultdict
from pathlib import Path


def rows(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1048576), b""):
            h.update(chunk)
    return h.hexdigest()


def fasta(path):
    result = {}
    name = None
    with Path(path).open() as handle:
        for line in handle:
            if line.startswith(">"):
                name = line[1:].split()[0]
                if name in result:
                    raise ValueError(f"duplicate FASTA ID {name}: {path}")
                result[name] = []
            elif name is not None:
                result[name].append(line.strip().upper())
    return {key: "".join(value).rstrip("*") for key, value in result.items()}


def source_map(project, sources, genome_map):
    src = {r["source_file"]: r for r in rows(sources)}
    mapping = {r["genome"]: r["source_file"] for r in rows(genome_map)}
    if len(mapping) != len(rows(genome_map)):
        raise ValueError("duplicate genome source map")
    return {genome: (source_path(project, src[key]["proteome_project_path"]),
                     src[key]["proteome_sha256"],
                     source_path(project, src[key]["gff_project_path"]),
                     src[key]["gff_sha256"])
            for genome, key in mapping.items()}


def position_index(original_order, reconstructed_order, supplemental):
    """Existing ordered positions from both pinned order routes, family-blind.

    Positions are those of the native backend's effective gene index: the frozen
    final gene order, with each reconstructed (genome, scaffold) replaced by its
    reconstructed order. An identity is admitted by route 1 when it is in the
    final order, or by route 2 when the pinned supplemental-interval input names
    the same (genome, scaffold, protein) on a reconstructed scaffold, the two
    admission paths the backend consumes. A reconstructed row without that
    record, or a supplemental interval without a reconstructed row, is not a
    position.
    """
    recon = defaultdict(dict)
    for row in rows(reconstructed_order):
        recon[(row["genome"], row["scaffold"])][row["protein_id"]] = row["position"]
    index = defaultdict(set)
    with Path(original_order).open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = (row["genome"], row["protein_id"])
            frame = (row["genome"], row["chrom"])
            if frame in recon:
                if row["protein_id"] in recon[frame]:
                    index[key].add(("final_order", row["chrom"], recon[frame][row["protein_id"]]))
                else:
                    index.setdefault(key, set())  # superseded by its scaffold's reconstructed order
            else:
                index[key].add(("final_order", row["chrom"], row["position"]))
    anchors = {(r["genome"], r["scaffold"], r["protein_id"]) for r in rows(supplemental)}
    if len(anchors) != len(rows(supplemental)):
        raise ValueError("duplicate supplemental protein/scaffold segment")
    final_ids = set(index)
    for genome, scaffold, protein in sorted(anchors):
        position = recon.get((genome, scaffold), {}).get(protein)
        if position is None:
            continue
        if (genome, protein) in final_ids:
            raise ValueError(f"identity admitted by both order routes: {(genome, protein)}")
        index[(genome, protein)].add(("reconstructed_order", scaffold, position))
    return {k: v for k, v in index.items() if v}


def alias_frame(positions, genome, protein, gff_path, gff_sha, curated_scaffold):
    if digest(gff_path) != gff_sha:
        raise ValueError(f"annotation pin mismatch: {gff_path}")
    gff_frames = set()
    with gff_path.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] not in {"mRNA", "transcript"}:
                continue
            attrs = dict(part.split("=", 1) for part in fields[8].split(";") if "=" in part)
            if attrs.get("ID") == protein:
                gff_frames.add(fields[0])
    order_frames = {frame for _, frame, _ in positions}
    if len(gff_frames) != 1 or len(order_frames) != 1 or gff_frames != order_frames:
        raise ValueError(f"alias annotation/order frame ambiguity: {genome}/{protein}")
    if curated_scaffold and curated_scaffold not in gff_frames:
        raise ValueError(f"alias accepted/annotation frame mismatch: {genome}/{protein}")


def check_member_order(members, member_order, index):
    """The member-layer order is a second view of the same positions: where it
    lists an accepted identity, it must agree with the shared position index."""
    listed = defaultdict(set)
    for row in rows(member_order):
        listed[(row["genome"], row["protein_id"])].add((row["scaffold"], row["position"]))
    accepted = {(r["genome"], r["protein_id"]) for r in members}
    if set(listed) - accepted:
        raise ValueError(f"member order lists non-accepted identities: {sorted(set(listed) - accepted)[:5]}")
    for key, positions in listed.items():
        derived = {(frame, pos) for _, frame, pos in index.get(key, set())}
        if positions != derived:
            raise ValueError(f"member order disagrees with shared position index: {key}: {positions} vs {derived}")
    return len(listed)


def resolve(members, original_order, member_order, accepted, sources, reconstructed_order, supplemental):
    accepted_ids = {(r["genome"], r["protein_id"]) for r in members}
    if len(accepted_ids) != len(members):
        raise ValueError("duplicate accepted member")
    by_protein = defaultdict(set)
    for r in members:
        by_protein[r["protein_id"]].add(r["family"])
    if any(len(v) != 1 for v in by_protein.values()):
        raise ValueError("cross-family accepted identity collision")
    index = position_index(original_order, reconstructed_order, supplemental)
    check_member_order(members, member_order, index)
    results = []
    cache = {}
    for r in members:
        key = (r["genome"], r["protein_id"])
        positions = index.get(key, set())
        if positions:
            chosen, method, issue = r["protein_id"], "exact_order_id", ""
        else:
            seq = accepted.get(r["unit_id"], "") if r["unit_id"] else ""
            chosen, method = "", ""
            if not seq:
                issue = "accepted_sequence_missing" if r["unit_id"] else "no_existing_order_position"
            elif r["genome"] not in sources:
                issue = "pinned_proteome_missing"
            else:
                if r["genome"] not in cache:
                    path, expected, _, _ = sources[r["genome"]]
                    if digest(path) != expected:
                        raise ValueError(f"proteome pin mismatch: {path}")
                    cache[r["genome"]] = fasta(path)
                hits = sorted(pid for pid, aa in cache[r["genome"]].items()
                              if aa == seq and index.get((r["genome"], pid)))
                if len(hits) > 1:
                    raise ValueError(f"ambiguous exact-AA ordered models: {key}: {hits}")
                chosen = hits[0] if hits else ""
                method = "exact_aa_existing_order" if chosen else ""
                issue = "" if chosen else "no_exact_aa_existing_order_model"
                if chosen:
                    _, _, gff_path, gff_sha = sources[r["genome"]]
                    alias_frame(index[(r["genome"], chosen)], r["genome"], chosen, gff_path, gff_sha,
                                r["curated_scaffold"])
        routes = sorted({p[0] for p in index.get((r["genome"], chosen), set())}) if chosen else []
        results.append({"genome": r["genome"], "family": r["family"],
                        "accepted_id": r["protein_id"], "unit_id": r["unit_id"],
                        "order_id": chosen, "method": method, "issue": issue,
                        "order_segments": len(index[(r["genome"], chosen)]) if chosen else 0,
                        "order_route": ";".join(routes)})
    chosen_ids = [(x["genome"], x["order_id"]) for x in results if x["order_id"]]
    if len(chosen_ids) != len(set(chosen_ids)):
        raise ValueError("two accepted copies resolve to one ordered model")
    return results


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--project", type=Path, required=True)
    p.add_argument("--members", type=Path, required=True)
    p.add_argument("--order", type=Path, required=True)
    p.add_argument("--member-order", type=Path, required=True)
    p.add_argument("--accepted-fasta", type=Path, required=True)
    p.add_argument("--accepted-sha256", required=True)
    p.add_argument("--sources", type=Path, required=True)
    p.add_argument("--genome-map", type=Path, required=True)
    p.add_argument("--order-sources", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    if digest(args.accepted_fasta) != args.accepted_sha256:
        raise ValueError("accepted sequence pin mismatch")
    pinned = {r["role"]: r for r in rows(args.order_sources)}
    route2 = {}
    for role in ("reconstructed_order", "missing_order"):
        path = source_path(args.project, pinned[role]["project_relative_path"])
        if digest(path) != pinned[role]["sha256"]:
            raise ValueError(f"order-route pin mismatch: {role}")
        route2[role] = path
    result = resolve(rows(args.members), args.order, args.member_order, fasta(args.accepted_fasta),
                     source_map(args.project, args.sources, args.genome_map),
                     route2["reconstructed_order"], route2["missing_order"])
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(result[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(result)
    from collections import Counter
    print(f"accepted={len(result)} disposition={dict(Counter('eligible' if r['order_id'] else r['issue'] for r in result))}")


if __name__ == "__main__":
    main()
