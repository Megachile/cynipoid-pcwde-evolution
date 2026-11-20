#!/usr/bin/env python3
"""Join the four frozen PCWDE sequence sets without assuming locus equivalence.

Identity and pinned crosswalks take precedence. Unlinked source models retain
their own rows. Sequence identity or species/ordinal similarity is never a link.
"""
import argparse
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path

FAMILIES = ("GH5_2", "PL1_8", "PL4_4")
SETS = ("locus_tree", "selection", "catalytic_alignment", "origin_tree")
# Public inquiline labels mirror the accepted Table S7 species alias authority.
PUBLIC_TIP_ALIASES = {
    "Ceroptres_(DI_inquiline)": "Ceroptres sp. ex Druon ignotum",
    "Synergus_(BK7_inquiline)": "Synergus sp. ex Belonocnema kinseyi",
    "Synergus_(DI_inquiline)": "Synergus sp. ex Druon ignotum",
    "Synergus_sp._(D._quercuslanigerum_inquiline)": "Synergus sp. ex Druon quercuslanigerum",
}


def unique_reciprocal_overlaps(intervals):
    """Return unique reciprocal CDS overlaps on one contig and strand.

    Input records contain family, source_id, side (older/current), contig,
    strand, start and end. This tier is usable only with pinned coordinates.
    """
    candidates = defaultdict(set)
    for i, left in enumerate(intervals):
        for right in intervals[i + 1:]:
            if left["side"] == right["side"] or any(left[k] != right[k] for k in ("family", "contig", "strand")):
                continue
            if max(int(left["start"]), int(right["start"])) <= min(int(left["end"]), int(right["end"])):
                candidates[left["source_id"]].add(right["source_id"])
                candidates[right["source_id"]].add(left["source_id"])
    return {source: next(iter(hits)) for source, hits in candidates.items()
            if len(hits) == 1 and candidates[next(iter(hits))] == {source}}


def read(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({field: row.get(field, "") for field in fields} for row in rows)


class Resolver:
    def __init__(self):
        self.rows = {}
        self.aliases = defaultdict(set)
        self.source_owner = {}
        self.pending = []

    def add_row(self, family, canonical, **meta):
        key = (family, canonical)
        if key not in self.rows:
            self.rows[key] = dict(family=family, canonical_id=canonical,
                                  public_species_name="", source="", source_accession="",
                                  bk_chromosome="", link_status="not_linked", link_note="",
                                  exclusion_note="", **{f"in_{s}": "false" for s in SETS},
                                  **{f"source_ids_{s}": "" for s in SETS})
        self.rows[key].update({k: v for k, v in meta.items() if v})
        self.aliases[key].add(canonical)
        return key

    def alias(self, key, value):
        if value:
            self.aliases[key].add(value)

    def resolve(self, family, source_id):
        matches = [key for key, ids in self.aliases.items()
                   if key[0] == family and source_id in ids]
        if len(matches) > 1:
            raise ValueError(f"ambiguous pinned ID {family} {source_id}: {matches}")
        return matches[0] if matches else None

    def candidates(self, family, source_id):
        return {key for key, ids in self.aliases.items()
                if key[0] == family and source_id in ids}

    def collect(self, family, source_id, set_name, meta=None):
        self.pending.append((family, source_id, set_name, meta or {}))

    def assign_collected(self):
        """Assign a complete set of source records without first-source priority."""
        records = sorted(self.pending, key=lambda x: (x[0], x[2], x[1]))
        source_sets = defaultdict(set)
        for family, source_id, set_name, _ in records:
            source_sets[family, source_id].add(set_name)
        for (family, source_id), member_sets in sorted(source_sets.items()):
            if not self.candidates(family, source_id):
                self.add_row(family, source_id,
                             link_status="source_exact" if len(member_sets) > 1 else "not_linked")
        keys = [(f, s, n) for f, s, n, _ in records]
        if len(keys) != len(set(keys)):
            raise ValueError("duplicate producer source ID")
        candidates = {}
        for family, source_id, set_name, _ in records:
            candidates[(family, set_name, source_id)] = self.candidates(family, source_id)
        # A row cannot hold two distinct IDs from the same producer. Reject all
        # participants, including the first source in lexical/input order.
        users = defaultdict(set)
        for record, matches in candidates.items():
            for match in matches:
                users[(record[0], record[1], match)].add(record)
        collisions = {record for members in users.values() if len(members) > 1
                      for record in members}
        for family, source_id, set_name, meta in records:
            record = (family, set_name, source_id)
            matches = candidates[record]
            if len(matches) == 1 and record not in collisions:
                key = next(iter(matches))
            else:
                canonical = f"{set_name}:{source_id}"
                key = self.add_row(family, canonical, link_status="not_linked")
                why = "ambiguous pinned alias" if len(matches) > 1 else "many-to-one source collision" if record in collisions else meta.get("link_note", "No unique pinned sequence link")
                self.rows[key]["link_note"] = why + ("; candidates: " + ";".join(sorted(k[1] for k in matches)) if matches else "")
                if len(matches) == 1 and record in collisions:
                    candidate = self.rows[next(iter(matches))]
                    for field in ("public_species_name", "source", "source_accession", "bk_chromosome"):
                        self.rows[key][field] = candidate[field]
            row = self.rows[key]
            for field, value in meta.items():
                if field in row and field != "link_note" and not row[field] and value:
                    row[field] = value
            row[f"source_ids_{set_name}"] = source_id
            row[f"in_{set_name}"] = "true"
            self.source_owner[record] = key
        self.pending.clear()

    def mark(self, family, source_id, set_name, meta=None):
        source_key = (family, set_name, source_id)
        if source_key in self.source_owner:
            raise ValueError(f"duplicate source ID {source_key}")
        key = self.source_owner.get((family, "catalytic_alignment", source_id)) if set_name == "origin_tree" else None
        key = key or self.resolve(family, source_id)
        if key is not None and self.rows[key][f"in_{set_name}"] == "true":
            # Two distinct producer records cannot occupy one sequence row.
            # The shared source model is therefore not a unique sequence link.
            key = None
        if key is None:
            meta = meta or {}
            canonical = source_id
            if (family, canonical) in self.rows:
                canonical = f"{set_name}:{source_id}"
            key = self.add_row(family, canonical, link_note=meta.get("link_note", "No unique pinned sequence link"))
            self.alias(key, source_id)
            self.rows[key].update({k: v for k, v in meta.items() if k in self.rows[key] and v})
        row = self.rows[key]
        for field, value in (meta or {}).items():
            if field in row and field != "link_note" and not row[field] and value:
                row[field] = value
        field = f"source_ids_{set_name}"
        ids = set(filter(None, row[field].split(";")))
        ids.add(source_id)
        row[field] = ";".join(sorted(ids))
        row[f"in_{set_name}"] = "true"
        self.source_owner[source_key] = key


def build(args):
    root, build = args.data_root / "repository_inputs", args.build
    resolver = Resolver()
    for copy in read(build / "members/copy_placements.tsv"):
        if copy["family"] not in FAMILIES:
            continue
        proteins = [s.strip() for s in copy["protein_ids"].split(",") if s.strip()]
        for protein in proteins:
            key = resolver.add_row(copy["family"], protein,
                                   public_species_name=copy["species"], source="genome model",
                                   bk_chromosome=copy["bk_chromosome"], link_status="current_copy")
            if len(proteins) == 1:
                resolver.alias(key, copy["locus_protein_id"])
    locus = read(build / "members/locus_tips.tsv")
    repeated_locus_proteins = Counter((t["family"], t["source_protein_id"]) for t in locus)
    for tip in locus:
        family, source_id = tip["family"], tip["source_protein_id"]
        if repeated_locus_proteins[family, source_id] > 1:
            current = resolver.resolve(family, source_id)
            if current:
                resolver.rows[current]["link_status"] = "not_linked"
                resolver.rows[current]["link_note"] = "Shared source protein for multiple tree tips; sequence equivalence not established"
        key = None if repeated_locus_proteins[family, source_id] > 1 else resolver.resolve(family, source_id)
        if key is None:
            key = resolver.add_row(family, tip["display_tree_tip"] if repeated_locus_proteins[family, source_id] > 1 else source_id, public_species_name=tip["species"],
                                   source="Trinity transcript" if tip["source_kind"] == "transcript_assembly" else "genome model",
                                   link_status="not_linked" if repeated_locus_proteins[family, source_id] > 1 else "locus_tip",
                                   link_note="Shared source protein for multiple tree tips; sequence equivalence not established" if repeated_locus_proteins[family, source_id] > 1 else "Transcript tip outside current copy table" if tip["source_kind"] == "transcript_assembly" else "")
        if repeated_locus_proteins[family, source_id] == 1:
            resolver.alias(key, source_id)
        if repeated_locus_proteins[family, source_id] == 1:
            resolver.alias(key, tip["raw_tip"])
        resolver.alias(key, tip["display_tree_tip"])
        resolver.collect(family, tip["display_tree_tip"], "locus_tree")
    # Crosswalk the complete union, including rows created only by locus tips.
    # A crosswalk target can itself be ambiguous; attach it to every candidate
    # and let the symmetric assignment gate keep the source unlinked.
    for item in read(root / "tables/authority/table_s2/ordinary_identifier_crosswalk.tsv"):
        family = item["family"]
        if family not in FAMILIES:
            continue
        targets = resolver.candidates(family, item["stable_protein_id"])
        if not targets and "TRINITY" in item["stable_protein_id"]:
            # Locus tips namespace the pinned Trinity protein with its source
            # code (for example AF|AF_TRINITY_...), while the older selection
            # lock uses the unprefixed stable protein ID.
            namespaced = item["source_species_code"] + "|" + item["stable_protein_id"]
            targets = resolver.candidates(family, namespaced)
        for key in targets:
            resolver.alias(key, item["alias_id"])
            resolver.alias(key, item["stable_transcript_id"])
    selection = read(build / "table_s7/selection_membership_and_public_source.tsv")
    for item in selection:
        family, source_id = item["gene_family"], item["protein_id"]
        if family not in FAMILIES:
            continue
        resolver.collect(family, source_id, "selection", {
            "public_species_name": item["species"], "source": item["source"],
            "source_accession": item["source_accession"],
            "link_note": f"Selection source {item['source_accession']}; no pinned link to current copy or locus"})
    catalytic = {(r["family"], r["sequence_id"])
                 for r in read(build / "table_s7/catalytic_residue_by_sequence.tsv")
                 if r["alignment_set"] == "de_novo" and r["site_role"] == "catalytic"}
    for family, source_id in sorted(catalytic):
        resolver.collect(family, source_id, "catalytic_alignment")
    origin = read(root / "deposit/base/sidecars/origin_context_full_tip_mapping.tsv")
    public_origin_labels = {}
    for family in FAMILIES:
        for item in read(root / f"deposit/base/sidecars/origin_context_{family}_tip_labels.tsv"):
            public_origin_labels[(family, item["source_sequence_id"])] = item["base_display_label"]
    outgroups = Counter()
    for item in origin:
        family, source_id = item["family"], item["alignment_sequence_id"]
        if family not in FAMILIES:
            continue
        known = any(f == family and s == source_id for f, s, _, _ in resolver.pending)
        if not known and not resolver.candidates(family, source_id):
            if source_id.startswith("cynipoid|") or source_id.split("|", 1)[0].endswith("_ENA"):
                raise ValueError(f"unmapped origin cynipoid ID {family} {source_id}")
            outgroups[family] += 1
            continue
        resolver.collect(family, source_id, "origin_tree", {
            "public_species_name": public_origin_labels.get((family, source_id), "")})
    resolver.assign_collected()
    counts = Counter((family, name) for family, name, _ in resolver.source_owner)
    catalytic_summary = read(build / "catalytic_current/catalytic_site_conservation_summary.tsv")
    catalytic_totals = {}
    for family in FAMILIES:
        values = {int(r["n_total_sequences"]) for r in catalytic_summary
                  if r["family"] == family and r["alignment_set"] == "de_novo"
                  and r["site_role"] == "catalytic"}
        if len(values) != 1:
            raise ValueError(f"ambiguous current catalytic total for {family}: {values}")
        catalytic_totals[family] = values.pop()
    expected = {"locus_tree": (255, 294, 358), "selection": (257, 300, 354),
                "catalytic_alignment": tuple(catalytic_totals[f] for f in FAMILIES)}
    for name, values in expected.items():
        for family, value in zip(FAMILIES, values):
            if counts[family, name] != value:
                raise ValueError(f"producer count mismatch {family} {name}: {counts[family,name]} != {value}")
    if tuple(counts[f, "origin_tree"] + outgroups[f] for f in FAMILIES) != (427, 364, 444):
        raise ValueError(f"origin accounting failed: {counts}, {outgroups}")
    removed = read(root / "deposit/base/sidecars/removed_cb_origin_bp_tips.tsv")
    removed_origin = Counter()
    for item in removed:
        for family in FAMILIES:
            if item["tree"] == f"trees/origin_context_{family}_full.treefile":
                removed_origin[family] += 1
    if tuple(removed_origin[f] for f in FAMILIES) != (0, 1, 2):
        raise ValueError(f"origin removal record changed: {removed_origin}")
    rows = [resolver.rows[key] for key in sorted(resolver.rows)]
    rows = [row for row in rows if any(row[f"in_{name}"] == "true" for name in SETS)]
    for item in removed:
        for family in FAMILIES:
            if item["tree"] != f"trees/origin_context_{family}_full.treefile":
                continue
            matched = [r for r in rows if r["family"] == family and
                       (r["canonical_id"] == item["tip"] or
                        item["tip"] in r["source_ids_catalytic_alignment"].split(";"))]
            if not matched and (family, item["tip"]) not in catalytic:
                # C's current catalytic producer has already removed this tip.
                # The removal sidecar retains its historical inference record.
                continue
            if len(matched) != 1:
                raise ValueError(f"removed origin tip not uniquely accounted for: {family} {item['tip']}")
            matched[0]["exclusion_note"] = "Removed from deposited origin tree after inference under the recorded Bassettia/Callirhytis rule"
    tip_species = defaultdict(set)
    for row in rows:
        if row["public_species_name"]:
            for label in row["source_ids_locus_tree"].split(";"):
                match = re.fullmatch(r"(cynipoid\|.+)_\d+", label)
                if match:
                    tip_species[row["family"], match.group(1)].add(row["public_species_name"])
    for row in rows:
        if row["public_species_name"]:
            continue
        match = re.fullmatch(r"(cynipoid\|.+)_\d+", row["canonical_id"])
        if match and len(tip_species[row["family"], match.group(1)]) == 1:
            row["public_species_name"] = next(iter(tip_species[row["family"], match.group(1)]))
        elif match:
            row["public_species_name"] = PUBLIC_TIP_ALIASES.get(match.group(1).split("|", 1)[1], "")
    unnamed = [(row["family"], row["canonical_id"]) for row in rows if not row["public_species_name"]]
    if unnamed:
        raise ValueError(f"blank public_species_name in {len(unnamed)} roster rows: {unnamed[:12]}")
    for family in FAMILIES:
        for name in SETS:
            n = sum(r["family"] == family and r[f"in_{name}"] == "true" for r in rows)
            if n != counts[family, name]:
                raise ValueError(f"many-to-one mapping within {family} {name}: {n} rows for {counts[family,name]} IDs")
    for row in rows:
        flags = [row[f"in_{name}"] for name in SETS]
        if len(set(flags)) > 1 and row["link_status"] == "not_linked":
            row["link_note"] = row["link_note"] or "Model version difference; locus equivalence not established"
    fields = list(rows[0])
    write(args.out / "pcwde_sequence_set_membership.tsv", rows, fields)
    summary = []
    for family in FAMILIES:
        family_rows = [r for r in rows if r["family"] == family]
        for name in SETS:
            selected = [r for r in family_rows if r[f"in_{name}"] == "true"]
            summary.append(dict(family=family, metric=name, total=len(selected),
                                linked_by_pinned_id=sum(r["link_status"] != "not_linked" for r in selected),
                                linked_by_genomic_overlap=0,
                                not_linked=sum(r["link_status"] == "not_linked" for r in selected),
                                outgroup_count=outgroups[family] if name == "origin_tree" else 0))
        summary.append(dict(family=family, metric="origin_inference_before_recorded_pruning",
                            total=counts[family, "origin_tree"] + outgroups[family] + removed_origin[family],
                            linked_by_pinned_id="", linked_by_genomic_overlap="",
                            not_linked="", outgroup_count=outgroups[family]))
        summary.append(dict(family=family, metric="origin_removed_after_inference",
                            total=removed_origin[family], linked_by_pinned_id="",
                            linked_by_genomic_overlap="", not_linked="", outgroup_count=""))
        for i, left in enumerate(SETS):
            for right in SETS[i + 1:]:
                linked = [r for r in family_rows if r["link_status"] != "not_linked"]
                summary.append(dict(family=family, metric=f"{left}_vs_{right}_membership_differs",
                                    total=sum(r[f"in_{left}"] != r[f"in_{right}"] for r in linked),
                                    linked_by_pinned_id="", linked_by_genomic_overlap="",
                                    not_linked="", outgroup_count=""))
    write(args.out / "pcwde_sequence_set_summary.tsv", summary,
          ["family", "metric", "total", "linked_by_pinned_id", "linked_by_genomic_overlap", "not_linked", "outgroup_count"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--build", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    build(parser.parse_args())
