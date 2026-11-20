"""Supplement prose facts computed from the same intermediates the tables and figures use.

build_numbers.py adds the returned dictionary to numbers.json under "supplement";
the supplement template (authority/supplement/supplement_source.md) reads these
keys as {{supplement.<key>}} tokens. Every value is counted here from a build or
deposit file; nothing is typed in.

Environment: project Python 3 (standard library only).
"""
import csv
import gzip
import json
import math
import re
import struct
from collections import Counter
from parent_proteins import origin_parent

# BioProjects of the RNA-seq data first presented in this study (Data availability
# statement). Transcriptome assemblies from any other source are previously published.
STUDY_RNASEQ_BIOPROJECTS = ("PRJNA1502425", "PRJNA1504442")


def _rows(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _one(values, what):
    values = set(values)
    if len(values) != 1:
        raise SystemExit(f"supplement numbers: expected one value for {what}, found {sorted(values)}")
    return values.pop()


def check_like_unit(left_unit, left_count, right_unit, right_count, label):
    """Compare counts only when their declared biological units match."""
    if left_unit != right_unit:
        raise SystemExit(f"supplement numbers: unlike units for {label}: {left_unit} vs {right_unit}")
    if left_count != right_count:
        raise SystemExit(f"supplement numbers: {label} {left_unit} count differs: {left_count} vs {right_count}")


def panel_facts(source_panel):
    """Genome panel accounting (Table S1 source panel)."""
    rows = _rows(source_panel)
    genomes = [r for r in rows if r["source_type"] == "genome"]
    retained = [r for r in genomes if r["status"] == "retained"]
    excluded = [r for r in genomes if r["status"] == "excluded"]
    outgroup = [r for r in genomes if r["superfamily"] == "Platygastroidea"]

    def reason(r):
        text = r["exclusion_reason"]
        if text.startswith("Below 40% BUSCO completeness"):
            return "below_busco"
        if text.startswith("Redundant or replaced") or text.startswith("superseded"):
            return "redundant"
        if text.startswith("Other assembly-quality failure"):
            return "other_quality"
        raise SystemExit(f"supplement numbers: unclassified exclusion reason: {text!r}")

    reasons = Counter(reason(r) for r in excluded)
    floors = {float(m) for r in excluded for m in re.findall(r"below (\d+(?:\.\d+)?)%", r["exclusion_reason"])}
    source = Counter(r["source_and_assembly"].split(" / ")[0] for r in genomes)
    ibarrai = [r for r in genomes if r["species"] == "Cecinothofagus ibarrai"]
    epri = [r for r in retained if r["species"] == "Euceroptres primus"]
    transcriptomes = [r for r in rows if r["source_type"] == "transcriptome"]
    t_busco = [float(r["proteome_busco_C"]) for r in transcriptomes]
    family = Counter(r["family"] for r in retained if r["superfamily"] != "Platygastroidea")
    tribes = Counter(r["tribe"] for r in transcriptomes)
    if len(ibarrai) != 1 or len(epri) != 1 or len(outgroup) != 1:
        raise SystemExit("supplement numbers: panel lacks one C. ibarrai, retained E. primus or outgroup genome")
    return {
        "genomes_considered": len(genomes),
        "cynipoid_genomes_considered": len(genomes) - len(outgroup),
        "outgroup_genomes": len(outgroup),
        "ncbi_cynipoid_assemblies": sum(r["source_and_assembly"].startswith("NCBI / ")
                                        and r["superfamily"] != "Platygastroidea" for r in genomes),
        "ena_read_assemblies": source["ENA (assembled)"],
        "ena_read_origin_assemblies": sum(r["source_and_assembly"].startswith("ENA (assembled)") and r["accession_or_project"].startswith("PRJEB") for r in genomes),
        "ncbi_read_origin_assemblies": sum(r["source_and_assembly"].startswith("ENA (assembled)") and r["accession_or_project"].startswith("PRJNA") for r in genomes),
        "ncbi_genbank_later_assemblies": source["NCBI GenBank"],
        "reference_proteomes": sum(r["source_type"] == "reference_proteome" for r in rows),
        "busco_floor_pct": _one(floors, "BUSCO floor"),
        "excluded_below_busco": reasons["below_busco"],
        "excluded_redundant": reasons["redundant"],
        "excluded_other_quality": reasons["other_quality"],
        "excluded_after_busco": reasons["redundant"] + reasons["other_quality"],
        "retained_genomes": len(retained),
        "retained_diplolepididae": family["Diplolepididae"],
        "retained_figitidae": family["Figitidae"],
        "retained_cynipidae": family["Cynipidae"],
        "c_ibarrai_genome_busco_pct": float(ibarrai[0]["assembly_busco_C"]),
        "e_primus_genome_busco_pct": float(epri[0]["assembly_busco_C"]),
        "e_primus_proteins": int(epri[0]["sequence_count"]),
        "e_primus_assembly_mb": float(epri[0]["assembly_size_mb"]),
        "e_primus_n50_mb": int(epri[0]["n50_bp"]) / 1e6,
        "transcriptome_synergini": tribes["Synergini"],
        "transcriptome_ceroptresini": tribes["Ceroptresini"],
        "transcriptome_cynipini": tribes["Cynipini"],
        "transcriptome_assemblies": len(transcriptomes),
        "transcriptome_new": sum(r["accession_or_project"] in STUDY_RNASEQ_BIOPROJECTS for r in transcriptomes),
        "transcriptome_published": sum(r["accession_or_project"] not in STUDY_RNASEQ_BIOPROJECTS for r in transcriptomes),
        "retained_genome_proteins": sum(int(r["sequence_count"]) for r in retained),
        "transcriptome_protein_busco_min_pct": min(t_busco),
        "transcriptome_protein_busco_max_pct": max(t_busco),
    }


def bassettia_attribution_facts(removed_tips):
    """Identity and patristic distance of the CB-attributed B. pallida tips (removal record)."""
    rows = _rows(removed_tips)
    if not rows or not all(r["rule"].startswith("removed after inference") for r in rows):
        raise SystemExit("supplement numbers: unexpected B. pallida removal record")
    identity = [float(r["aligned_protein_identity_to_callirhytis"]) * 100 for r in rows]
    distance = [float(r["distance_to_nearest_callirhytis_subst_per_site"]) for r in rows]
    return {"bassettia_identity_min_pct": min(identity), "bassettia_identity_max_pct": max(identity),
            "bassettia_distance_min": min(distance), "bassettia_distance_max": max(distance)}


def e_primus_facts(facts_dir):
    """E. primus sequencing, assembly and annotation facts from the pinned source files."""
    pbi = gzip.open(facts_dir / "e_primus_hifi_reads.bam.pbi").read()
    magic, _version, _flags, n_reads = struct.unpack("<4sIHI", pbi[:14])
    if magic != b"PBI\x01":
        raise SystemExit("supplement numbers: E. primus HiFi index is not a PacBio BAM index")
    start = 32 + 4 * n_reads  # header, then int32 read-group ids, qStart and qEnd columns
    q_start = struct.unpack(f"<{n_reads}i", pbi[start:start + 4 * n_reads])
    q_end = struct.unpack(f"<{n_reads}i", pbi[start + 4 * n_reads:start + 8 * n_reads])
    hifi_bp = sum(e - s for s, e in zip(q_start, q_end))

    lengths = {}
    for line in (facts_dir / "e_primus_assembly.fasta.fai").read_text().splitlines():
        name, length = line.split("\t")[:2]
        lengths[name] = int(length)
    total = sum(lengths.values())
    chromosomes = [n for n in lengths if n.startswith("Chromosome")]
    contigs = [n for n in lengths if n.startswith("ctg")]
    if len(chromosomes) + len(contigs) + ("Mitogenome" in lengths) != len(lengths):
        raise SystemExit("supplement numbers: unexpected E. primus sequence names")
    ordered = sorted(lengths.values(), reverse=True)

    def nx(fraction):
        running = 0
        for index, length in enumerate(ordered, 1):
            running += length
            if running >= fraction * total:
                return length, index

    n90, l90 = nx(0.9)
    _n50, l50 = nx(0.5)
    contig_bp = sum(lengths[n] for n in contigs)

    busco = json.loads((facts_dir / "e_primus_genome_busco_hymenoptera_odb12.json").read_text())["results"]
    return {
        "e_primus_hifi_reads": n_reads,
        "e_primus_hifi_gb": hifi_bp / 1e9,
        "e_primus_sequences": len(lengths),
        "e_primus_chromosomes": len(chromosomes),
        "e_primus_l50": l50,
        "e_primus_n90_mb": n90 / 1e6,
        "e_primus_l90": l90,
        "e_primus_unplaced_contigs": len(contigs),
        "e_primus_unplaced_contig_kb": contig_bp / 1e3,
        "e_primus_unplaced_contig_pct": 100 * contig_bp / total,
        "e_primus_busco_single_pct": float(busco["Single copy percentage"]),
        "e_primus_busco_duplicated_pct": float(busco["Multi copy percentage"]),
        "e_primus_busco_fragmented_pct": float(busco["Fragmented percentage"]),
        "e_primus_busco_missing_pct": float(busco["Missing percentage"]),
        "e_primus_busco_n": int(busco["n_markers"]),
        "e_primus_busco_complete_pct": float(busco["Complete percentage"]),
    }


def later_proteome_facts(facts_dir):
    """Count the nine later-added genome proteomes (NCBI GenBank, July 2026) searched against the CynOG catalog."""
    rows = _rows(facts_dir / "current_proteomes.tsv")
    expected = {"AHIE", "AMIN", "CQUE", "DRUB", "PCAN", "PLAM", "Pe_aceris", "SHIS", "TPHE"}
    if {r["species_code"] for r in rows} != expected or len(rows) != len(expected):
        raise SystemExit("supplement numbers: later-added proteome roster differs from pinned record")
    return {"later_added_proteomes": len(rows)}


def cazy_lineage_group_facts(lineage_groups):
    """Genome counts per Table S7 lineage summary group, from the lineage table Table S7 reads."""
    counts = Counter(r["lineage_group"] for r in _rows(lineage_groups))
    names = {"Cynipini": "cynipini", "non-Cynipini Cynipidae": "other_cynipidae",
             "Diplolepididae": "diplolepididae", "E. primus": "e_primus"}
    if set(counts) - set(names) - {"outside lineage summaries"} or any(counts[n] == 0 for n in names):
        raise SystemExit("supplement numbers: unexpected Table S7 lineage groups")
    return {"lineage_group_genomes": {key: counts[name] for name, key in names.items()}}


def cynog_catalog_facts(facts_dir):
    """CynOGv3 catalog counts from the catalog build records."""
    genomes = _rows(facts_dir / "cynog_phase0_genome_summary.tsv")
    families = _rows(facts_dir / "cynog_summary.tsv")
    log = (facts_dir / "cynog_hmm_build.log").read_text()
    built = re.findall(r"^\s*HMMs built: (\d+)$", log, re.M)
    scan = re.findall(r"^\s*Self-scan: (\d+)/(\d+) correct", log, re.M)
    if len(built) != 1 or len(scan) != 1:
        raise SystemExit("supplement numbers: CynOG HMM build log lacks one profile count and self-scan line")
    return {
        "catalog_input_proteomes": len(genomes),
        "catalog_input_proteins": sum(int(r["clean_proteins"]) for r in genomes),
        "catalog_families": len(families),
        "catalog_family_proteins": sum(int(r["size"]) for r in families),
        "catalog_hmm_profiles": int(built[0]),
        "catalog_selfscan_correct": int(scan[0][0]),
        "catalog_selfscan_total": int(scan[0][1]),
    }


def screen_facts(facts_dir, domain_summary):
    """B-lectin outgroup screen composition and the raw-genome tblastn assembly count."""
    blectin = next(r for r in _rows(domain_summary) if r["family_id"] == "B_lectin")
    genomes = _rows(facts_dir / "raw_genome_tblastn_genomes.tsv")
    if len({r["accession"] for r in genomes}) != len(genomes):
        raise SystemExit("supplement numbers: duplicate assemblies in the tblastn genome manifest")
    outgroup = {taxon: int(blectin["screen_bacteria"] if taxon == "bacterium" else blectin[f"screen_{taxon}"])
                for taxon in ("bacterium", "plant", "fungus", "metazoa", "rotifer")}
    return {"b_lectin_screen": {"outgroup_taxa": outgroup,
                                "eukaryotic_outgroup": int(blectin["screen_eukaryote"])},
            "tblastn_assemblies": len(genomes)}


def expression_facts(expression_table):
    """Tissue and species rosters of the expression summary table (Table S13)."""
    rows = _rows(expression_table)
    adult_gland = {r["species"] for r in rows if r["life_stage"] == "adult" and r["tissue"] in ("venom_gland", "accessory_gland")}
    dissected = {"gut", "salivary_gland", "malpighian_tubule"}
    larval_tissue = {r["species"] for r in rows if r["life_stage"] == "larva" and r["tissue"] in dissected}
    return {
        "expression_adult_gland_species": len(adult_gland),
        "expression_larval_tissue_species": len(larval_tissue),
        "expression_larval_tissues": ", ".join(sorted({r["tissue"] for r in rows if r["life_stage"] == "larva"})),
        "expression_species": len({r["species"] for r in rows}),
    }


def catalytic_domain_facts(domain_summary):
    """GH30_3 and PI-PLC-X domain sets (accessory-domain tree inputs; F = cynipoid, R = outgroup)."""
    out = {}
    for row in _rows(domain_summary):
        family = row["family_id"]
        if family not in ("GH30_3", "PI_PLC_X"):
            continue
        out[family] = {"domains": int(row["focal_count"]) + int(row["outgroup_count"]),
                       "cynipoid": int(row["focal_count"]), "outgroup": int(row["outgroup_count"]),
                       "outgroup_taxa": {taxon: int(row["screen_bacteria"] if taxon == "bacterium" else row[f"screen_{taxon}"])
                                          for taxon in ("bacterium", "plant", "fungus", "metazoa", "rotifer")},
                       "eukaryotic_outgroup": int(row["screen_eukaryote"])}
    if set(out) != {"GH30_3", "PI_PLC_X"}:
        raise SystemExit("supplement numbers: missing catalytic domain summary")
    return out


def sequence_set_facts(family_membership, sequence_sets, table_s2, b_lectin_manifest):
    membership = _rows(family_membership)
    summary = {(r["family"], r["metric"]): r for r in _rows(sequence_sets)}
    table_rows = {r["family_id"]: r for r in _rows(table_s2)}
    medians = {}
    for row in table_rows.values():
        value = row["median_trimmed_length_aa"]
        if value:
            median = float(value)
            if not median.is_integer():
                raise SystemExit(f"supplement numbers: nonintegral trimmed median for {row['family_id']}: {value}")
            medians[row["family_id"]] = int(median)
    out = {}
    outgroup_total = 0
    outgroup_protein_total = 0
    for family in ("GH5_2", "PL1_8", "PL4_4"):
        inferred = summary[(family, "origin_inference_before_recorded_pruning")]
        total, outgroup = int(inferred["total"]), int(inferred["outgroup_count"])
        domains = [r for r in membership if r["family_id"] == family and r["unit_type"] == "domain"]
        groups = Counter(r["source_group"] for r in domains)
        outgroup_taxa = Counter(r["sequence_id"].split("|")[0] for r in domains if r["source_group"] != "focal_genome")
        if groups["bacterial_outgroup"] + groups["eukaryotic_outgroup"] != outgroup:
            raise SystemExit(f"supplement numbers: {family} outgroup membership differs from the inference set")
        parent_taxa = {}
        for taxon in outgroup_taxa:
            parent_taxa[taxon] = len({origin_parent(r['sequence_id']) for r in domains
                                     if r['source_group'] != 'focal_genome'
                                     and r['sequence_id'].split('|')[0] == taxon})
        table = table_rows[family]
        check_like_unit("domain", outgroup, "domain", int(table["outgroup_count"]), family)
        proteins = int(table["outgroup_proteins"])
        check_like_unit("protein", proteins, "protein",
                        int(table["screen_bacteria_proteins"]) + int(table["screen_eukaryote_proteins"]),
                        f"{family} parent-protein subtotal")
        if proteins > outgroup:
            raise SystemExit(f"supplement numbers: {family} has more outgroup proteins than domains")
        outgroup_total += outgroup
        outgroup_protein_total += proteins
        out[family] = {
            "analysis_domains": total, "analysis_cynipoid": total - outgroup, "outgroup": outgroup,
            "median_aa": medians[family],
            "displayed_domains": len(domains), "removed_after_inference": int(summary[(family, "origin_removed_after_inference")]["total"]),
            "bacterial_outgroup": groups["bacterial_outgroup"], "eukaryotic_outgroup": groups["eukaryotic_outgroup"],
            "outgroup_domain_tips": outgroup, "outgroup_parent_proteins": proteins,
            "screen_bacterial_parent_proteins": int(table["screen_bacteria_proteins"]),
            "screen_eukaryotic_parent_proteins": int(table["screen_eukaryote_proteins"]),
            "outgroup_taxa": dict(outgroup_taxa),
            "outgroup_parent_taxa": parent_taxa,
            "locus_tree_tips": int(summary[(family, "locus_tree")]["total"]),
        }
    blectin = [r for r in membership if r["family_id"] == "B_lectin"]
    proteins = {r["sequence_id"] for r in blectin if r["unit_type"] == "protein"}
    focal_domains = [r for r in blectin if r["unit_type"] == "domain" and r["source_group"] == "focal_genome"]
    outgroup_domains = [r for r in blectin if r["unit_type"] == "domain" and r["source_group"] != "focal_genome"]
    displayed = _rows(b_lectin_manifest)
    out["B_lectin"] = {
        "proteins": len(proteins),
        "domain_units": len({r["parent_locus_id"] for r in focal_domains}),
        "analysis_domains": len(focal_domains) + len(outgroup_domains),
        "analysis_outgroup_domains": len(outgroup_domains),
        "displayed_domains": len(displayed),
        "displayed_cynipoid_domains": sum(r["tip_class"] == "focal" for r in displayed),
        "displayed_units": len({r["biological_unit_id"] for r in displayed if r["tip_class"] == "focal"}),
    }
    out["B_lectin"]["max_domains_per_protein"] = max(Counter(r["parent_locus_id"] for r in focal_domains).values())
    out["B_lectin"]["removed_after_inference"] = out["B_lectin"]["analysis_domains"] - out["B_lectin"]["displayed_domains"]
    out["outgroup_sequences_total"] = outgroup_total
    out["outgroup_domain_tips_total"] = outgroup_total
    out["outgroup_parent_proteins_total"] = outgroup_protein_total
    out["attributed_transcripts_removed"] = sum(out[f]["removed_after_inference"] for f in ("GH5_2", "PL1_8", "PL4_4", "B_lectin"))
    return out


def gc3_facts(gc3):
    """Section F GC3 summaries over the intron producers' GC3 statistics (numbers["gc3"]).

    The per-family and B-lectin values themselves are read directly as
    {{gc3.pcwde.<family>.*}} and {{gc3.blectin.all.*}} tokens; this adds only the
    cross-family sum, ranges and P-value bound that the prose states.
    """
    families = ("GH5_2", "PL1_8", "PL4_4")
    pcwde = gc3["pcwde"]
    if set(pcwde) != set(families):
        raise SystemExit(f"supplement numbers: GC3 families {sorted(pcwde)}")
    host = [float(pcwde[f]["host_distance_pp"]) for f in families]
    bacterial = [float(pcwde[f]["bacterial_distance_pp"]) for f in families]
    p_values = [float(pcwde[f]["wilcoxon_bacteria_minus_host_distance_p"]) for f in families]
    if not all(0 < p < 1 for p in p_values):
        raise SystemExit("supplement numbers: GC3 Wilcoxon P outside (0, 1)")
    # Smallest integer exponent e with every family P strictly below 10**e.
    exponent = math.floor(math.log10(max(p_values))) + 1
    if not max(p_values) < 10 ** exponent:
        raise SystemExit("supplement numbers: GC3 P bound is not strict")
    return {
        "pcwde_copies": sum(int(pcwde[f]["copies"]) for f in families),
        "pcwde_host_distance_pp_min": min(host),
        "pcwde_host_distance_pp_max": max(host),
        "pcwde_bacterial_distance_pp_min": min(bacterial),
        "pcwde_bacterial_distance_pp_max": max(bacterial),
        "pcwde_wilcoxon_p_upper_exponent": exponent,
    }
