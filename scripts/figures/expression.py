#!/usr/bin/env python3
"""Render the accepted current Fig. 7 tissue-parity design.

DC: the eight reviewed Figure 2 transcript tips are drawn as individual tissue
profiles and labeled with their paired BK neighbor label plus an asterisk. All
unpaired DC family members remain visible as one unlabeled gray pooled group.
BK: F003 gut, Malpighian tubule, salivary gland plus the existing adult venom
gland, carcass, and accessory-gland series. Whole-body/whole-larva, head, and
cuticle/pelt are excluded.

For the three focal PCWDE families, the exact reviewed Figure 2 label manifest
is the single display authority: eight DC tips inherit the already paired BK
neighbor label with an asterisk, and nine BK genomic loci retain their own
chromosome labels. Plot each selected BK locus across all six tissues; place
all remaining BK genes in one ``other`` column with gene-by-replicate points
left separate (not summed). Direct labels are used instead of mapping lines or
a new tree-neighbor calculation.
"""

import importlib.util
import os
from collections import defaultdict
from pathlib import Path

import matplotlib.colors as mcolors
import matplotlib as mpl
import numpy as np
import pandas as pd


mpl.rcParams["svg.hashsalt"] = "accepted-figure8-f003-parity"

DATA_ROOT_VALUE = os.environ.get("PCWDE_FIGURE_DATA_ROOT", "")
if not DATA_ROOT_VALUE:
    raise SystemExit("PCWDE_FIGURE_DATA_ROOT is required")
DATA_ROOT = Path(DATA_ROOT_VALUE).resolve()
if not DATA_ROOT.is_dir():
    raise SystemExit(f"PCWDE_FIGURE_DATA_ROOT is not a directory: {DATA_ROOT}")
FIGURE_INPUT = DATA_ROOT / "expression"
SHARED_INPUT = DATA_ROOT / "shared"
OUTPUT_VALUE = os.environ.get("PCWDE_EXPRESSION_OUTPUT_DIR", "")
if not OUTPUT_VALUE:
    raise SystemExit("PCWDE_EXPRESSION_OUTPUT_DIR is required")
OUT = Path(OUTPUT_VALUE).resolve()
SOURCE = Path(__file__).with_name("expression_base.py")
F003_LONG = FIGURE_INPUT / "f003_expression_long.tsv"
REGION_STATUS = FIGURE_INPUT / "region_status.tsv"
FIG2_FOCAL = FIGURE_INPUT / "reviewed_labels.tsv"
FIG2_LABEL_CROSSWALK = FIGURE_INPUT / "raw_tip_to_normalized_display_label.tsv"

# Bind the path-adapted base scaffold before importing it. All evidence inputs
# therefore resolve below the one Figshare-style public data root.
os.environ["PCWDE_EXPRESSION_INPUT_DIR"] = str(FIGURE_INPUT)
os.environ["PCWDE_FIGURE_SHARED_INPUT_DIR"] = str(SHARED_INPUT)
os.environ["PCWDE_EXPRESSION_MEMBER_METADATA"] = str(DATA_ROOT.parent / "tables/authority/table_s2/metadata.tsv")

TOP_PARALOG_MEAN_TPM_MIN = 100.0
FOCAL_FAMILIES = {"PL1_8", "PL4_4", "GH5_2"}
BK_GENOME_KEY = "Belonocnema_kinseyi_GCF"

CM_TO_CHR = {
    "CM021338.1": "chr1", "CM021339.1": "chr2", "CM021340.1": "chr3",
    "CM021341.1": "chr4", "CM021342.1": "chr5", "CM021343.1": "chr6",
    "CM021344.1": "chr7", "CM021345.1": "chr8", "CM021346.1": "chr9",
    "CM021347.1": "chr10",
}


def chromosome(scaffold):
    scaffold = str(scaffold)
    for accession, label in CM_TO_CHR.items():
        if accession in scaffold:
            return label
    return "frag"


def chromosome_sort_key(label):
    if label.startswith("chr") and label[3:].isdigit():
        return (0, int(label[3:]))
    return (1, label)


def build_bk_synteny_labels():
    """Reproduce the chr#letter naming rule from 07_family_summary_matrix.py."""
    status = pd.read_csv(REGION_STATUS, sep="\t", keep_default_na=False)
    protein_to_label = {}
    protein_to_paralog_label = {}
    protein_to_kind = {}
    protein_to_region = {}
    audit = []
    for family in sorted(FOCAL_FAMILIES):
        fam = status[status["family"].eq(family)].copy()
        items = defaultdict(list)

        blocks = fam[(fam["region_type"].eq("block")) & fam["block_id"].ne("SINGLETON")]
        for block_id, group in blocks.groupby("block_id"):
            present = group[group["status"].eq("PRESENT")]
            owner = present.iloc[0] if not present.empty else group.iloc[0]
            items[chromosome(owner["scaffold"])].append(("block", block_id))

        singletons = fam[fam["region_type"].eq("singleton")]
        for singleton_id, group in singletons.groupby("block_id"):
            present = group[group["status"].eq("PRESENT")]
            owner = present.iloc[0] if not present.empty else group.iloc[0]
            items[chromosome(owner["scaffold"])].append(("singleton", singleton_id))

        id_to_label = {}
        for chrom in sorted(items, key=chromosome_sort_key):
            ordered = sorted(items[chrom], key=lambda x: (0 if x[0] == "block" else 1, x[1]))
            for index, (_kind, region_id) in enumerate(ordered):
                if index >= 26:
                    raise ValueError(f"More than 26 {family} regions on {chrom}")
                id_to_label[region_id] = f"{chrom}{chr(ord('a') + index)}"

        bk = fam[(fam["genome"].eq(BK_GENOME_KEY)) & fam["status"].eq("PRESENT")]
        for row in bk.itertuples(index=False):
            if not row.target_gene_ids:
                continue
            label = id_to_label[row.block_id]
            region_kind = row.region_type
            target_ids = row.target_gene_ids.split(",")
            for copy_index, protein_id in enumerate(target_ids, start=1):
                prior_kind = protein_to_kind.get(protein_id)
                # region_status currently duplicates the BK PL4_4 chr1 locus as
                # both supported BLOCK_0017 and a singleton. Prefer the
                # cross-species block, but retain both rows in the audit.
                use_assignment = prior_kind is None or (
                    prior_kind == "singleton" and region_kind == "block"
                )
                if use_assignment:
                    protein_to_label[protein_id] = label
                    protein_to_paralog_label[protein_id] = (
                        label if len(target_ids) == 1 else f"{label}.{copy_index}"
                    )
                    protein_to_kind[protein_id] = region_kind
                    protein_to_region[protein_id] = row.block_id
                audit.append((family, protein_id, row.block_id, region_kind, label,
                              label if len(target_ids) == 1 else f"{label}.{copy_index}"))
    audit_df = pd.DataFrame(
        audit, columns=["family", "protein_id", "region_id", "region_type",
                        "synteny_label", "paralog_label"]
    )
    audit_df["used_for_figure"] = [
        protein_to_region.get(row.protein_id) == row.region_id
        and protein_to_label.get(row.protein_id) == row.synteny_label
        for row in audit_df.itertuples(index=False)
    ]
    return protein_to_label, protein_to_paralog_label, audit_df


spec = importlib.util.spec_from_file_location("canonical_fig7", SOURCE)
fig7 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fig7)
fig7.OUTPUT_DIR = OUT

fig7.DC_TISSUES = {
    "Gut": fig7.DC_TISSUES["Gut"],
    "Malp.": fig7.DC_TISSUES["Malp."],
    "Saliv.": fig7.DC_TISSUES["Saliv."],
    "Venom": fig7.DC_TISSUES["Venom"],
    "Carc.": fig7.DC_TISSUES["Carc."],
    "Acc.": fig7.DC_TISSUES["Acc."],
}
fig7.DC_ORDER = list(fig7.DC_TISSUES)
fig7.DC_LARVAL_N = 3

f003 = pd.read_csv(F003_LONG, sep="\t", keep_default_na=False)
f003 = f003[f003["family"].isin(FOCAL_FAMILIES | {"GH1_e72", "CBM13"})].copy()
f003_samples = f003[["sample_id", "tissue_type"]].drop_duplicates()
f003_tissues = [
    ("Gut", "gut"),
    ("Malp.", "malpighian_tubule"),
    ("Saliv.", "salivary_gland"),
]
fig7.BK_TISSUES = {
    label: f003_samples.loc[f003_samples["tissue_type"].eq(tissue), "sample_id"].tolist()
    for label, tissue in f003_tissues
}
fig7.BK_TISSUES.update({
    "Venom": ["tpm_BKS_VGX_02", "tpm_BKS_VGX_03", "tpm_BKS_VGX_04", "tpm_BKS_VGX_06", "tpm_BKS_VGX_08"],
    "Carc.": ["tpm_BKS_CRX_02", "tpm_BKS_CRX_03", "tpm_BKS_CRX_04", "tpm_BKS_CRX_06", "tpm_BKS_CRX_08"],
    "Acc.": ["tpm_BKS_AGX_02", "tpm_BKS_AGX_03", "tpm_BKS_AGX_04", "tpm_BKS_AGX_06", "tpm_BKS_AGX_08"],
})
fig7.BK_ORDER = list(fig7.BK_TISSUES)
fig7.BK_LARVAL_N = 3

# One reviewed label authority is shared with Figure 2. The tree-tip manifest
# is joined to the canonical tip/protein crosswalk rather than translated by a
# second hand-maintained dictionary.
focal_labels = pd.read_csv(FIG2_FOCAL, sep="\t", keep_default_na=False)
tip_crosswalk = pd.read_csv(FIG2_LABEL_CROSSWALK, sep="\t", keep_default_na=False)
focal_labels = focal_labels.merge(
    tip_crosswalk[["family", "raw_tree_tip", "tree_protein_id"]],
    on=["family", "raw_tree_tip"], how="left", validate="one_to_one",
)
if focal_labels["tree_protein_id"].eq("").any() or len(focal_labels) != 17:
    raise ValueError("Figure 2 focal-label manifest/crosswalk join failed")
focal_labels["short_label"] = focal_labels["display_label"].str.split(" — ").str[-1]
selected_bk = {
    family: set(group["short_label"])
    for family, group in focal_labels[focal_labels["highlight_group"].eq("BK")].groupby("family")
}
selected_dc = {}
pair_rows = []
for family, group in focal_labels.groupby("family", sort=False):
    pending_dc = None
    selected_dc[family] = {}
    for row in group.itertuples(index=False):
        if row.highlight_group == "DC":
            if pending_dc is not None:
                raise ValueError(f"Unpaired consecutive DC Figure 2 labels in {family}")
            pending_dc = row
        elif row.highlight_group == "BK" and pending_dc is not None:
            transferred = f"{row.short_label}*"
            selected_dc[family][pending_dc.tree_protein_id] = transferred
            pair_rows.append({
                "family": family,
                "dc_tree_protein_id": pending_dc.tree_protein_id,
                "figure2_dc_label": pending_dc.short_label,
                "paired_bk_label": row.short_label,
                "figure7_dc_display_label": transferred,
                "pair_authority": "adjacent reviewed pair in selected_tip_labels.tsv",
            })
            pending_dc = None
    if pending_dc is not None:
        raise ValueError(f"Final DC Figure 2 label lacks a paired BK row in {family}")
if len(pair_rows) != 8 or sum(len(labels) for labels in selected_dc.values()) != 8:
    raise ValueError(f"Expected eight reviewed DC/BK pairs, observed {len(pair_rows)}")

f003_pivot = f003.pivot_table(
    index=["family", "paper_protein_id"], columns="sample_id", values="TPM",
    aggfunc="sum", fill_value=0,
)
synteny_map, paralog_label_map, synteny_audit = build_bk_synteny_labels()

original_load = fig7.load_data


def load_with_f003_and_synteny():
    family_data = original_load()
    family_key = {"nonHGT": "GH1_e72", "CBM13": "CBM13"}
    for key, frame in family_data.items():
        family = family_key.get(key, key)
        frame = frame.copy()
        # The combined DC/BK frame carries the BK-only `_pid_stripped` column;
        # DC rows are null in that column. Normalize the populated per-row
        # protein identifier directly so both species retain exact member IDs.
        protein_ids = frame["protein_id"].astype(str).str.replace(r"\.\d+$", "", regex=True)
        for sample_id in f003_samples["sample_id"]:
            frame[sample_id] = [
                float(f003_pivot.loc[(family, pid), sample_id])
                if species == "BK" and (family, pid) in f003_pivot.index else 0.0
                for species, pid in zip(frame["_species"], protein_ids)
            ]
        frame["_display_family"] = family
        frame["_display_protein"] = protein_ids.values
        frame["_synteny_label"] = protein_ids.map(synteny_map).fillna("unmapped")
        frame["_paralog_label"] = protein_ids.map(paralog_label_map).fillna("unmapped")
        frame["_paralog_unit"] = frame.get("gene", frame["_display_protein"]).fillna("")
        frame.loc[frame["_paralog_unit"].eq(""), "_paralog_unit"] = frame["_display_protein"]
        # RefSeq may expose multiple protein/transcript isoforms for one genomic
        # gene while the synteny pipeline retains one representative protein.
        # Propagate the representative's region label to its gene-mates.
        bk_labeled = frame[
            frame["_species"].eq("BK") & frame["_synteny_label"].ne("unmapped")
        ]
        gene_labels = bk_labeled.groupby("_paralog_unit")["_synteny_label"].agg(
            lambda values: sorted(set(values))
        )
        conflicting = gene_labels[gene_labels.map(len).gt(1)]
        if not conflicting.empty:
            raise ValueError(f"Conflicting synteny labels within BK genes: {conflicting.to_dict()}")
        gene_label_map = gene_labels.map(lambda values: values[0]).to_dict()
        bk_copy_labeled = frame[
            frame["_species"].eq("BK") & frame["_paralog_label"].ne("unmapped")
        ]
        gene_copy_labels = bk_copy_labeled.groupby("_paralog_unit")["_paralog_label"].agg(
            lambda values: sorted(set(values))
        )
        if gene_copy_labels.map(len).gt(1).any():
            raise ValueError(
                f"Conflicting within-block labels: "
                f"{gene_copy_labels[gene_copy_labels.map(len).gt(1)].to_dict()}"
            )
        gene_copy_label_map = gene_copy_labels.map(lambda values: values[0]).to_dict()
        missing = frame["_species"].eq("BK") & frame["_synteny_label"].eq("unmapped")
        frame.loc[missing, "_synteny_label"] = frame.loc[missing, "_paralog_unit"].map(
            gene_label_map
        ).fillna("unmapped")
        frame.loc[missing, "_paralog_label"] = frame.loc[missing, "_paralog_unit"].map(
            gene_copy_label_map
        ).fillna("unmapped")
        family_data[key] = frame
    return family_data


winner_rows = []
selected_paralog_rows = []
rendered_label_rows = []


def label_offset(index, count):
    """Keep endpoint labels inside their axes while preserving short leaders."""
    if index == 0:
        return (4, 10), "left"
    if index == count - 1:
        return (-4, 10), "right"
    return (0, 10), "center"


def color_shades(base_color, count):
    rgb = np.array(mcolors.to_rgb(base_color))
    if count == 1:
        return [base_color]
    return [mcolors.to_hex(rgb * (0.55 + 0.45 * i / (count - 1))) for i in range(count)]


def draw_bk_tissue_winners(ax, frame, tissue_map, tissue_order, family_color):
    """Show the reviewed Figure 2 BK loci plus a split-point other group."""
    bk = frame[frame["_species"].eq("BK")].copy()
    sample_columns = sorted({
        col for columns in tissue_map.values() for col in columns if col in bk.columns
    })
    membership = (
        bk.groupby(["_paralog_unit", "_synteny_label", "_paralog_label"], as_index=False)
        .agg(
            _display_family=("_display_family", "first"),
            _display_proteins=("_display_protein", lambda values: ",".join(sorted(set(values)))),
        )
    )
    expression = bk.groupby(
        ["_paralog_unit", "_synteny_label", "_paralog_label"], as_index=False
    )[sample_columns].sum()
    bk = expression.merge(
        membership, on=["_paralog_unit", "_synteny_label", "_paralog_label"],
        how="left", validate="one_to_one"
    )
    selected_units = {}
    for tissue in tissue_order:
        cols = [col for col in tissue_map[tissue] if col in bk.columns]
        if not cols or bk.empty:
            continue
        means = bk[cols].mean(axis=1)
        maxima = bk[cols].max(axis=1)
        mean_winner_idx = means.idxmax()
        max_winner_idx = maxima.idxmax()
        row = bk.loc[mean_winner_idx]
        winner_mean = float(means.loc[mean_winner_idx])
        retained = winner_mean > TOP_PARALOG_MEAN_TPM_MIN
        if retained:
            prior = selected_units.get(row["_paralog_unit"])
            if prior is None or winner_mean > prior["peak_mean_tpm"]:
                selected_units[row["_paralog_unit"]] = {
                    "peak_tissue": tissue,
                    "peak_mean_tpm": winner_mean,
                }
        winner_rows.append({
            "family": row["_display_family"], "tissue": tissue,
            "mean_winner_paralog_unit": row["_paralog_unit"],
            "winner_member_proteins": row["_display_proteins"],
            "synteny_label": row["_synteny_label"],
            "winner_mean_tpm": winner_mean,
            "winner_max_tpm": float(maxima.loc[mean_winner_idx]),
            "max_based_winner_paralog_unit": bk.loc[max_winner_idx, "_paralog_unit"],
            "mean_max_winner_agree": mean_winner_idx == max_winner_idx,
            "retained_over_100_mean_tpm": retained,
        })

    # Replace the former per-tissue winner display rule with the exact reviewed
    # Figure 2 locus set. Tissue winners are still calculated above as an audit,
    # but do not silently change which named loci appear in the paper figure.
    family = str(bk["_display_family"].iloc[0])
    requested_labels = selected_bk[family]
    requested = bk[bk["_paralog_label"].isin(requested_labels)].copy()
    observed_labels = set(requested["_paralog_label"])
    if observed_labels != requested_labels:
        raise ValueError(
            f"Figure 7 BK focal-label mismatch for {family}: "
            f"requested={sorted(requested_labels)} observed={sorted(observed_labels)}"
        )
    selected_units = {}
    for _, row in requested.iterrows():
        tissue_means = {}
        for tissue in tissue_order:
            cols = [col for col in tissue_map[tissue] if col in requested.columns]
            tissue_means[tissue] = float(row[cols].astype(float).mean())
        peak_tissue = max(tissue_order, key=lambda tissue: tissue_means[tissue])
        selected_units[row["_paralog_unit"]] = {
            "peak_tissue": peak_tissue,
            "peak_mean_tpm": tissue_means[peak_tissue],
        }

    selected = bk[bk["_paralog_unit"].isin(selected_units)].copy()
    selected["_first_winning_tissue"] = selected["_paralog_unit"].map(
        lambda unit: tissue_order.index(selected_units[unit]["peak_tissue"])
    )
    selected = selected.sort_values(["_first_winning_tissue", "_synteny_label", "_paralog_unit"])
    other = bk[~bk["_paralog_unit"].isin(selected_units)].copy()
    count = len(selected) + (0 if other.empty else 1)
    colors = color_shades(family_color, max(len(selected), 1))
    span = 0.76
    offsets = [0] if count == 1 else np.linspace(-span / 2, span / 2, count)
    width = 0.54 if count == 1 else span / count * 0.80
    for (_, row), color, offset in zip(selected.iterrows(), colors, offsets):
        positions = [index + 1 + offset for index in range(len(tissue_order))]
        profile = []
        for tissue in tissue_order:
            cols = [col for col in tissue_map[tissue] if col in selected.columns]
            profile.append(row[cols].astype(float).to_numpy())
        bp = ax.boxplot(
            profile, positions=positions, widths=width, patch_artist=True,
            medianprops=dict(color=color, linewidth=1.1),
            whiskerprops=dict(color="#aaa", linewidth=0.7),
            capprops=dict(color="#aaa", linewidth=0.7), flierprops=dict(marker=""),
        )
        for patch in bp["boxes"]:
            patch.set_facecolor(fig7.lighten(color, 0.68))
            patch.set_edgecolor(color)
            patch.set_linewidth(0.7)
        for position, values in zip(positions, profile):
            x = np.random.normal(position, 0.012, size=len(values))
            ax.scatter(x, values, color=color, s=12, alpha=0.8,
                       edgecolor="white", linewidth=0.2, zorder=3)
        peak_tissue = selected_units[row["_paralog_unit"]]["peak_tissue"]
        peak_mean_tpm = selected_units[row["_paralog_unit"]]["peak_mean_tpm"]
        peak_index = tissue_order.index(peak_tissue)
        peak_values = profile[peak_index]
        peak_top = float(np.max(peak_values))
        offset_points, horizontal = label_offset(peak_index, len(tissue_order))
        if row["_paralog_label"] == "chr7a.4":
            offset_points = (offset_points[0], offset_points[1] + 7)
        ax.annotate(
            row["_paralog_label"], xy=(positions[peak_index], peak_top),
            xytext=offset_points, textcoords="offset points", ha=horizontal,
            va="bottom", fontsize=9.5, color=color, fontweight="bold",
            annotation_clip=True, zorder=8,
        )
        selected_paralog_rows.append({
            "family": row["_display_family"],
            "paralog_unit": row["_paralog_unit"],
            "member_proteins": row["_display_proteins"],
            "synteny_label": row["_synteny_label"],
            "paralog_label": row["_paralog_label"],
            "peak_tissue": peak_tissue,
            "peak_mean_tpm": peak_mean_tpm,
        })
        rendered_label_rows.append({
            "family": row["_display_family"], "species": "BK",
            "stable_unit": row["_paralog_unit"],
            "display_label": row["_paralog_label"],
            "peak_tissue": peak_tissue,
        })

    if not other.empty:
        offset = offsets[-1]
        positions = [index + 1 + offset for index in range(len(tissue_order))]
        profile = []
        for tissue in tissue_order:
            cols = [col for col in tissue_map[tissue] if col in other.columns]
            profile.append(other[cols].astype(float).to_numpy().ravel())
        other_color = "#8a8a8a"
        bp = ax.boxplot(
            profile, positions=positions, widths=width, patch_artist=True,
            medianprops=dict(color=other_color, linewidth=1.0),
            whiskerprops=dict(color="#aaa", linewidth=0.7),
            capprops=dict(color="#aaa", linewidth=0.7), flierprops=dict(marker=""),
        )
        for patch in bp["boxes"]:
            patch.set_facecolor("#e1e1e1")
            patch.set_edgecolor(other_color)
            patch.set_linewidth(0.7)
        for position, values in zip(positions, profile):
            x = np.random.normal(position, 0.012, size=len(values))
            ax.scatter(x, values, color=other_color, s=9, alpha=0.55,
                       edgecolor="white", linewidth=0.15, zorder=3)
    ax.set_xticks(range(1, len(tissue_order) + 1))
    ax.set_xticklabels(tissue_order)


original_simple = fig7.draw_simple_panel
original_grouped = fig7.draw_grouped_panel


def draw_dc_focal_overlay(ax, frame, tissue_map, tissue_order, family_color):
    """Split reviewed DC proteins from an unlabeled gray remainder, as for BK."""
    family = str(frame["_display_family"].iloc[0])
    requested = selected_dc[family]
    dc = frame[frame["_species"].eq("DC")].copy()
    sample_columns = sorted({
        col for columns in tissue_map.values() for col in columns if col in dc.columns
    })
    dc = dc.groupby("_display_protein", as_index=False)[sample_columns].sum()
    observed = set(dc.loc[dc["_display_protein"].isin(requested), "_display_protein"])
    if observed != set(requested):
        raise ValueError(
            f"Figure 7 DC focal-label mismatch for {family}: "
            f"requested={sorted(requested)} observed={sorted(observed)}"
        )
    selected = dc[dc["_display_protein"].isin(requested)].copy()
    selected["_display_label"] = selected["_display_protein"].map(requested)
    peak_rows = []
    for _, row in selected.iterrows():
        means = {}
        for tissue in tissue_order:
            cols = [col for col in tissue_map[tissue] if col in selected.columns]
            means[tissue] = float(row[cols].astype(float).mean())
        peak_tissue = max(tissue_order, key=lambda tissue: means[tissue])
        peak_rows.append((row["_display_protein"], peak_tissue, means[peak_tissue]))
    peak_map = {protein: (tissue, value) for protein, tissue, value in peak_rows}
    selected["_peak_order"] = selected["_display_protein"].map(
        lambda protein: tissue_order.index(peak_map[protein][0])
    )
    selected = selected.sort_values(["_peak_order", "_display_label", "_display_protein"])
    other = dc[~dc["_display_protein"].isin(requested)].copy()
    count = len(selected) + (0 if other.empty else 1)
    colors = color_shades(family_color, max(len(selected), 1))
    span = 0.76
    offsets = [0] if count == 1 else np.linspace(-span / 2, span / 2, count)
    width = 0.54 if count == 1 else span / count * 0.80

    for (_, row), color, offset in zip(selected.iterrows(), colors, offsets):
        positions = [index + 1 + offset for index in range(len(tissue_order))]
        profile = []
        for tissue in tissue_order:
            cols = [col for col in tissue_map[tissue] if col in selected.columns]
            profile.append(row[cols].astype(float).to_numpy())
        bp = ax.boxplot(
            profile, positions=positions, widths=width, patch_artist=True,
            medianprops=dict(color=color, linewidth=1.1),
            whiskerprops=dict(color="#aaa", linewidth=0.7),
            capprops=dict(color="#aaa", linewidth=0.7), flierprops=dict(marker=""),
        )
        for patch in bp["boxes"]:
            patch.set_facecolor(fig7.lighten(color, 0.68))
            patch.set_edgecolor(color)
            patch.set_linewidth(0.7)
        for position, values in zip(positions, profile):
            x = np.random.normal(position, 0.012, size=len(values))
            ax.scatter(x, values, color=color, s=12, alpha=0.8,
                       edgecolor="white", linewidth=0.2, zorder=3)
        peak_tissue, peak_mean = peak_map[row["_display_protein"]]
        peak_index = tissue_order.index(peak_tissue)
        peak_top = float(np.max(profile[peak_index]))
        horizontal_offset, horizontal = label_offset(peak_index, len(tissue_order))
        if row["_display_label"] == "chr10b.2*":
            # Extra clearance above the DC venom bar/point cloud.
            horizontal_offset = (horizontal_offset[0], horizontal_offset[1] + 14)
        ax.annotate(
            row["_display_label"], xy=(positions[peak_index], peak_top),
            xytext=horizontal_offset, textcoords="offset points", ha=horizontal,
            va="bottom", fontsize=9.5, color=color, fontweight="bold",
            annotation_clip=True, zorder=8,
        )
        rendered_label_rows.append({
            "family": family, "species": "DC", "stable_unit": row["_display_protein"],
            "display_label": row["_display_label"], "peak_tissue": peak_tissue,
        })

    if not other.empty:
        positions = [index + 1 + offsets[-1] for index in range(len(tissue_order))]
        profile = []
        for tissue in tissue_order:
            cols = [col for col in tissue_map[tissue] if col in other.columns]
            profile.append(other[cols].astype(float).to_numpy().ravel())
        other_color = "#8a8a8a"
        bp = ax.boxplot(
            profile, positions=positions, widths=width, patch_artist=True,
            medianprops=dict(color=other_color, linewidth=1.0),
            whiskerprops=dict(color="#aaa", linewidth=0.7),
            capprops=dict(color="#aaa", linewidth=0.7), flierprops=dict(marker=""),
        )
        for patch in bp["boxes"]:
            patch.set_facecolor("#e1e1e1")
            patch.set_edgecolor(other_color)
            patch.set_linewidth(0.7)
        for position, values in zip(positions, profile):
            x = np.random.normal(position, 0.012, size=len(values))
            ax.scatter(x, values, color=other_color, s=9, alpha=0.55,
                       edgecolor="white", linewidth=0.15, zorder=3)
    ax.set_xticks(range(1, len(tissue_order) + 1))
    ax.set_xticklabels(tissue_order)
    # Leave room for direct labels above the highest selected profile. The
    # broken GH5 upper axis is finalized separately by the scaffold below.
    y_lo, y_hi = ax.get_ylim()
    ax.set_ylim(y_lo, y_lo + 1.18 * (y_hi - y_lo))


def simple_overlay(ax, frame, species, tissue_map, tissue_order, family_color):
    if species == "BK" and frame["_display_family"].iloc[0] in FOCAL_FAMILIES:
        return draw_bk_tissue_winners(ax, frame, tissue_map, tissue_order, family_color)
    if species == "DC" and frame["_display_family"].iloc[0] in FOCAL_FAMILIES:
        return draw_dc_focal_overlay(ax, frame, tissue_map, tissue_order, family_color)
    return original_simple(ax, frame, species, tissue_map, tissue_order, family_color)


def grouped_overlay(ax, frame, species, tissue_map, tissue_order, loci, **kwargs):
    if species == "BK" and frame["_display_family"].iloc[0] == "GH5_2":
        return draw_bk_tissue_winners(ax, frame, tissue_map, tissue_order, "#1B7837")
    if species == "DC" and frame["_display_family"].iloc[0] == "GH5_2":
        return draw_dc_focal_overlay(ax, frame, tissue_map, tissue_order, "#1B7837")
    return original_grouped(ax, frame, species, tissue_map, tissue_order, loci, **kwargs)


fig7.load_data = load_with_f003_and_synteny
fig7.draw_simple_panel = simple_overlay
fig7.draw_grouped_panel = grouped_overlay
fig7.draw_nonhgt_legend = lambda ax: None
np.random.seed(7)
fig7.main()

OUT.mkdir(parents=True, exist_ok=True)
winners = pd.DataFrame(winner_rows).sort_values(["family", "tissue"])
winners.to_csv(OUT / "bk_top_paralog_per_tissue.tsv", sep="\t", index=False)
pd.DataFrame(selected_paralog_rows).drop_duplicates().sort_values(
    ["family", "synteny_label", "paralog_unit"]
).to_csv(OUT / "bk_selected_paralogs.tsv", sep="\t", index=False)
synteny_audit.to_csv(OUT / "bk_protein_synteny_labels.tsv", sep="\t", index=False)
rendered_labels = pd.DataFrame(rendered_label_rows).drop_duplicates()
expected_labels = {
    (family, "DC", display_label)
    for family, labels in selected_dc.items() for display_label in labels.values()
} | {
    (family, "BK", display_label)
    for family, labels in selected_bk.items() for display_label in labels
}
observed_labels = set(zip(
    rendered_labels["family"], rendered_labels["species"], rendered_labels["display_label"]
))
if observed_labels != expected_labels or len(rendered_labels) != 17:
    raise SystemExit(
        f"Figure 2/Figure 7 exact-label gate failed: "
        f"missing={sorted(expected_labels-observed_labels)} "
        f"extra={sorted(observed_labels-expected_labels)}"
    )
rendered_labels.sort_values(["family", "species", "display_label"]).to_csv(
    OUT / "figure2_figure7_exact_label_audit.tsv", sep="\t", index=False
)
pd.DataFrame(pair_rows).to_csv(
    OUT / "figure2_neighbor_pair_transfer.tsv", sep="\t", index=False
)

# Machine-readable accounting for the exact focal display units. The reviewed
# DC proteins are individual units; unpaired proteins are retained as an
# unlabeled pooled remainder. BK retains its native genomic-gene units.
accounting_rows = []
rendered_data = load_with_f003_and_synteny()
for key in sorted(FOCAL_FAMILIES):
    frame = rendered_data[key].copy()
    for species, tissue_map in (("DC", fig7.DC_TISSUES), ("BK", fig7.BK_TISSUES)):
        sub = frame[frame["_species"].eq(species)].copy()
        if species == "DC":
            groups = []
            for protein_id, display_label in selected_dc[key].items():
                groups.append((display_label, protein_id,
                               sub[sub["_display_protein"].eq(protein_id)]))
            unpaired = sub[~sub["_display_protein"].isin(selected_dc[key])]
            for protein_id, group in unpaired.groupby("_display_protein", sort=True):
                groups.append(("unpaired_pool", protein_id, group))
        else:
            groups = [(unit, unit, group) for unit, group in
                      sub.groupby("_paralog_unit", sort=True, dropna=False)]
        for display_unit, point_unit, group in groups:
            members = ",".join(sorted(set(group["_display_protein"].astype(str))))
            for tissue, samples in tissue_map.items():
                for sample in samples:
                    if sample not in group.columns:
                        raise ValueError(f"Missing accounting sample: {key} {species} {sample}")
                    accounting_rows.append({
                        "family": key,
                        "species": species,
                        "display_unit": display_unit,
                        "point_unit": point_unit,
                        "member_protein_ids": members,
                        "tissue": tissue,
                        "replicate": sample,
                        "tpm": float(group[sample].sum()),
                    })
accounting = pd.DataFrame(accounting_rows).sort_values(
    ["family", "species", "display_unit", "tissue", "replicate"]
)
accounting.to_csv(OUT / "figure7_expression_accounting.tsv", sep="\t", index=False)
accounting[accounting["species"].eq("DC")].to_csv(
    OUT / "dc_split_profile_accounting.tsv", sep="\t", index=False
)

if (winners["synteny_label"] == "unmapped").any():
    raise SystemExit("A retained/considered BK focal protein lacks a synteny label")
print(f"PASS: {OUT}")
print(winners.to_string(index=False))
