#!/usr/bin/env python3
"""Parse current PCWDE HyPhy results and reproduce the legacy summary contract."""

from __future__ import annotations

import csv
import json
import argparse
from pathlib import Path

FAMILIES = ("GH5_2", "PL1_8", "PL4_4")
PP = 0.9


def read_tsv(path: Path) -> list[dict]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    if fields is None:
        fields = list(rows[0])
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def load_json(path: Path) -> dict:
    data = json.loads(path.read_text())
    if data.get("runtime") != "2.5.85":
        raise SystemExit(f"{path}: expected HyPhy runtime 2.5.85, got {data.get('runtime')}")
    return data


def fubar_rows(family: str, path: Path) -> tuple[list[dict], dict]:
    data = load_json(path)
    headers = [row[0] for row in data["MLE"]["headers"]]
    content = data["MLE"]["content"]["0"]
    required = ["alpha", "beta", "Prob[alpha>beta]", "Prob[alpha<beta]"]
    indices = {name: headers.index(name) for name in required}
    bf_index = headers.index("BayesFactor[alpha<beta]")
    rows = []
    omegas = []
    counts = {"purifying": 0, "positive": 0, "neutral": 0}
    for site, values in enumerate(content, 1):
        alpha = values[indices["alpha"]]
        beta = values[indices["beta"]]
        pp_purifying = values[indices["Prob[alpha>beta]"]]
        pp_positive = values[indices["Prob[alpha<beta]"]]
        if pp_purifying > PP:
            classification = "purifying"
        elif pp_positive > PP:
            classification = "positive"
        else:
            classification = "neutral"
        counts[classification] += 1
        omega = beta / alpha if alpha > 0.001 else None
        if omega is not None:
            omegas.append(omega)
        rows.append({
            "family": family, "site": site, "alpha": alpha, "beta": beta,
            "beta_minus_alpha": beta - alpha, "pp_purifying": pp_purifying,
            "pp_positive": pp_positive, "bayes_factor_positive": values[bf_index],
            "omega_if_alpha_gt_0.001": "" if omega is None else omega,
            "classification_pp_gt_0.9": classification,
        })
    omegas.sort()
    median = omegas[len(omegas) // 2] if omegas else 0
    summary = {
        "sites": len(rows), "purifying": counts["purifying"],
        "positive": counts["positive"], "neutral": counts["neutral"],
        "median_omega": median,
    }
    return rows, summary


def absrel_rows(family: str, path: Path, tip_ids: set[str]) -> tuple[list[dict], dict]:
    data = load_json(path)
    attrs = data["branch attributes"]["0"]
    rows = []
    for branch, values in attrs.items():
        original = values.get("original name", branch)
        is_tip = branch in tip_ids or original in tip_ids
        corrected = values.get("Corrected P-value")
        rates = values.get("Rate Distributions", [])
        rows.append({
            "family": family, "branch": branch, "original_name": original,
            "branch_type": "terminal" if is_tip else "internal",
            "corrected_p_value": corrected,
            "uncorrected_p_value": values.get("Uncorrected P-value"),
            "significant_corrected_p_lt_0.05": corrected is not None and corrected < 0.05,
            "lrt": values.get("LRT"), "rate_classes": values.get("Rate classes"),
            "rate_distributions_json": json.dumps(rates, separators=(",", ":")),
        })
    tips = [row for row in rows if row["branch_type"] == "terminal"]
    summary = {
        "absrel_significant_terminal": sum(row["significant_corrected_p_lt_0.05"] for row in tips),
        "absrel_terminal_tested": len(tips),
        "absrel_significant_all": sum(row["significant_corrected_p_lt_0.05"] for row in rows),
        "absrel_all_tested": len(rows),
    }
    return rows, summary


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Summarize the reported FUBAR and aBSREL analyses."
    )
    parser.add_argument(
        "--input-dir", required=True, type=Path,
        help="Directory containing selection_membership_lock.tsv and outputs/<family> JSON files.",
    )
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    input_dir = args.input_dir.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=False)
    membership = read_tsv(input_dir / "selection_membership_lock.tsv")
    all_fubar = []
    all_absrel = []
    summary_rows = []
    for family in FAMILIES:
        tip_ids = {
            row["analysis_id"] for row in membership
            if row["family"] == family and row["selection_status"] == "included"
        }
        root = input_dir / "outputs" / family
        fubar, fubar_summary = fubar_rows(family, root / f"{family}_fubar.json")
        absrel, absrel_summary = absrel_rows(family, root / f"{family}_absrel.json", tip_ids)
        if absrel_summary["absrel_terminal_tested"] != len(tip_ids):
            raise SystemExit(
                f"{family}: aBSREL terminal branches {absrel_summary['absrel_terminal_tested']} != tips {len(tip_ids)}"
            )
        summary_rows.append({
            "family": family, **fubar_summary, **absrel_summary,
            "selection_sequences": len(tip_ids),
        })
        all_fubar.extend(fubar)
        all_absrel.extend(absrel)

    write_tsv(output_dir / "selection_family_summary.tsv", summary_rows)
    write_tsv(output_dir / "fubar_site_results.tsv", all_fubar)
    write_tsv(output_dir / "absrel_branch_results.tsv", all_absrel)
    print(json.dumps(summary_rows, indent=2))


if __name__ == "__main__":
    main()
