#!/usr/bin/env python3
"""
L4.13 Analysis -- Phase 1: Strain Inventory
===========================================

Scans the SPDI database (external/spdi_db/L4.13.1/ and L4.13.2/) for strains with a valid spdi.txt,
extracts QC metrics from TBannotator report.json, and writes
data/strains.csv with one row per strain.

Columns produced (locally available):
    strain_name, sub_lineage, n_spdis, n_missing_genes, n_missing_rd,
    mean_depth, covered_bases_percent, mean_mapq, mean_baseq, gc_content

Columns left empty (to be filled by phase1b_ena_metadata.py):
    country, collection_date, host, bioproject, biosample, platform

Usage:
    python phase1_strain_inventory.py
"""

import csv
import json
import sys
from collections import Counter
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent
REPO_ROOT = PROJECT_DIR.parent

BDD = Path(__file__).resolve().parent.parent / "external" / "spdi_db"  # per-strain SPDI profiles; see external/README.md
SUB_CLADES = ["L4.13.1", "L4.13.2"]
DATA_DIR = PROJECT_DIR / "data"
RESULTS_DIR = PROJECT_DIR / "résultats"

REFERENCE = "NC_000962.3"


def list_valid_strains(directory):
    valid = []
    if not directory.is_dir():
        return valid
    for entry in directory.iterdir():
        if entry.is_dir():
            if (entry / REFERENCE / "spdi.txt").is_file():
                valid.append(entry.name)
    return sorted(valid)


def count_spdis(spdi_file):
    with open(spdi_file) as f:
        return sum(1 for line in f if line.strip())


def load_report_qc(report_path):
    qc = {
        "n_missing_genes": "",
        "n_missing_rd": "",
        "mean_depth": "",
        "covered_bases_percent": "",
        "mean_mapq": "",
        "mean_baseq": "",
        "gc_content": "",
    }
    if not report_path.is_file():
        return qc
    try:
        with open(report_path) as f:
            d = json.load(f)
    except (json.JSONDecodeError, OSError):
        return qc

    ms = d.get("mapping_stats", {}) or {}
    qc["mean_depth"] = ms.get("mean_depth", "")
    qc["covered_bases_percent"] = ms.get("covered_bases_percent", "")
    qc["mean_mapq"] = ms.get("mean_mapq", "")
    qc["mean_baseq"] = ms.get("mean_baseq", "")

    af = (d.get("quality", {}) or {}).get("after_filtering", {}) or {}
    qc["gc_content"] = af.get("gc_content", "")

    qc["n_missing_genes"] = len(d.get("missing_genes", []) or [])
    qc["n_missing_rd"] = len(d.get("missing_rd", []) or [])
    return qc


def main():
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("L4.13 Analysis -- Phase 1: Strain Inventory")
    print("=" * 70)

    rows = []
    n_total = 0
    n_with_report = 0

    for sub in SUB_CLADES:
        sub_dir = BDD / sub
        print(f"\n[scan] {sub_dir}")
        strains = list_valid_strains(sub_dir)
        print(f"  {len(strains)} strain(s) with spdi.txt")
        n_total += len(strains)

        for sra in strains:
            spdi_file = sub_dir / sra / REFERENCE / "spdi.txt"
            report_file = sub_dir / sra / REFERENCE / "report.json"
            n_spdis = count_spdis(spdi_file)
            qc = load_report_qc(report_file)
            if report_file.is_file():
                n_with_report += 1
            rows.append({
                "strain_name": sra,
                "sub_lineage": sub,
                "n_spdis": n_spdis,
                **qc,
                "country": "",
                "collection_date": "",
                "host": "",
                "bioproject": "",
                "biosample": "",
                "platform": "",
            })

    print(f"\n[inventory] {n_total} total strains, {n_with_report} with report.json")

    fieldnames = [
        "strain_name", "sub_lineage", "n_spdis",
        "n_missing_genes", "n_missing_rd",
        "mean_depth", "covered_bases_percent", "mean_mapq", "mean_baseq",
        "gc_content",
        "country", "collection_date", "host",
        "bioproject", "biosample", "platform",
    ]

    csv_path = DATA_DIR / "strains.csv"
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        for row in sorted(rows, key=lambda r: (r["sub_lineage"], r["strain_name"])):
            w.writerow(row)
    print(f"[write] {csv_path} ({len(rows)} rows)")

    # Summary stats
    print("\n[stats]")
    by_sub = Counter(r["sub_lineage"] for r in rows)
    for sub, n in by_sub.most_common():
        print(f"  {sub}: {n}")

    spdi_counts = [r["n_spdis"] for r in rows if isinstance(r["n_spdis"], int)]
    if spdi_counts:
        spdi_counts.sort()
        n = len(spdi_counts)
        print(f"\n  SPDI counts (all): min={spdi_counts[0]}, "
              f"median={spdi_counts[n // 2]}, max={spdi_counts[-1]}, "
              f"mean={sum(spdi_counts) / n:.1f}")

    for sub in SUB_CLADES:
        sub_counts = [r["n_spdis"] for r in rows if r["sub_lineage"] == sub]
        if sub_counts:
            sub_counts.sort()
            n = len(sub_counts)
            print(f"  SPDI counts ({sub}): min={sub_counts[0]}, "
                  f"median={sub_counts[n // 2]}, max={sub_counts[-1]}, "
                  f"mean={sum(sub_counts) / n:.1f}")

    depths = [r["mean_depth"] for r in rows
              if isinstance(r["mean_depth"], (int, float))]
    if depths:
        depths.sort()
        n = len(depths)
        print(f"  mean_depth: min={depths[0]:.1f}, median={depths[n // 2]:.1f}, "
              f"max={depths[-1]:.1f}")

    covs = [r["covered_bases_percent"] for r in rows
            if isinstance(r["covered_bases_percent"], (int, float))]
    if covs:
        covs.sort()
        n = len(covs)
        print(f"  covered_bases_percent: min={covs[0]:.4f}, "
              f"median={covs[n // 2]:.4f}, max={covs[-1]:.4f}")

    summary_path = RESULTS_DIR / "phase1_summary.txt"
    with open(summary_path, "w") as f:
        f.write("L4.13 Analysis -- Phase 1 Summary\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Total strains: {n_total}\n")
        f.write(f"With report.json: {n_with_report}\n\n")
        for sub, n in by_sub.most_common():
            f.write(f"  {sub}: {n}\n")
        if spdi_counts:
            f.write(f"\nSPDI counts (all): min={spdi_counts[0]}, "
                    f"median={spdi_counts[len(spdi_counts) // 2]}, "
                    f"max={spdi_counts[-1]}, "
                    f"mean={sum(spdi_counts) / len(spdi_counts):.1f}\n")
    print(f"[write] {summary_path}")

    print("\n" + "=" * 70)
    print("Phase 1 complete.")
    print("=" * 70)


if __name__ == "__main__":
    main()
