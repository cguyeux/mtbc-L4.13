#!/usr/bin/env python3
"""
L4.13 Analysis -- Phase 1b: ENA metadata retrieval
===================================================

Reads the 258 SRA accessions from data/strains.csv, queries the ENA Portal
API in batches, and writes data/ena_metadata.csv. Also merges the new
columns back into data/strains.csv (country, collection_date, host,
host_scientific_name, isolation_source, bioproject, biosample,
instrument_platform, instrument_model, library_strategy, read_count,
base_count, first_public).

Usage:
    python phase1b_ena_metadata.py
"""

import csv
import sys
import time
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent
DATA_DIR = PROJECT_DIR / "data"
RESULTS_DIR = PROJECT_DIR / "résultats"

STRAINS_CSV = DATA_DIR / "strains.csv"
ENA_CSV = DATA_DIR / "ena_metadata.csv"

ENA_FIELDS = [
    "run_accession", "sample_accession", "study_accession",
    "country", "collection_date", "host", "host_scientific_name",
    "isolation_source",
    "instrument_platform", "instrument_model", "library_strategy",
    "read_count", "base_count", "first_public",
]

ENA_URL = "https://www.ebi.ac.uk/ena/portal/api/search"
BATCH_SIZE = 50


def query_batch(accessions):
    """POST-style query to ENA search; returns list of dict rows."""
    clauses = " OR ".join(f'run_accession="{a}"' for a in accessions)
    params = {
        "result": "read_run",
        "query": clauses,
        "fields": ",".join(ENA_FIELDS),
        "format": "tsv",
        "limit": "0",
    }
    data = urllib.parse.urlencode(params).encode("utf-8")
    req = urllib.request.Request(ENA_URL, data=data, method="POST")
    with urllib.request.urlopen(req, timeout=60) as resp:
        text = resp.read().decode("utf-8")
    lines = text.strip().split("\n")
    if not lines:
        return []
    header = lines[0].split("\t")
    rows = []
    for line in lines[1:]:
        cells = line.split("\t")
        if len(cells) < len(header):
            cells += [""] * (len(header) - len(cells))
        rows.append(dict(zip(header, cells)))
    return rows


def main():
    if not STRAINS_CSV.is_file():
        print(f"ERROR: {STRAINS_CSV} not found. Run phase1 first.")
        sys.exit(1)

    print("=" * 70)
    print("L4.13 Analysis -- Phase 1b: ENA metadata retrieval")
    print("=" * 70)

    # Load strains
    with open(STRAINS_CSV) as f:
        reader = csv.DictReader(f)
        strains = list(reader)
        base_fieldnames = reader.fieldnames
    accessions = [s["strain_name"] for s in strains]
    print(f"\n[load] {len(accessions)} accessions from {STRAINS_CSV}")

    # Query ENA in batches
    ena_rows = {}
    total_batches = (len(accessions) + BATCH_SIZE - 1) // BATCH_SIZE
    for i in range(0, len(accessions), BATCH_SIZE):
        batch = accessions[i:i + BATCH_SIZE]
        batch_idx = i // BATCH_SIZE + 1
        print(f"[query] batch {batch_idx}/{total_batches} "
              f"({len(batch)} accessions)...", end=" ")
        try:
            results = query_batch(batch)
            for r in results:
                ena_rows[r["run_accession"]] = r
            print(f"got {len(results)} records")
        except Exception as e:
            print(f"ERROR: {e}")
        time.sleep(0.3)

    print(f"\n[merge] {len(ena_rows)}/{len(accessions)} accessions with ENA data")

    # Write ena_metadata.csv
    with open(ENA_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=ENA_FIELDS, extrasaction="ignore")
        w.writeheader()
        for acc in accessions:
            row = ena_rows.get(acc, {"run_accession": acc})
            w.writerow({k: row.get(k, "") for k in ENA_FIELDS})
    print(f"[write] {ENA_CSV}")

    # Merge into strains.csv (update existing columns only)
    col_map = {
        "country": "country",
        "collection_date": "collection_date",
        "host": "host",
        "bioproject": "study_accession",
        "biosample": "sample_accession",
        "platform": "instrument_platform",
    }
    for s in strains:
        acc = s["strain_name"]
        ena = ena_rows.get(acc, {})
        for strains_col, ena_col in col_map.items():
            if strains_col in base_fieldnames and ena_col in ena:
                s[strains_col] = ena.get(ena_col, "")

    with open(STRAINS_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=base_fieldnames, extrasaction="ignore")
        w.writeheader()
        for s in strains:
            w.writerow(s)
    print(f"[update] {STRAINS_CSV} with 6 ENA columns merged")

    # Summary
    print("\n[stats]")
    c = Counter((s.get("country") or "(unknown)").split(":")[0].strip()
                for s in strains)
    print(f"\nCountry distribution ({len([x for x in c if x != '(unknown)'])} known):")
    for country, n in c.most_common(25):
        print(f"  {country or '(empty)':<40s} {n}")

    years = Counter()
    for s in strains:
        d = s.get("collection_date") or ""
        # Crude year extraction
        if len(d) >= 4 and d[:4].isdigit():
            years[d[:4]] += 1
        else:
            years["(unknown)"] += 1
    print(f"\nYear distribution (top 15):")
    for y, n in sorted(years.items()):
        print(f"  {y}: {n}")

    hosts = Counter((s.get("host") or "(unknown)") for s in strains)
    print(f"\nHost distribution:")
    for h, n in hosts.most_common(10):
        print(f"  {h or '(empty)':<40s} {n}")

    bps = Counter(s.get("bioproject") or "(unknown)" for s in strains)
    print(f"\nTop 10 BioProjects:")
    for bp, n in bps.most_common(10):
        print(f"  {bp or '(empty)':<20s} {n}")

    # Per-sub-clade country breakdown
    for sub in ["L4.13.1", "L4.13.2"]:
        c_sub = Counter((s.get("country") or "(unknown)").split(":")[0].strip()
                        for s in strains if s.get("sub_lineage") == sub)
        print(f"\n{sub} top 10 countries:")
        for country, n in c_sub.most_common(10):
            print(f"  {country or '(empty)':<40s} {n}")

    # Write summary file
    summary_path = RESULTS_DIR / "phase1b_summary.txt"
    with open(summary_path, "w") as f:
        f.write("L4.13 -- Phase 1b: ENA metadata summary\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Total strains: {len(strains)}\n")
        f.write(f"With ENA record: {len(ena_rows)}\n\n")
        f.write("Country distribution:\n")
        for country, n in c.most_common():
            f.write(f"  {country or '(empty)':<40s} {n}\n")
        f.write("\nYear distribution:\n")
        for y, n in sorted(years.items()):
            f.write(f"  {y}: {n}\n")
        f.write("\nHost distribution:\n")
        for h, n in hosts.most_common():
            f.write(f"  {h or '(empty)':<40s} {n}\n")
        f.write("\nBioProject distribution (top 30):\n")
        for bp, n in bps.most_common(30):
            f.write(f"  {bp or '(empty)':<20s} {n}\n")
    print(f"\n[write] {summary_path}")

    print("\n" + "=" * 70)
    print("Phase 1b complete.")
    print("=" * 70)


if __name__ == "__main__":
    main()
