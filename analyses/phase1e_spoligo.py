#!/usr/bin/env python3
"""
L4.13 Analysis -- Phase 1e: in silico spoligotyping
====================================================

Uses the ESP1..ESP43 coverage data already present in TBannotator's
`report.json["known_coverage"]` to reconstruct the 43-bit spoligotype
pattern of each L4.13 strain and encode it as an octal SITVIT-style
code.

Rule for presence / absence of an ESP locus:
    present  -> known_coverage[ESP_i]['covered_bases_percent'] >= 0.5
    absent   -> below 0.5

Encoding:
    43-bit binary:  ESP1 ESP2 ... ESP43
    octal code:  14 groups of 3 bits + 1 group of 1 bit => 15 digits
                 (SpolDB4 / SITVIT convention)

Usage:
    python phase1e_spoligo.py
"""

import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent
REPO_ROOT = PROJECT_DIR.parent

BDD = Path(__file__).resolve().parent.parent / "external" / "spdi_db"  # per-strain SPDI profiles; see external/README.md
SUB_CLADES = ["L4.13.1", "L4.13.2"]
DATA_DIR = PROJECT_DIR / "data"
RESULTS_DIR = PROJECT_DIR / "résultats"
REFERENCE = "NC_000962.3"
THRESHOLD = 0.5  # coverage threshold for "present"


def binary_to_octal43(bits):
    """Convert 43-bit list of 0/1 to SITVIT-style octal code (15 digits).

    43 = 14 * 3 + 1. Classic encoding groups the first 42 bits as 14 octal
    digits and adds the last bit as a binary suffix (0 or 1). We join into
    a 15-char string.
    """
    if len(bits) != 43:
        raise ValueError(f"expected 43 bits, got {len(bits)}")
    octal = []
    for i in range(0, 42, 3):
        v = bits[i] * 4 + bits[i + 1] * 2 + bits[i + 2]
        octal.append(str(v))
    octal.append(str(bits[42]))
    return "".join(octal)


def extract_pattern(known_coverage):
    """Return (list_43_bits, list_of_missing_ESP_keys)."""
    bits = []
    missing_key = []
    for i in range(1, 44):
        key = f"ESP{i}"
        entry = known_coverage.get(key)
        if entry is None:
            bits.append(0)
            missing_key.append(key)
            continue
        pct = entry.get("covered_bases_percent", 0)
        bits.append(1 if pct >= THRESHOLD else 0)
    return bits, missing_key


def pattern_vis(bits):
    """Return ASCII: filled square for 1, dot for 0."""
    return "".join("■" if b else "□" for b in bits)


def main():
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("L4.13 Analysis -- Phase 1e: in silico spoligotyping")
    print("=" * 70)

    rows = []
    missing_kc = 0

    for sub in SUB_CLADES:
        sub_dir = BDD / sub
        strains = sorted(p.name for p in sub_dir.iterdir()
                         if p.is_dir() and (p / REFERENCE / "spdi.txt").is_file())
        print(f"\n[scan] {sub}: {len(strains)} strains")
        for sra in strains:
            rp = sub_dir / sra / REFERENCE / "report.json"
            if not rp.is_file():
                continue
            with open(rp) as f:
                d = json.load(f)
            kc = d.get("known_coverage")
            if not kc:
                missing_kc += 1
                continue
            bits, missing_keys = extract_pattern(kc)
            octal = binary_to_octal43(bits)
            rows.append({
                "strain_name": sra,
                "sub_lineage": sub,
                "pattern_binary": "".join(str(b) for b in bits),
                "pattern_visual": pattern_vis(bits),
                "octal_sitvit": octal,
                "n_present": sum(bits),
                "n_absent": 43 - sum(bits),
                "n_missing_ESP_keys": len(missing_keys),
            })

    print(f"\n[summary] {len(rows)} strains with known_coverage, "
          f"{missing_kc} without")

    # Write per-strain CSV
    csv_path = DATA_DIR / "spoligotypes.csv"
    fieldnames = ["strain_name", "sub_lineage", "octal_sitvit",
                  "pattern_binary", "pattern_visual",
                  "n_present", "n_absent", "n_missing_ESP_keys"]
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        for row in sorted(rows, key=lambda r: (r["sub_lineage"], r["strain_name"])):
            w.writerow(row)
    print(f"[write] {csv_path}")

    # Aggregate: unique octals per sub-clade
    for sub in SUB_CLADES:
        sub_rows = [r for r in rows if r["sub_lineage"] == sub]
        octs = Counter(r["octal_sitvit"] for r in sub_rows)
        print(f"\n--- {sub} ({len(sub_rows)} strains) ---")
        print(f"  Unique octal patterns: {len(octs)}")
        for oct_code, n in octs.most_common(10):
            # Find a representative pattern visual
            sample = next(r for r in sub_rows if r["octal_sitvit"] == oct_code)
            print(f"    {oct_code}  n={n:<4d}  {sample['pattern_visual']}")

    # All patterns combined
    all_octs = Counter(r["octal_sitvit"] for r in rows)
    print(f"\n--- All L4.13 ({len(rows)} strains) ---")
    print(f"  Unique octal patterns: {len(all_octs)}")
    print(f"  Top 10:")
    for oct_code, n in all_octs.most_common(10):
        sample = next(r for r in rows if r["octal_sitvit"] == oct_code)
        pct = n / len(rows) * 100
        print(f"    {oct_code}  n={n:<4d} ({pct:4.1f}%)  {sample['pattern_visual']}")

    # Ural signature check:
    # - Ural-1 (SIT-262): spacers 29-31 absent (i.e. ESP29,30,31 = 0), rest present
    # - Ural-2 (SIT-35):  spacers 29-32 absent + possibly 33
    # - Characteristic: 33-36 present (distinguishes from LAM which lacks 33-36)
    # - Characteristic: 1-28 mostly present except sometimes 19
    print("\n--- Ural signature check ---")
    print("  Ural families typically have:")
    print("    - spacers 29-31 or 29-32 ABSENT (deletion of the IS6110-RD_ural region)")
    print("    - spacers 33-36 PRESENT (distinguishes from LAM)")
    print("    - most of spacers 1-28 PRESENT")
    print()
    n_ural_like = 0
    n_29_31_absent = 0
    n_33_36_present = 0
    for r in rows:
        bits = [int(b) for b in r["pattern_binary"]]
        # ESP29..31 absent  (indices 28..30 in 0-based)
        cond_29_31 = sum(bits[28:31]) == 0
        cond_33_36 = all(bits[i] for i in (32, 33, 34, 35))  # ESP33..36
        if cond_29_31:
            n_29_31_absent += 1
        if cond_33_36:
            n_33_36_present += 1
        if cond_29_31 and cond_33_36:
            n_ural_like += 1
    print(f"  Strains with spacers 29-31 absent: {n_29_31_absent}/{len(rows)}")
    print(f"  Strains with spacers 33-36 present: {n_33_36_present}/{len(rows)}")
    print(f"  Strains matching BOTH (Ural-compatible): {n_ural_like}/{len(rows)}")

    # Per sub-clade
    for sub in SUB_CLADES:
        sub_rows = [r for r in rows if r["sub_lineage"] == sub]
        n_u = 0
        for r in sub_rows:
            bits = [int(b) for b in r["pattern_binary"]]
            if sum(bits[28:31]) == 0 and all(bits[i] for i in (32, 33, 34, 35)):
                n_u += 1
        print(f"  {sub}: {n_u}/{len(sub_rows)} Ural-compatible")

    # Write summary
    summary_path = RESULTS_DIR / "phase1e_spoligo_summary.txt"
    with open(summary_path, "w") as f:
        f.write("L4.13 -- Phase 1e: in silico spoligotyping\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Total strains typed: {len(rows)}\n")
        f.write(f"Threshold for presence: covered_bases_percent >= {THRESHOLD}\n\n")
        f.write("=== Unique octal patterns (global) ===\n")
        for oct_code, n in all_octs.most_common():
            pct = n / len(rows) * 100
            f.write(f"  {oct_code}  n={n:<4d} ({pct:4.1f}%)\n")
        f.write(f"\n=== Ural signature ===\n")
        f.write(f"Ural-compatible (ESP29-31 absent AND ESP33-36 present): "
                f"{n_ural_like}/{len(rows)} "
                f"({n_ural_like/len(rows)*100:.1f}%)\n")
        for sub in SUB_CLADES:
            sub_rows = [r for r in rows if r["sub_lineage"] == sub]
            n_u = 0
            for r in sub_rows:
                bits = [int(b) for b in r["pattern_binary"]]
                if sum(bits[28:31]) == 0 and all(bits[i] for i in (32, 33, 34, 35)):
                    n_u += 1
            f.write(f"  {sub}: {n_u}/{len(sub_rows)} "
                    f"({n_u/len(sub_rows)*100:.1f}%)\n")
    print(f"\n[write] {summary_path}")

    print("\n" + "=" * 70)
    print("Phase 1e complete.")
    print("=" * 70)


if __name__ == "__main__":
    main()
