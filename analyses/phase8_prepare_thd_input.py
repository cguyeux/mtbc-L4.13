"""
L4.13 Phase 8 -- Prepare THD haplotype matrix input
====================================================

Builds a binary haplotype matrix (rows = 258 L4.13 strains, columns =
SPDIs polymorphic within L4.13) for Rasigade's Time-scaled Haplotypic
Density (THD) analysis.

Output:
    data/thd_haplotypes.csv -- id, sub_lineage, S1, S2, ..., Sm
         where S_i are 0/1 markers at polymorphic SPDI positions.

Only intra-L4.13 polymorphic sites are used (present in 2..N-2 of 258
strains). These are the same sites used for the MK Tier 4 (phase 3b).
"""

import csv
from pathlib import Path
from collections import Counter

ROOT = Path(__file__).resolve().parent.parent
BDD = Path(__file__).resolve().parent.parent / "external" / "spdi_db"  # per-strain SPDI profiles; see external/README.md
DATA = ROOT / "data"
REFERENCE = "NC_000962.3"
SUB_CLADES = ["L4.13.1", "L4.13.2"]

OUT = DATA / "thd_haplotypes.csv"


def list_strains(d):
    if not d.is_dir(): return []
    return sorted(p.name for p in d.iterdir()
                  if p.is_dir() and (p / REFERENCE / "spdi.txt").is_file())


def load_spdis(d, sra):
    with (d / sra / REFERENCE / "spdi.txt").open() as f:
        return {line.strip() for line in f if line.strip()}


def main():
    # Load strains
    strains = []  # (id, group, set)
    for sub in SUB_CLADES:
        for sra in list_strains(BDD / sub):
            strains.append((sra, sub, load_spdis(BDD / sub, sra)))
    n = len(strains)
    print(f"[load] {n} L4.13 strains")

    # Union of SPDIs
    all_spdis = set()
    for _, _, s in strains:
        all_spdis |= s

    # Count frequency
    counts = {s: 0 for s in all_spdis}
    for _, _, sset in strains:
        for s in sset:
            counts[s] += 1

    # Keep polymorphic intra-L4.13
    poly = sorted(s for s, c in counts.items() if 2 <= c <= n - 2)
    print(f"[filter] {len(poly)} polymorphic sites (2..{n-2}/258)")

    # Build wide CSV: id, sub_lineage, S1..Sm
    # Column names: use index (S1, S2, ...) to avoid huge SPDI headers
    col_names = [f"S{i+1}" for i in range(len(poly))]
    with OUT.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["id", "sub_lineage"] + col_names)
        for sid, grp, sset in strains:
            row = [sid, grp] + [1 if s in sset else 0 for s in poly]
            w.writerow(row)
    print(f"[write] {OUT} ({n} rows x {len(poly)} markers)")


if __name__ == "__main__":
    main()
