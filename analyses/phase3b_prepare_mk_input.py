"""
L4.13 Phase 3b -- Prepare tier-annotated CSV for the MK test
=============================================================

Builds the tier-annotated CSV expected by `/mk-ascertainment`:

  Tier 1 -- Divergence (core-exclusive to L4.13): 187 SPDIs
            present in >= 95% of L4.13, absent in sister L4.x.
  Tier 4 -- Polymorphism (intra-L4.13, no cross-lineage filter):
            present in 2..N-2 of the 258 L4.13 strains.
            NOTE: the skill explicitly asks for Tier 4 WITHOUT
            cross-lineage filter, so this mimics the true polymorphism
            pool used to control ascertainment bias in the MK test.

Annotates Tier 4 via the spdi-annotation script and merges with the
already-annotated core-exclusive Tier 1.

Output:
    data/mk_input_annotated.csv -- columns SPDI, Effect, Tier, ...

Usage:
    python phase3b_prepare_mk_input.py
"""

from __future__ import annotations

import csv
import subprocess
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BDD = Path(__file__).resolve().parent.parent / "external" / "spdi_db"  # per-strain SPDI profiles; see external/README.md
DATA = ROOT / "data"
REFERENCE = "NC_000962.3"
SUB_CLADES = ["L4.13.1", "L4.13.2"]

SKILL_ANNOTATE = Path(__file__).resolve().parent.parent / "external" / "annotate_spdis.py"
GFF3 = Path(__file__).resolve().parent.parent / "external" / "NC_000962.3.gff3"
GB = Path(__file__).resolve().parent.parent / "external" / "NC_000962.3.gb"


def list_strains(d: Path) -> list[str]:
    if not d.is_dir():
        return []
    return sorted(p.name for p in d.iterdir()
                  if p.is_dir() and (p / REFERENCE / "spdi.txt").is_file())


def load_spdis(d: Path, sra: str) -> set[str]:
    with (d / sra / REFERENCE / "spdi.txt").open() as f:
        return {line.strip() for line in f if line.strip()}


def main():
    print("=" * 70)
    print("L4.13 Phase 3b -- Prepare MK test input")
    print("=" * 70)

    # Load all L4.13 strains
    all_spdis_by_strain = {}
    for sub in SUB_CLADES:
        d = BDD / sub
        for sra in list_strains(d):
            all_spdis_by_strain[f"{sub}_{sra}"] = load_spdis(d, sra)
    n_total = len(all_spdis_by_strain)
    print(f"\n[load] {n_total} L4.13 strains")

    # Build Tier 4: polymorphic within L4.13 (2..N-2), no cross-lineage filter
    counts = Counter()
    for spdis in all_spdis_by_strain.values():
        for s in spdis:
            counts[s] += 1
    tier4_spdis = [s for s, c in counts.items() if 2 <= c <= n_total - 2]
    print(f"[tier 4] polymorphic intra-L4.13 (2..{n_total-2}): "
          f"{len(tier4_spdis)} SPDIs")

    # Load existing Tier 1 (the 187 core-exclusive L4.13)
    tier1_file = DATA / "core_exclusive_L4.13.txt"
    tier1_spdis = [line.strip() for line in tier1_file.read_text().splitlines()
                   if line.strip()]
    print(f"[tier 1] core-exclusive L4.13: {len(tier1_spdis)} SPDIs")

    # Drop Tier 1 SPDIs from Tier 4 (no double-counting)
    tier1_set = set(tier1_spdis)
    tier4_clean = [s for s in tier4_spdis if s not in tier1_set]
    print(f"[tier 4 - tier 1] after excluding Tier 1: {len(tier4_clean)} SPDIs")

    # Write combined input (to annotate via spdi-annotation script)
    combined_input = DATA / "mk_input_raw.csv"
    with combined_input.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["SPDI", "Tier"])
        for spdi in tier1_spdis:
            w.writerow([spdi, 1])
        for spdi in tier4_clean:
            w.writerow([spdi, 4])
    print(f"[write] {combined_input} ({len(tier1_spdis)+len(tier4_clean)} SPDIs)")

    # Annotate via the skill script
    print(f"\n[annotate] running local annotation on combined set...")
    combined_out = DATA / "mk_input_annotated.csv"
    cmd = [
        sys.executable, str(SKILL_ANNOTATE),
        str(combined_input),
        "--gff3", str(GFF3),
        "--genbank", str(GB),
        "--skip-tbannotator",
        "-o", str(combined_out),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    print(proc.stdout[-500:] if proc.stdout else "")
    if proc.returncode != 0:
        print(proc.stderr[-500:])
        sys.exit(1)
    print(f"[write] {combined_out}")

    # Quick counts
    tier_ns = Counter()
    tier_s = Counter()
    with combined_out.open() as f:
        for row in csv.DictReader(f):
            t = int(row["Tier"])
            eff = row["Effect"]
            if eff in ("missense_variant", "stop_gained", "stop_lost"):
                tier_ns[t] += 1
            elif eff == "synonymous_variant":
                tier_s[t] += 1
    print("\n[counts] (NS / S)")
    for t in sorted(set(tier_ns) | set(tier_s)):
        print(f"  Tier {t}: NS={tier_ns[t]:4d}  S={tier_s[t]:4d}  "
              f"ratio={tier_ns[t]/tier_s[t] if tier_s[t] else float('inf'):.2f}")

    print("\n" + "=" * 70)
    print("Phase 3b complete. Ready for MK test.")
    print("=" * 70)


if __name__ == "__main__":
    main()
