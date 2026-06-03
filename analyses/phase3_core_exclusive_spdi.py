"""
L4.13 Analysis -- Phase 3: Core-exclusive SPDIs per clade
==========================================================

For each of L4.13 (258), L4.13.1 (61), L4.13.2 (197), extract:
    - Core SPDIs: present in >= MIN_FREQ (default 95%) of the target clade
    - Exclusive: absent from all tested sister sub-lineages (up to N sampled
      strains per sister, default 50)

Sister sub-lineages tested (up to SAMPLE_SIZE strains each, seeded):
    L4.1, L4.1.1, L4.1.2, L4.1_proto, L4.2.1, L4.2.2,
    L4.3.1, L4.3.2, L4.3.3, L4.3.4, L4.4.2,
    L4.6.1, L4.6.2, L4.8, L4.9, L4.14, L4.15

Also, inside L4.13:
    - Core-exclusive of L4.13.1 : present >=95% in L4.13.1, absent in L4.13.2
    - Core-exclusive of L4.13.2 : present >=95% in L4.13.2, absent in L4.13.1

Outputs (per target):
    data/core_exclusive_L4.13.txt     (one SPDI per line)
    data/core_exclusive_L4.13.1.txt
    data/core_exclusive_L4.13.2.txt
    résultats/phase3_summary.txt

Usage:
    python phase3_core_exclusive_spdi.py [--min-freq 0.95] [--sample 50]
"""

from __future__ import annotations

import argparse
import random
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BDD = Path(__file__).resolve().parent.parent / "external" / "spdi_db"  # per-strain SPDI profiles; see external/README.md
DATA_DIR = ROOT / "data"
RES_DIR = ROOT / "résultats"
REFERENCE = "NC_000962.3"

TARGET_SUBCLADES = ["L4.13.1", "L4.13.2"]

SISTER_LINEAGES = [
    "L4.1", "L4.1.1", "L4.1.2", "L4.1_proto",
    "L4.2.1", "L4.2.2",
    "L4.3.1", "L4.3.2", "L4.3.3", "L4.3.4",
    "L4.4.2",
    "L4.6.1", "L4.6.2",
    "L4.8", "L4.9",
    "L4.14", "L4.15",
]


def list_strains(d: Path) -> list[str]:
    if not d.is_dir():
        return []
    return sorted(p.name for p in d.iterdir()
                  if p.is_dir() and (p / REFERENCE / "spdi.txt").is_file())


def load_spdis(d: Path, sra: str) -> set[str]:
    with (d / sra / REFERENCE / "spdi.txt").open() as f:
        return {line.strip() for line in f if line.strip()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-freq", type=float, default=0.95,
                    help="Minimum frequency in target clade (default 0.95)")
    ap.add_argument("--sample", type=int, default=50,
                    help="Max strains sampled per sister sub-lineage")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    random.seed(args.seed)
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    RES_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("L4.13 Phase 3 -- core-exclusive SPDIs")
    print("=" * 70)

    # Load L4.13 strains (ingroup)
    ingroup = {}  # label -> set
    for sub in TARGET_SUBCLADES:
        d = BDD / sub
        for sra in list_strains(d):
            ingroup[f"{sub}_{sra}"] = (sub, load_spdis(d, sra))
    print(f"\n[load] ingroup: {len(ingroup)} strains")

    # Load sister strains
    sister_spdis = []  # list of (label, sub, set)
    for sub in SISTER_LINEAGES:
        d = BDD / sub
        strains = list_strains(d)
        if not strains:
            continue
        if len(strains) > args.sample:
            strains = random.sample(strains, args.sample)
        for sra in strains:
            sister_spdis.append((f"{sub}_{sra}", sub, load_spdis(d, sra)))
    print(f"[load] sister: {len(sister_spdis)} strains across "
          f"{len(SISTER_LINEAGES)} sub-lineages")

    # Pan-SPDI (union)
    all_ingroup_spdis = set()
    for _, s in ingroup.values():
        all_ingroup_spdis |= s
    print(f"[pan-SPDI ingroup] {len(all_ingroup_spdis)} SPDIs")

    # Frequency per SPDI in each target
    def freq(spdis_by_label, target_sub=None):
        """Return dict spdi -> count in subset."""
        c = {}
        for _, (sub, s) in spdis_by_label.items():
            if target_sub is not None and sub != target_sub:
                continue
            for spdi in s:
                c[spdi] = c.get(spdi, 0) + 1
        return c

    n_l413 = len(ingroup)
    n_l131 = sum(1 for _, (s, _) in ingroup.items() if s == "L4.13.1")
    n_l132 = sum(1 for _, (s, _) in ingroup.items() if s == "L4.13.2")
    print(f"\n  L4.13={n_l413}, L4.13.1={n_l131}, L4.13.2={n_l132}")

    # Counts
    f_l413 = freq(ingroup, None)
    f_l131 = freq(ingroup, "L4.13.1")
    f_l132 = freq(ingroup, "L4.13.2")

    # Sister presence count (any sister strain)
    sister_contains = {}  # spdi -> count
    for _, _, s in sister_spdis:
        for spdi in s:
            sister_contains[spdi] = sister_contains.get(spdi, 0) + 1

    n_sisters = len(sister_spdis)

    # Core-exclusive for L4.13 (whole clade)
    thr_l413 = int(args.min_freq * n_l413)
    core_excl_l413 = sorted(
        spdi for spdi, c in f_l413.items()
        if c >= thr_l413 and sister_contains.get(spdi, 0) == 0
    )

    # Core-exclusive for L4.13.1 (present in most L4.13.1, absent in L4.13.2 + absent in sisters)
    thr_l131 = int(args.min_freq * n_l131)
    core_excl_l131 = sorted(
        spdi for spdi, c in f_l131.items()
        if c >= thr_l131
        and f_l132.get(spdi, 0) == 0
        and sister_contains.get(spdi, 0) == 0
    )

    # Core-exclusive for L4.13.2
    thr_l132 = int(args.min_freq * n_l132)
    core_excl_l132 = sorted(
        spdi for spdi, c in f_l132.items()
        if c >= thr_l132
        and f_l131.get(spdi, 0) == 0
        and sister_contains.get(spdi, 0) == 0
    )

    # Write SPDIs
    (DATA_DIR / "core_exclusive_L4.13.txt").write_text(
        "\n".join(core_excl_l413) + "\n")
    (DATA_DIR / "core_exclusive_L4.13.1.txt").write_text(
        "\n".join(core_excl_l131) + "\n")
    (DATA_DIR / "core_exclusive_L4.13.2.txt").write_text(
        "\n".join(core_excl_l132) + "\n")

    # Summary
    print(f"\n[result]")
    print(f"  L4.13 core-exclusive: {len(core_excl_l413)} SPDIs")
    print(f"  L4.13.1 core-exclusive: {len(core_excl_l131)} SPDIs")
    print(f"  L4.13.2 core-exclusive: {len(core_excl_l132)} SPDIs")

    # Write report
    lines = []
    lines.append("L4.13 -- Phase 3: core-exclusive SPDIs\n")
    lines.append("=" * 50 + "\n\n")
    lines.append(f"Parameters:\n")
    lines.append(f"  min_freq in target clade: {args.min_freq}\n")
    lines.append(f"  sample size per sister sub-lineage: {args.sample}\n")
    lines.append(f"  seed: {args.seed}\n\n")
    lines.append(f"Ingroup: L4.13 = {n_l413} (L4.13.1={n_l131}, L4.13.2={n_l132})\n")
    lines.append(f"Sister set: {n_sisters} strains from {len(SISTER_LINEAGES)} sub-lineages\n\n")
    lines.append(f"=== Core-exclusive SPDIs ===\n")
    lines.append(f"  L4.13   (whole clade): {len(core_excl_l413)}\n")
    lines.append(f"  L4.13.1 (vs L4.13.2 + sisters): {len(core_excl_l131)}\n")
    lines.append(f"  L4.13.2 (vs L4.13.1 + sisters): {len(core_excl_l132)}\n\n")

    lines.append(f"=== Top-15 L4.13 core-exclusive SPDIs (lowest position) ===\n")
    for spdi in core_excl_l413[:15]:
        lines.append(f"  {spdi}\n")

    if len(core_excl_l413) > 15:
        lines.append(f"  ... (+{len(core_excl_l413)-15} more)\n")

    lines.append(f"\n=== L4.13.1 core-exclusive (all {len(core_excl_l131)}) ===\n")
    for spdi in core_excl_l131:
        lines.append(f"  {spdi}\n")

    lines.append(f"\n=== L4.13.2 core-exclusive (all {len(core_excl_l132)}) ===\n")
    for spdi in core_excl_l132:
        lines.append(f"  {spdi}\n")

    summary_path = RES_DIR / "phase3_summary.txt"
    summary_path.write_text("".join(lines))
    print(f"\n[write] {summary_path}")
    print(f"[write] {DATA_DIR / 'core_exclusive_L4.13.txt'}")
    print(f"[write] {DATA_DIR / 'core_exclusive_L4.13.1.txt'}")
    print(f"[write] {DATA_DIR / 'core_exclusive_L4.13.2.txt'}")

    print("\n" + "=" * 70)
    print("Phase 3 complete.")
    print("=" * 70)


if __name__ == "__main__":
    main()
