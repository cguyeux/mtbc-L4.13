#!/usr/bin/env python3
"""
L4.13 Analysis -- Phase 1d: Validation of candidate barcoding markers
======================================================================

Validates three candidate SPDI markers proposed for the L4.13 clade and
its two sub-clades, after reassignment of the 1/2 convention:

    L4.13     -- NC_000962.3:3222055:G:A    (new parent marker)
    L4.13.1   -- NC_000962.3:187910:T:C     (majority sub-clade, 197 strains)
    L4.13.2   -- NC_000962.3:1129237:G:A    (minority sub-clade, 61 strains)

For each marker, two tests are run:

1. INCLUSION test -- the marker must be present in (almost) all strains
   of the target clade. We report penetrance per sub-clade.

2. EXCLUSION test -- the marker must be (almost) absent in sister
   sub-lineages L4.1.*, L4.2.*, L4.3.*, L4.4.*, L4.6.*, L4.8, L4.9,
   L4.14, L4.15. We sample up to N strains per sub-lineage.

Outputs:
    résultats/phase1d_marker_validation.txt
    résultats/phase1d_marker_validation.tsv

Usage:
    python phase1d_validate_markers.py [--sample-size 100]
"""

import argparse
import random
import sys
from collections import Counter, defaultdict
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent
REPO_ROOT = PROJECT_DIR.parent

BDD = Path(__file__).resolve().parent.parent / "external" / "spdi_db"  # per-strain SPDI profiles; see external/README.md
RESULTS_DIR = PROJECT_DIR / "résultats"
REFERENCE = "NC_000962.3"

CANDIDATES = [
    ("L4.13",   "NC_000962.3:3222055:G:A"),
    ("L4.13.1", "NC_000962.3:187910:T:C"),
    ("L4.13.2", "NC_000962.3:1129237:G:A"),
]

# Target clade (after convention flip)
TARGET_SUBCLADES = {
    # physical dir -> new label after convention flip
    "L4.13.1": "L4.13.2",  # 61 strains -> L4.13.2 (minority)
    "L4.13.2": "L4.13.1",  # 197 strains -> L4.13.1 (majority)
}

# Sister lineages -- all L4.x that are NOT in L4.13
SISTER_LINEAGES = [
    "L4.1", "L4.1.1", "L4.1.2", "L4.1_proto",
    "L4.2.1", "L4.2.2",
    "L4.3.1", "L4.3.2", "L4.3.3", "L4.3.4",
    "L4.4.2",
    "L4.6.1", "L4.6.2",
    "L4.8", "L4.9",
    "L4.14", "L4.15",
]


def list_strains(sub):
    d = BDD / sub
    if not d.is_dir():
        return []
    return [p.name for p in d.iterdir()
            if p.is_dir() and (p / REFERENCE / "spdi.txt").is_file()]


def spdi_set(sub, sra):
    f = BDD / sub / sra / REFERENCE / "spdi.txt"
    with open(f) as fh:
        return {line.strip() for line in fh if line.strip()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample-size", type=int, default=100,
                    help="Max strains sampled per sister sub-lineage")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    random.seed(args.seed)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("L4.13 Analysis -- Phase 1d: Candidate marker validation")
    print("=" * 70)

    # ---- INCLUSION TEST ----
    print("\n[inclusion] testing 3 candidate markers on 258 L4.13 strains")
    inclusion_counts = {lab: Counter() for lab, _ in CANDIDATES}
    inclusion_totals = Counter()
    inclusion_misses = defaultdict(list)  # (candidate_label, physical_dir) -> [strains missing]

    for physical_dir, new_label in TARGET_SUBCLADES.items():
        strains = list_strains(physical_dir)
        inclusion_totals[physical_dir] = len(strains)
        for sra in strains:
            sample = spdi_set(physical_dir, sra)
            for lab, spdi in CANDIDATES:
                if spdi in sample:
                    inclusion_counts[lab][physical_dir] += 1
                else:
                    inclusion_misses[(lab, physical_dir)].append(sra)

    # ---- EXCLUSION TEST ----
    print(f"[exclusion] sampling up to {args.sample_size} strains per "
          f"sister sub-lineage ({len(SISTER_LINEAGES)} sub-lineages)")
    exclusion_counts = {lab: defaultdict(lambda: [0, 0])  # [n_present, n_sampled]
                        for lab, _ in CANDIDATES}

    for sub in SISTER_LINEAGES:
        strains = list_strains(sub)
        if not strains:
            continue
        if len(strains) > args.sample_size:
            strains = random.sample(strains, args.sample_size)
        for sra in strains:
            sample = spdi_set(sub, sra)
            for lab, spdi in CANDIDATES:
                exclusion_counts[lab][sub][1] += 1
                if spdi in sample:
                    exclusion_counts[lab][sub][0] += 1

    # ---- REPORT ----
    lines = []
    lines.append("L4.13 -- Phase 1d: Candidate marker validation\n")
    lines.append("=" * 60 + "\n\n")
    lines.append(f"Random seed: {args.seed}, sample size per sister: "
                 f"{args.sample_size}\n\n")

    lines.append("CONVENTION NOTE\n")
    lines.append("-" * 60 + "\n")
    lines.append("  After convention flip, the 197-strain sub-clade\n")
    lines.append("  becomes the new majority sub-clade L4.13.1; and\n")
    lines.append("  the 61-strain sub-clade becomes L4.13.2.\n\n")

    # Inclusion summary
    lines.append("INCLUSION TEST -- marker penetrance in L4.13\n")
    lines.append("-" * 60 + "\n")
    for lab, spdi in CANDIDATES:
        lines.append(f"\n  Candidate: {lab} = {spdi}\n")
        for physical_dir in ["L4.13.1", "L4.13.2"]:
            n_tot = inclusion_totals[physical_dir]
            n_pres = inclusion_counts[lab][physical_dir]
            pct = n_pres / n_tot * 100.0 if n_tot else 0.0
            new = TARGET_SUBCLADES[physical_dir]
            lines.append(f"    physical {physical_dir} (new {new}, n={n_tot}): "
                         f"present in {n_pres}/{n_tot} ({pct:.1f}%)\n")
        # All 258
        total_all = sum(inclusion_totals.values())
        present_all = sum(inclusion_counts[lab].values())
        lines.append(f"    all 258 L4.13 strains: "
                     f"{present_all}/{total_all} ({present_all/total_all*100:.1f}%)\n")
        # Missing examples
        for (l, d), misses in inclusion_misses.items():
            if l == lab and misses:
                shown = ", ".join(misses[:5])
                extra = f" (+{len(misses)-5} more)" if len(misses) > 5 else ""
                lines.append(f"    absent in {d}: {shown}{extra}\n")

    # Exclusion summary
    lines.append("\n\nEXCLUSION TEST -- marker in sister L4 sub-lineages\n")
    lines.append("-" * 60 + "\n")
    for lab, spdi in CANDIDATES:
        lines.append(f"\n  Candidate: {lab} = {spdi}\n")
        total_pres = 0
        total_samp = 0
        for sub in SISTER_LINEAGES:
            n_pres, n_samp = exclusion_counts[lab][sub]
            if n_samp == 0:
                continue
            total_pres += n_pres
            total_samp += n_samp
            pct = n_pres / n_samp * 100.0
            flag = "  <-- LEAK" if n_pres > 0 else ""
            lines.append(f"    {sub:<15s} {n_pres}/{n_samp} ({pct:.1f}%){flag}\n")
        if total_samp:
            pct = total_pres / total_samp * 100.0
            lines.append(f"    TOTAL sister lineages: {total_pres}/{total_samp} ({pct:.2f}%)\n")

    # Verdict
    lines.append("\n\nVERDICT\n")
    lines.append("-" * 60 + "\n")
    for lab, spdi in CANDIDATES:
        in_all = sum(inclusion_counts[lab].values())
        tot_in = sum(inclusion_totals.values())
        out_pres = sum(exclusion_counts[lab][s][0] for s in SISTER_LINEAGES)
        out_samp = sum(exclusion_counts[lab][s][1] for s in SISTER_LINEAGES)

        sensitivity = in_all / tot_in * 100.0 if tot_in else 0.0
        leak_rate = out_pres / out_samp * 100.0 if out_samp else 0.0
        spec = 100.0 - leak_rate

        if lab == "L4.13":
            expected = tot_in  # should be in all 258
        elif lab == "L4.13.1":
            expected = inclusion_totals["L4.13.2"]  # majority = physical L4.13.2
        else:
            expected = inclusion_totals["L4.13.1"]  # minority = physical L4.13.1

        lines.append(f"\n  {lab} ({spdi}):\n")
        lines.append(f"    inclusion (all L4.13): {in_all}/{tot_in} "
                     f"({sensitivity:.1f}%)\n")
        lines.append(f"    exclusion (sister L4): {out_pres}/{out_samp} "
                     f"leaking ({leak_rate:.2f}%), spec={spec:.2f}%\n")
        status = "VALID" if (sensitivity >= 95.0 and leak_rate <= 1.0) else "REVIEW"
        lines.append(f"    status: {status}\n")

    report = "".join(lines)
    print("\n" + report)

    txt_path = RESULTS_DIR / "phase1d_marker_validation.txt"
    with open(txt_path, "w") as f:
        f.write(report)
    print(f"\n[write] {txt_path}")

    # TSV report
    tsv_path = RESULTS_DIR / "phase1d_marker_validation.tsv"
    with open(tsv_path, "w") as f:
        f.write("candidate_label\tspdi\tscope\tsub_lineage\tn_present\tn_total\tpct\n")
        for lab, spdi in CANDIDATES:
            for physical_dir in ["L4.13.1", "L4.13.2"]:
                n_tot = inclusion_totals[physical_dir]
                n_pres = inclusion_counts[lab][physical_dir]
                pct = n_pres / n_tot * 100.0 if n_tot else 0.0
                f.write(f"{lab}\t{spdi}\tinclusion\t{physical_dir}\t{n_pres}\t{n_tot}\t{pct:.2f}\n")
            for sub in SISTER_LINEAGES:
                n_pres, n_samp = exclusion_counts[lab][sub]
                if n_samp == 0:
                    continue
                pct = n_pres / n_samp * 100.0
                f.write(f"{lab}\t{spdi}\texclusion\t{sub}\t{n_pres}\t{n_samp}\t{pct:.2f}\n")
    print(f"[write] {tsv_path}")

    print("\n" + "=" * 70)
    print("Phase 1d complete.")
    print("=" * 70)


if __name__ == "__main__":
    main()
