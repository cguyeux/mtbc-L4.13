"""
L4.13 Phase 3e -- Stratified McDonald-Kreitman test
====================================================

Follows up on phase 3c (global MK test, p=0.72, neutrality) and
phase 3d (pathway-level enrichment: PDIM OR=6.15 p=9e-6, ESX OR=4.95
p=0.006). Here we formalise the focal selection signal by running
MK tests **stratified by functional pathway**.

Hypothesis: the global MK neutrality hides a focal positive selection
on the PDIM biosynthesis pathway and the T7SS ESX-2/4 accessory modules.
If true, the stratified MK on these pathways should show Dn/Ds >> Pn/Ps
with significant Fisher p-values and strongly positive DoS, while the
rest-of-genome stratum should remain neutral.

Input:
    data/mk_input_annotated.csv (Tier 1 + Tier 4, annotated)

Strata tested (two regexes per pathway: strict = named genes; extended
= broader gene-name prefix):

    PDIM_strict : pks1, pks5, pks10, mas, ppsD, tesA,
                  fadD8, fadD9, fadD21, fadD29, fadD30
    PDIM_extended: starts with pks | mas | ppsD | tesA | fadD (any)
    ESX_strict  : eccB2, eccC2, eccD2, eccC4, espK
    ESX_extended: starts with ecc | esp | esx
    PE_PPE      : starts with PE | PPE (strict)
    Rest        : everything else

Output:
    résultats/phase3e_mk_stratified.txt

Usage:
    python phase3e_mk_stratified.py
"""

from __future__ import annotations

import csv
import re
from collections import Counter
from pathlib import Path

import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
RES = ROOT / "résultats"
INPUT = DATA / "mk_input_annotated.csv"
OUT = RES / "phase3e_mk_stratified.txt"

NS_EFFECTS = {"missense_variant", "stop_gained", "stop_lost"}
S_EFFECTS = {"synonymous_variant"}

PDIM_STRICT = {
    "pks1", "pks5", "pks10", "mas", "ppsD", "tesA",
    "fadD8", "fadD9", "fadD21", "fadD29", "fadD30",
}

ESX_STRICT = {
    "eccB2", "eccC2", "eccD2", "eccC4", "espK",
}


def classify_gene(gene_name: str) -> list[str]:
    """Return list of stratum labels a gene belongs to."""
    g = (gene_name or "").strip()
    labels = []
    if not g:
        return labels

    if g in PDIM_STRICT:
        labels.append("PDIM_strict")
    if g in ESX_STRICT:
        labels.append("ESX_strict")
    # Extended patterns on lowercase start
    gl = g.lower()
    if (gl.startswith("pks") or g == "mas" or gl.startswith("ppsd")
            or gl.startswith("tesa") or gl.startswith("fadd")):
        labels.append("PDIM_extended")
    if gl.startswith("ecc") or gl.startswith("esp") or gl.startswith("esx"):
        labels.append("ESX_extended")
    if re.match(r"^PE[_\d]", g) or re.match(r"^PPE\d", g):
        labels.append("PE_PPE")
    return labels


def dos(dn, ds, pn, ps):
    d = dn / (dn + ds) if (dn + ds) else 0
    p = pn / (pn + ps) if (pn + ps) else 0
    return d - p


def bootstrap_dos_ci(dn, ds, pn, ps, n=10000, seed=42):
    rng = np.random.default_rng(seed)
    dt, pt = dn + ds, pn + ps
    if dt == 0 or pt == 0:
        return (float("nan"), float("nan"))
    d_p = dn / dt
    p_p = pn / pt
    d_boot = rng.binomial(dt, d_p, n) / dt
    p_boot = rng.binomial(pt, p_p, n) / pt
    ci = np.percentile(d_boot - p_boot, [2.5, 97.5])
    return float(ci[0]), float(ci[1])


def mk_stats(dn, ds, pn, ps):
    """Return dict of MK statistics for one stratum."""
    if dn + ds + pn + ps == 0:
        return None
    table = np.array([[dn, ds], [pn, ps]])
    odds, p_two = stats.fisher_exact(table, alternative="two-sided")
    p_greater = stats.fisher_exact(table, alternative="greater").pvalue
    d_ratio = dn / ds if ds else float("inf")
    p_ratio = pn / ps if ps else float("inf")
    d_val = dos(dn, ds, pn, ps)
    ci_lo, ci_hi = bootstrap_dos_ci(dn, ds, pn, ps)
    ni = (pn / ps) / (dn / ds) if ds and ps and dn else float("nan")
    alpha = 1 - ni if ni == ni else float("nan")
    return {
        "Dn": dn, "Ds": ds, "Pn": pn, "Ps": ps,
        "Dn/Ds": d_ratio, "Pn/Ps": p_ratio,
        "odds": odds, "p_two": p_two, "p_greater": p_greater,
        "DoS": d_val, "CI_lo": ci_lo, "CI_hi": ci_hi,
        "NI": ni, "alpha": alpha,
    }


def main():
    print("=" * 70)
    print("L4.13 Phase 3e -- stratified MK test")
    print("=" * 70)

    counts = {}  # stratum -> {tier -> Counter(NS/S)}
    strata = ["PDIM_strict", "PDIM_extended",
              "ESX_strict", "ESX_extended", "PE_PPE",
              "Rest_of_genome"]
    for s in strata:
        counts[s] = {1: Counter(), 4: Counter()}
    counts["All"] = {1: Counter(), 4: Counter()}

    gene_hits_tier1 = {s: set() for s in strata[:-1]}

    # Track which genes are assigned to which stratum across all input
    # (to list them in the report)
    with INPUT.open() as f:
        for row in csv.DictReader(f):
            t = int(row["Tier"])
            if t not in (1, 4):
                continue
            eff = row["Effect"]
            is_ns = eff in NS_EFFECTS
            is_s = eff in S_EFFECTS
            if not (is_ns or is_s):
                continue
            gene = row.get("Gene_name", "") or row.get("Locus_tag", "")
            labels = classify_gene(gene)

            # Record in each applicable stratum
            recorded_any = False
            for lab in labels:
                if is_ns:
                    counts[lab][t]["NS"] += 1
                else:
                    counts[lab][t]["S"] += 1
                if t == 1 and lab in gene_hits_tier1:
                    gene_hits_tier1[lab].add(gene)
                recorded_any = True

            # Rest of genome: everything not in PDIM_extended, ESX_extended,
            # PE_PPE (these three are the non-overlapping "focal" sets)
            in_focal = any(x in labels for x in
                           ("PDIM_extended", "ESX_extended", "PE_PPE"))
            if not in_focal:
                if is_ns:
                    counts["Rest_of_genome"][t]["NS"] += 1
                else:
                    counts["Rest_of_genome"][t]["S"] += 1

            if is_ns:
                counts["All"][t]["NS"] += 1
            else:
                counts["All"][t]["S"] += 1

    # Report
    lines = []
    lines.append("L4.13 -- Stratified MK test by functional pathway\n")
    lines.append("=" * 60 + "\n\n")
    lines.append("Strata tested (NS = missense/stop_gained/stop_lost, "
                 "S = synonymous_variant):\n\n")

    header = (f"{'Stratum':<18s} {'Dn':>4s} {'Ds':>4s}  {'Pn':>5s} {'Ps':>5s}  "
              f"{'Dn/Ds':>6s} {'Pn/Ps':>6s}  "
              f"{'p_Fisher':>10s} {'p_greater':>10s}  "
              f"{'DoS':>7s} {'CI95':>18s}\n")
    lines.append(header)
    lines.append("-" * len(header) + "\n")

    for stratum in ["All"] + strata:
        c1 = counts[stratum][1]
        c4 = counts[stratum][4]
        st = mk_stats(c1["NS"], c1["S"], c4["NS"], c4["S"])
        if st is None:
            continue
        ci = f"[{st['CI_lo']:+.3f},{st['CI_hi']:+.3f}]"
        lines.append(
            f"{stratum:<18s} {st['Dn']:>4d} {st['Ds']:>4d}  "
            f"{st['Pn']:>5d} {st['Ps']:>5d}  "
            f"{st['Dn/Ds']:>6.2f} {st['Pn/Ps']:>6.2f}  "
            f"{st['p_two']:>10.4g} {st['p_greater']:>10.4g}  "
            f"{st['DoS']:>+7.3f} {ci:>18s}\n"
        )

    lines.append("\n")
    lines.append("Notes:\n")
    lines.append("  - p_Fisher = two-sided Fisher exact test on 2x2 table\n")
    lines.append("  - p_greater = one-sided test for NS excess in divergence\n")
    lines.append("  - DoS = Direction of Selection = Dn/(Dn+Ds) - Pn/(Pn+Ps)\n")
    lines.append("  - CI95 = 95% bootstrap CI on DoS (10,000 replicates)\n")
    lines.append("  - Rest_of_genome excludes PDIM_extended, ESX_extended, PE_PPE\n\n")

    lines.append("=" * 60 + "\n")
    lines.append("Genes hit in Tier 1 (core-exclusive L4.13) per stratum\n")
    lines.append("=" * 60 + "\n")
    for stratum, genes in gene_hits_tier1.items():
        lines.append(f"\n{stratum} ({len(genes)} unique genes in Tier 1):\n")
        for g in sorted(genes):
            lines.append(f"  {g}\n")

    OUT.write_text("".join(lines))
    print("".join(lines))
    print(f"\n[write] {OUT}")


if __name__ == "__main__":
    main()
