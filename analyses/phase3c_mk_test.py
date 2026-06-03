"""
L4.13 Phase 3c -- McDonald-Kreitman test
=========================================

Tests whether the excess of non-synonymous variants observed among the
core-exclusive SPDIs of L4.13 (Tier 1) is statistically distinguishable
from the intra-lineage polymorphism pool (Tier 4, no cross-lineage
filter). If the test is significant with DoS > 0, the excess cannot be
explained by the ascertainment bias of the core-exclusive filter alone.

Input: data/mk_input_annotated.csv (from phase 3b)

Reports:
    - 2x2 contingency table (Dn, Ds, Pn, Ps)
    - Fisher's exact test p-value
    - Direction of Selection (DoS) = Dn/(Dn+Ds) - Pn/(Pn+Ps)
    - Bootstrap 95% CI on DoS (10 000 replicates)
    - alpha = 1 - (Ds * Pn) / (Dn * Ps)
    - Binomial test: p(Dn observed | p_NS = 2/3 under neutral CDS)

Usage:
    python phase3c_mk_test.py
"""

import csv
from collections import Counter
from pathlib import Path

import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
RESULTS = ROOT / "résultats"
INPUT = DATA / "mk_input_annotated.csv"


def load_counts(path):
    ns = {'missense_variant', 'stop_gained', 'stop_lost'}
    s = {'synonymous_variant'}
    by_tier_ns, by_tier_s = Counter(), Counter()
    ns_spdis = {1: [], 4: []}
    s_spdis = {1: [], 4: []}
    with path.open() as f:
        for row in csv.DictReader(f):
            t = int(row['Tier'])
            if t not in (1, 4):
                continue
            eff = row['Effect']
            if eff in ns:
                by_tier_ns[t] += 1
                ns_spdis[t].append(row['SPDI'])
            elif eff in s:
                by_tier_s[t] += 1
                s_spdis[t].append(row['SPDI'])
    return by_tier_ns, by_tier_s, ns_spdis, s_spdis


def dos(dn, ds, pn, ps):
    """Direction of Selection (Stoletzki & Eyre-Walker 2011)."""
    d_frac = dn / (dn + ds) if (dn + ds) > 0 else 0
    p_frac = pn / (pn + ps) if (pn + ps) > 0 else 0
    return d_frac - p_frac


def bootstrap_dos(dn, ds, pn, ps, n_boot=10000, seed=42):
    rng = np.random.default_rng(seed)
    d_total = dn + ds
    p_total = pn + ps
    d_p = dn / d_total if d_total else 0
    p_p = pn / p_total if p_total else 0
    d_boot = rng.binomial(d_total, d_p, n_boot) / d_total if d_total else np.zeros(n_boot)
    p_boot = rng.binomial(p_total, p_p, n_boot) / p_total if p_total else np.zeros(n_boot)
    dos_boot = d_boot - p_boot
    return np.percentile(dos_boot, [2.5, 97.5])


def main():
    print("=" * 70)
    print("L4.13 Phase 3c -- McDonald-Kreitman test")
    print("=" * 70)

    by_tier_ns, by_tier_s, _, _ = load_counts(INPUT)

    Dn = by_tier_ns[1]
    Ds = by_tier_s[1]
    Pn = by_tier_ns[4]
    Ps = by_tier_s[4]

    print(f"\n[contingency]")
    print(f"                NS    S")
    print(f"  Divergence  {Dn:4d}  {Ds:3d}   (Tier 1: core-exclusive L4.13)")
    print(f"  Polymorph.  {Pn:4d}  {Ps:3d}   (Tier 4: intra-L4.13, no x-lineage)")

    table = np.array([[Dn, Ds], [Pn, Ps]])
    odds, p_fisher = stats.fisher_exact(table, alternative='two-sided')
    p_greater = stats.fisher_exact(table, alternative='greater').pvalue

    d_val = dos(Dn, Ds, Pn, Ps)
    ci_lo, ci_hi = bootstrap_dos(Dn, Ds, Pn, Ps)

    # Neutrality Index = (Pn/Ps) / (Dn/Ds)
    ni = (Pn / Ps) / (Dn / Ds) if Ds > 0 and Ps > 0 else float('nan')
    alpha = 1 - ni if ni == ni else float('nan')

    d_ratio = Dn / Ds if Ds > 0 else float('inf')
    p_ratio = Pn / Ps if Ps > 0 else float('inf')

    print(f"\n[ratios]")
    print(f"  Dn/Ds = {d_ratio:.3f}  (divergence, Tier 1)")
    print(f"  Pn/Ps = {p_ratio:.3f}  (polymorphism, Tier 4)")

    print(f"\n[Fisher's exact test]")
    print(f"  odds ratio = {odds:.3f}")
    print(f"  p-value (two-sided) = {p_fisher:.4g}")
    print(f"  p-value (greater Dn excess) = {p_greater:.4g}")

    print(f"\n[Direction of Selection]")
    print(f"  DoS = {d_val:+.4f}")
    print(f"  95% bootstrap CI = [{ci_lo:+.4f}, {ci_hi:+.4f}]")

    print(f"\n[Neutrality Index / alpha]")
    print(f"  NI = (Pn/Ps) / (Dn/Ds) = {ni:.3f}")
    print(f"  alpha = 1 - NI = {alpha:+.4f}")

    # Interpretation
    print(f"\n[interpretation]")
    if p_fisher < 0.05 and d_val > 0:
        print(f"  => Significant (p={p_fisher:.3g}) excess of non-synonymous")
        print(f"     divergence vs polymorphism ratio.")
        print(f"  => DoS > 0 -> consistent with positive selection on the")
        print(f"     L4.13-fixed lineage-defining substitutions.")
    elif p_fisher < 0.05 and d_val < 0:
        print(f"  => Significant (p={p_fisher:.3g}) NEGATIVE DoS -> purifying")
        print(f"     selection dominant or ascertainment artifact.")
    else:
        print(f"  => NOT significant (p={p_fisher:.3g}).")
        print(f"  => The NS excess observed among core-exclusive SPDIs is")
        print(f"     compatible with the intra-lineage polymorphism ratio;")
        print(f"     ascertainment bias cannot be ruled out as the sole")
        print(f"     explanation without further analysis.")

    # Binomial: is Dn compatible with the null proportion p=0.6613 (Nei-Gojobori
    # typical for CDS) given Dn + Ds total?
    # Approximate with standard CDS NS:S ratio ~ 1.95 -> p_NS ~ 0.6613.
    total_cds = Dn + Ds
    if total_cds > 0:
        p_null = 0.6613
        binom = stats.binomtest(Dn, total_cds, p_null, alternative='greater')
        print(f"\n[binomial null (p_NS = {p_null:.4f}, typical CDS)]")
        print(f"  P(NS >= {Dn} | {total_cds}, p={p_null}) = {binom.pvalue:.4g}")

    # Write the result text
    RESULTS.mkdir(exist_ok=True)
    out = RESULTS / "phase3c_mk_test.txt"
    with out.open("w") as f:
        f.write(f"L4.13 -- McDonald-Kreitman test\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Contingency table:\n")
        f.write(f"                NS    S\n")
        f.write(f"  Divergence  {Dn:4d}  {Ds:3d}  (Tier 1 L4.13 core-exclusive)\n")
        f.write(f"  Polymorph.  {Pn:4d}  {Ps:3d}  (Tier 4 intra-L4.13)\n\n")
        f.write(f"Dn/Ds = {d_ratio:.3f}\n")
        f.write(f"Pn/Ps = {p_ratio:.3f}\n\n")
        f.write(f"Fisher's exact p-value (two-sided): {p_fisher:.4g}\n")
        f.write(f"Fisher's exact p-value (greater): {p_greater:.4g}\n")
        f.write(f"Direction of Selection (DoS): {d_val:+.4f}\n")
        f.write(f"95% bootstrap CI: [{ci_lo:+.4f}, {ci_hi:+.4f}]\n")
        f.write(f"Neutrality Index: {ni:.3f}\n")
        f.write(f"alpha: {alpha:+.4f}\n")

    print(f"\n[write] {out}")


if __name__ == "__main__":
    main()
