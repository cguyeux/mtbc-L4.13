"""
L4.13 Phase 5b -- Relative dating via external MTBC calibration
================================================================

The root-to-tip regression (phase 5) gave a non-interpretable signal
due to the BIN+G tree model not being suitable for nucleotide-scale
dating. Here we use an alternative: **relative dating anchored on
published MTBC node ages**.

Approach:
    1. Load the RAxML tree.
    2. Measure the branch length from the outgroup MRCA to the
       L4.13 MRCA (call it d_trunk).
    3. Measure the average distance from L4.13 MRCA to tips
       (call it d_crown).
    4. Use published MTBC reference dates to calibrate:
       - Comas 2013: MTBC MRCA ~ 6,000 ya (6 kya)
       - Pepperell 2013: L4 MRCA ~ 1,500-2,000 ya
       - Wirth 2008: MTBC MRCA estimation conservative (~ 40,000 ya)
    5. Compute L4.13 age estimates under several priors.

This is a **pragmatic scaling** approach, not a rigorous dating --
priors are shown explicitly so readers can judge the uncertainty.

Output:
    résultats/phase5b_relative_dating.txt

Usage:
    python phase5b_relative_dating.py
"""

from __future__ import annotations

from pathlib import Path
import numpy as np

from Bio import Phylo

ROOT = Path(__file__).resolve().parent.parent
TREE_DIR = ROOT / "résultats" / "phase2_tree"
RES_DIR = ROOT / "résultats"
SUPPORT_TREE = TREE_DIR / "L4.13.raxml.support"
STRAINS_TSV = TREE_DIR / "strains.tsv"


def load_groups():
    g = {}
    with STRAINS_TSV.open() as f:
        next(f)
        for line in f:
            lbl, sra, grp = line.rstrip().split("\t")
            g[lbl] = (grp, sra)
    return g


def main():
    groups = load_groups()
    tree = Phylo.read(str(SUPPORT_TREE), "newick")

    # Root on L4.4.2
    l442 = [l for l, (g, _) in groups.items() if g == "L4.4.2"]
    tree.root_with_outgroup(*l442)

    l442_mrca = tree.common_ancestor([t for t in tree.get_terminals()
                                       if t.name in l442])
    l413_tips = [t for t in tree.get_terminals()
                 if groups.get(t.name, (None,))[0] in ("L4.13.1", "L4.13.2")]
    l413_mrca = tree.common_ancestor(l413_tips)

    # MRCA of the whole ingroup + outgroup = tree root
    root = tree.root

    # Trunk = root -> L4.13 MRCA
    d_root_to_l413 = tree.distance(root, l413_mrca)

    # L4.4.2 branch from root
    d_root_to_l442 = tree.distance(root, l442_mrca)

    # Crown mean distance (L4.13 MRCA to tips)
    crown_dists = [tree.distance(l413_mrca, t) for t in l413_tips]
    d_crown_mean = float(np.mean(crown_dists))
    d_crown_max = float(np.max(crown_dists))
    d_crown_med = float(np.median(crown_dists))

    # L4.13.1 vs L4.13.2 sub-clade split depth
    l131_tips = [t for t in tree.get_terminals()
                 if groups.get(t.name, (None,))[0] == "L4.13.1"]
    l132_tips = [t for t in tree.get_terminals()
                 if groups.get(t.name, (None,))[0] == "L4.13.2"]
    l131_mrca = tree.common_ancestor(l131_tips)
    l132_mrca = tree.common_ancestor(l132_tips)
    d_131 = tree.distance(l413_mrca, l131_mrca)
    d_132 = tree.distance(l413_mrca, l132_mrca)

    lines = []
    lines.append("L4.13 -- Relative dating via external MTBC calibration\n")
    lines.append("=" * 60 + "\n\n")
    lines.append("Distances (binary BIN+G tree, subst/site):\n")
    lines.append(f"  root -> L4.13 MRCA:       {d_root_to_l413:.4f}  (trunk)\n")
    lines.append(f"  root -> L4.4.2 MRCA:      {d_root_to_l442:.4f}  (outgroup arm)\n")
    lines.append(f"  L4.13 MRCA -> tips (mean):{d_crown_mean:.4f}  (crown)\n")
    lines.append(f"  L4.13 MRCA -> tips (med): {d_crown_med:.4f}\n")
    lines.append(f"  L4.13 MRCA -> tips (max): {d_crown_max:.4f}\n")
    lines.append(f"  L4.13 MRCA -> L4.13.1 MRCA: {d_131:.4f}\n")
    lines.append(f"  L4.13 MRCA -> L4.13.2 MRCA: {d_132:.4f}\n\n")

    # Proportions
    total = d_root_to_l413 + d_crown_mean
    p_trunk = d_root_to_l413 / total
    p_crown = d_crown_mean / total
    lines.append(f"Ratio trunk / crown = {d_root_to_l413/d_crown_mean:.2f}\n")
    lines.append(f"Trunk fraction of total (root -> tips):"
                 f" {p_trunk*100:.1f}%\n")
    lines.append(f"Crown fraction: {p_crown*100:.1f}%\n\n")

    # Calibration scenarios
    # Scenario 1: outgroup MRCA = L4 MRCA assumed at 1750 ya (Pepperell 2013)
    # Tree is rooted at midpoint of outgroup arm, so root ~ MRCA(L4.4.2, ingroup)
    # which is close to L4 MRCA in principle.
    scenarios = [
        ("Pepperell 2013 (L4 MRCA ~ 1750 ya)", 1750),
        ("Bos 2014 (L4 MRCA ~ 1000 ya, conservative)", 1000),
        ("Menardo 2019 (L4 MRCA ~ 2500 ya, extended)", 2500),
    ]

    lines.append("Calibration scenarios\n")
    lines.append("-" * 60 + "\n")
    lines.append("Assumption: the tree root approximates the L4 MRCA,\n")
    lines.append("which is calibrated by literature priors.\n\n")

    # Total tree depth from root to mean tip
    total_depth = d_root_to_l413 + d_crown_mean
    for label, age_root in scenarios:
        # Rate = age_root / (total substitutions per site to tip)
        # (i.e. years per subst/site)
        # But since root is LMRCA(outgroup + ingroup), trunk covers L4 -> L4.13
        years_per_subst = age_root / total_depth
        # L4.13 MRCA age = crown_mean * years_per_subst
        age_l413 = d_crown_mean * years_per_subst
        age_131 = (d_crown_mean - d_131) * years_per_subst
        age_132 = (d_crown_mean - d_132) * years_per_subst
        # If negative, set to 0 (sub-clade recent)
        age_131 = max(0, age_131)
        age_132 = max(0, age_132)

        lines.append(f"\n  {label}:\n")
        lines.append(f"    assumed root age:       {age_root:>6} years ago\n")
        lines.append(f"    implied rate (subst/site/year): "
                     f"{total_depth/age_root:.2e}\n")
        lines.append(f"    L4.13 MRCA age:         {age_l413:>6.0f} years ago"
                     f" (~{2026-int(age_l413)} CE)\n")
        lines.append(f"    L4.13.1 MRCA age:       {age_131:>6.0f} years ago"
                     f" (~{2026-int(age_131)} CE)\n")
        lines.append(f"    L4.13.2 MRCA age:       {age_132:>6.0f} years ago"
                     f" (~{2026-int(age_132)} CE)\n")

    lines.append("\nImportant caveats\n")
    lines.append("-" * 60 + "\n")
    lines.append("  1. BIN+G tree distances are NOT directly comparable to\n")
    lines.append("     published nucleotide-scale substitution rates. The\n")
    lines.append("     above estimates assume the *relative* proportions on\n")
    lines.append("     the tree are accurate, which is broadly true but not\n")
    lines.append("     exact.\n")
    lines.append("  2. Root-to-tip regression (phase 5) showed R^2 = 0.09,\n")
    lines.append("     indicating weak temporal signal. A rigorous date\n")
    lines.append("     requires a nucleotide alignment + LSD2 or BEAST2.\n")
    lines.append("  3. Scenarios differ by >2x due to uncertainty on the\n")
    lines.append("     L4 MRCA age itself, which the MTBC phylogenetics\n")
    lines.append("     community has yet to converge on.\n")
    lines.append("  4. Therefore these numbers are *orders of magnitude*\n")
    lines.append("     rather than point estimates. Use as heuristic only.\n")

    txt = "".join(lines)
    print(txt)
    out = RES_DIR / "phase5b_relative_dating.txt"
    out.write_text(txt)
    print(f"\n[write] {out}")


if __name__ == "__main__":
    main()
