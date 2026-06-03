"""
L4.13 Analysis -- Phase 2b: Tree topology analysis
===================================================

Parses the RAxML support tree and reports:

1. Monophyly of L4.13, L4.13.1, L4.13.2, outgroup clades
2. Bootstrap support for the key bipartitions
3. Placement of the 2 outlier strains (ERR2517475, SRR33445629)
4. Detection of internal sub-clades within L4.13.1 and L4.13.2
5. Location of the 14 identical-sequence pairs
6. Root issue (why outgroup is not monophyletic)

Requires: ete3

Usage:
    python phase2b_analyze_tree.py
"""

from __future__ import annotations

from collections import Counter, defaultdict
from pathlib import Path

from ete3 import Tree

ROOT = Path(__file__).resolve().parent.parent
TREE_DIR = ROOT / "résultats" / "phase2_tree"
RESULTS_DIR = ROOT / "résultats"
OUTPUT = RESULTS_DIR / "phase2b_tree_analysis.txt"

SUPPORT_TREE = TREE_DIR / "L4.13.raxml.support"
STRAINS_TSV = TREE_DIR / "strains.tsv"

OUTLIER_SRAS = {"ERR2517475", "SRR33445629"}

def load_strain_groups():
    """Return dict label -> group ('L4.13.1', 'L4.13.2', 'L4.4.2', 'L4.2.1')."""
    groups = {}
    with open(STRAINS_TSV) as f:
        next(f)  # header
        for line in f:
            lbl, sra, g = line.rstrip("\n").split("\t")
            groups[lbl] = g
    return groups


def lines(s=""):
    out.append(s + "\n")


def is_monophyletic(tree, labels):
    """Check if the set of labels is monophyletic in the tree."""
    if len(labels) < 2:
        return True, None
    try:
        mrca = tree.get_common_ancestor(labels)
        desc = {leaf.name for leaf in mrca.get_leaves()}
        return desc == set(labels), mrca
    except Exception as e:
        return False, None


def main():
    print("=" * 70)
    print("L4.13 Phase 2b -- Tree topology analysis")
    print("=" * 70)

    global out
    out = []

    groups = load_strain_groups()
    by_group = defaultdict(list)
    for lbl, g in groups.items():
        by_group[g].append(lbl)

    lines("L4.13 Phase 2b -- Tree topology analysis")
    lines("=" * 60)
    lines()
    lines(f"Support tree: {SUPPORT_TREE.name}")
    lines(f"Total taxa: {sum(len(v) for v in by_group.values())}")
    for g in sorted(by_group):
        lines(f"  {g}: {len(by_group[g])}")
    lines()

    # Load tree -- format 0 means internal node names are support values
    t = Tree(str(SUPPORT_TREE), format=0)
    leaves = [l.name for l in t.get_leaves()]
    lines(f"Leaves in tree: {len(leaves)}")
    lines()

    # --- Root check ---
    lines("=" * 60)
    lines("ROOTING")
    lines("=" * 60)
    lines("RAxML reports: 'Outgroup is not monophyletic' -- investigate.")
    lines()
    outgroup_labels = by_group["L4.4.2"] + by_group["L4.2.1"]
    ok, mrca = is_monophyletic(t, outgroup_labels)
    lines(f"Outgroup (L4.4.2 + L4.2.1, n={len(outgroup_labels)}) monophyletic: {ok}")
    if not ok:
        # Check each sub-outgroup alone
        for sub in ["L4.4.2", "L4.2.1"]:
            ok_sub, mrca_sub = is_monophyletic(t, by_group[sub])
            lines(f"  - {sub} alone (n={len(by_group[sub])}): "
                  f"monophyletic = {ok_sub}")
            if ok_sub and mrca_sub is not None:
                desc = [l.name for l in mrca_sub.get_leaves()]
                lines(f"      MRCA descendants: {len(desc)}")
    # Reroot on L4.4.2 (sister of L4.13) if monophyletic
    l442 = by_group["L4.4.2"]
    ok_l442, mrca_l442 = is_monophyletic(t, l442)
    if ok_l442:
        t.set_outgroup(mrca_l442)
        lines("→ Re-rooted tree on L4.4.2 MRCA.")
    else:
        # Use midpoint as fallback
        t.set_outgroup(t.get_midpoint_outgroup())
        lines("→ Re-rooted on midpoint (outgroups not monophyletic).")
    lines()

    # --- Monophyly checks ---
    lines("=" * 60)
    lines("MONOPHYLY")
    lines("=" * 60)
    for g in ["L4.13.1", "L4.13.2", "L4.4.2", "L4.2.1"]:
        ok, mrca = is_monophyletic(t, by_group[g])
        status = "monophyletic" if ok else "NOT monophyletic"
        lines(f"  {g} (n={len(by_group[g])}): {status}")
        if ok and mrca is not None:
            sup = mrca.support
            lines(f"    MRCA bootstrap support: {sup:.0f}")
        if not ok and mrca is not None:
            desc = {l.name for l in mrca.get_leaves()}
            foreign = desc - set(by_group[g])
            lines(f"    MRCA contains {len(foreign)} foreign leaves: "
                  f"{sorted(foreign)[:10]}"
                  + (f" (+{len(foreign)-10} more)" if len(foreign) > 10 else ""))
    lines()

    # L4.13 = 1+2 together
    l413 = by_group["L4.13.1"] + by_group["L4.13.2"]
    ok, mrca = is_monophyletic(t, l413)
    lines(f"  L4.13 (= L4.13.1 ∪ L4.13.2, n={len(l413)}): "
          f"{'monophyletic' if ok else 'NOT monophyletic'}")
    if ok and mrca is not None:
        lines(f"    MRCA bootstrap support: {mrca.support:.0f}")
    lines()

    # --- Outlier placement ---
    lines("=" * 60)
    lines("OUTLIERS PLACEMENT")
    lines("=" * 60)
    for sra in OUTLIER_SRAS:
        matches = [l for l in leaves if sra in l]
        if not matches:
            lines(f"  {sra}: NOT FOUND in tree")
            continue
        lbl = matches[0]
        node = t & lbl
        lines(f"  {lbl}")
        # Ancestor groups: walk up
        anc = node.up
        depth = 0
        while anc is not None and depth < 5:
            leaves_below = [l.name for l in anc.get_leaves()]
            grps = Counter(groups.get(x, "?") for x in leaves_below)
            lines(f"    ancestor #{depth}: support={anc.support:.0f}, "
                  f"n_leaves={len(leaves_below)}, groups={dict(grps)}")
            anc = anc.up
            depth += 1
    lines()

    # --- Support for key clades ---
    lines("=" * 60)
    lines("KEY CLADE SUPPORTS")
    lines("=" * 60)
    lines("Bootstrap support at the MRCA of each group:")
    for g, labels in [
        ("L4.13.1", by_group["L4.13.1"]),
        ("L4.13.2", by_group["L4.13.2"]),
        ("L4.13",   l413),
    ]:
        ok, mrca = is_monophyletic(t, labels)
        if ok and mrca is not None:
            lines(f"  {g:<10s}: support = {mrca.support:5.0f}   "
                  f"(branch length to parent = {mrca.dist:.6f})")
        else:
            lines(f"  {g:<10s}: non-monophyletic, no single MRCA to score")
    lines()

    # --- Internal sub-structure within L4.13.2 ---
    lines("=" * 60)
    lines("INTERNAL SUB-STRUCTURE WITHIN L4.13.2 (197 strains)")
    lines("=" * 60)
    ok, mrca132 = is_monophyletic(t, by_group["L4.13.2"])
    if ok and mrca132 is not None:
        # Collect direct children that are themselves cladal
        internal_clades = []
        for child in mrca132.children:
            n = len(child.get_leaves())
            sup = child.support
            if n >= 5:
                internal_clades.append((n, sup, child))
        internal_clades.sort(key=lambda x: -x[0])
        lines(f"Direct children of L4.13.2 MRCA with >=5 leaves:")
        for n, sup, node in internal_clades[:10]:
            lines(f"    support={sup:5.0f}  n_leaves={n}")
        # Deeper: find all internal nodes with >=10 leaves and bootstrap >=80
        well_supported = []
        for node in mrca132.traverse():
            if node == mrca132 or node.is_leaf():
                continue
            n = len(node.get_leaves())
            if n >= 10 and node.support >= 80:
                well_supported.append((n, node.support, node))
        well_supported.sort(key=lambda x: -x[0])
        lines(f"\nWell-supported internal clades (n>=10, bootstrap>=80):")
        for n, sup, node in well_supported[:15]:
            # Sample some leaves
            leaves_here = [l.name.split("_", 1)[1] for l in node.get_leaves()][:3]
            lines(f"    support={sup:5.0f}  n={n:<4d}  e.g. {leaves_here}")
    lines()

    # --- Internal sub-structure within L4.13.1 ---
    lines("=" * 60)
    lines("INTERNAL SUB-STRUCTURE WITHIN L4.13.1 (61 strains)")
    lines("=" * 60)
    ok, mrca131 = is_monophyletic(t, by_group["L4.13.1"])
    if ok and mrca131 is not None:
        well_supported = []
        for node in mrca131.traverse():
            if node == mrca131 or node.is_leaf():
                continue
            n = len(node.get_leaves())
            if n >= 5 and node.support >= 80:
                well_supported.append((n, node.support, node))
        well_supported.sort(key=lambda x: -x[0])
        lines(f"Well-supported internal clades (n>=5, bootstrap>=80):")
        for n, sup, node in well_supported[:10]:
            leaves_here = [l.name.split("_", 1)[1] for l in node.get_leaves()][:3]
            lines(f"    support={sup:5.0f}  n={n:<4d}  e.g. {leaves_here}")
    lines()

    # --- Branch length between L4.13.1 and L4.13.2 ---
    lines("=" * 60)
    lines("BIFURCATION DEPTH (L4.13.1 vs L4.13.2)")
    lines("=" * 60)
    ok1, m1 = is_monophyletic(t, by_group["L4.13.1"])
    ok2, m2 = is_monophyletic(t, by_group["L4.13.2"])
    if ok1 and ok2:
        # Branch length from L4.13 MRCA to each sub-clade MRCA
        ok_all, mrca_all = is_monophyletic(t, l413)
        if ok_all:
            # distance from common ancestor down to each sub MRCA
            d1 = mrca_all.get_distance(m1)
            d2 = mrca_all.get_distance(m2)
            lines(f"  Distance L4.13 MRCA -> L4.13.1 MRCA: {d1:.5f}")
            lines(f"  Distance L4.13 MRCA -> L4.13.2 MRCA: {d2:.5f}")
            lines(f"  Ratio (L4.13.2/L4.13.1): {d2/d1 if d1>0 else '(div by 0)'}")
    lines()

    # Write
    OUTPUT.write_text("".join(out))
    print("".join(out))
    print(f"\n[write] {OUTPUT}")


if __name__ == "__main__":
    main()
