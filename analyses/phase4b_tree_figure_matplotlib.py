"""
L4.13 Phase 4b -- Tree figure via matplotlib + Biopython (cleaner)
===================================================================

Produces a cleaner rectangular PDF of the RAxML tree:
    - Leaves: one coloured square per sub-clade, no text
    - Internal nodes: bootstrap shown only if >= 80
    - Outliers and basal candidate highlighted
    - Branches colour-coded by resolved sub-clade

Output:
    article/figures/phase4_tree_rect.pdf
    article/figures/phase4_tree_rect.png
"""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from Bio import Phylo

ROOT = Path(__file__).resolve().parent.parent
TREE_DIR = ROOT / "résultats" / "phase2_tree"
FIG_DIR = ROOT / "article" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

SUPPORT_TREE = TREE_DIR / "L4.13.raxml.support"
STRAINS_TSV = TREE_DIR / "strains.tsv"
STRAINS_CSV = ROOT / "data" / "strains.csv"

COLOURS = {
    "L4.13.1": "#1f77b4",
    "L4.13.2": "#d62728",
    "L4.4.2":  "#555555",
    "L4.2.1":  "#999999",
}
OUTGROUP = {"L4.4.2", "L4.2.1"}

HIGHLIGHTS = {
    "SRR33894822": ("#ffbf00", "basal candidate"),
    "SRR33445629": ("#ff7f0e", "multi-system outlier"),
    "ERR2517475":  ("#ff7f0e", "multi-system outlier"),
}


def load_groups():
    g = {}
    with STRAINS_TSV.open() as f:
        next(f)
        for line in f:
            lbl, sra, grp = line.rstrip().split("\t")
            g[lbl] = (grp, sra)
    return g


def build_positions(tree):
    """Return dict clade -> (x, y) post-rooting."""
    x = {tree.root: 0.0}
    for clade in tree.find_clades(order="preorder"):
        for child in clade.clades:
            x[child] = x[clade] + (child.branch_length or 0)
    terminals = tree.get_terminals()
    y = {t: i for i, t in enumerate(terminals)}
    def resolve(c):
        if c in y:
            return y[c]
        ys = [resolve(ch) for ch in c.clades]
        y[c] = (min(ys) + max(ys)) / 2
        return y[c]
    resolve(tree.root)
    return x, y


def clade_group(clade, groups):
    """Return the sub-clade group if all descendants share it, else None."""
    terms = [t.name for t in clade.get_terminals()]
    gs = {groups[t][0] for t in terms if t in groups}
    return list(gs)[0] if len(gs) == 1 else None


def main():
    print("=" * 70)
    print("L4.13 Phase 4b -- matplotlib tree (clean)")
    print("=" * 70)

    groups = load_groups()
    tree = Phylo.read(str(SUPPORT_TREE), "newick")

    l442 = [l for l, (g, _) in groups.items() if g == "L4.4.2"]
    tree.root_with_outgroup(*l442)

    n_by_grp = Counter(g for _, (g, _) in groups.items())
    print(f"[count] {dict(n_by_grp)}")

    x, y = build_positions(tree)
    fig, ax = plt.subplots(figsize=(7.5, 18))

    # Draw branches
    for clade in tree.find_clades():
        parent = None
        for p in tree.find_clades():
            if clade in p.clades:
                parent = p
                break
        grp = clade_group(clade, groups)
        colour = COLOURS.get(grp, "#000000") if grp else "#000000"
        lw = 0.7
        if parent is not None:
            ax.plot([x[parent], x[clade]], [y[clade], y[clade]],
                    color=colour, linewidth=lw, solid_capstyle="butt")
        # Vertical connector between children
        if clade.clades:
            ys = [y[c] for c in clade.clades]
            ax.plot([x[clade]] * 2, [min(ys), max(ys)],
                    color=colour, linewidth=lw, solid_capstyle="butt")

    # Leaf markers
    xmax = max(x.values())
    strip1 = xmax * 1.01
    strip2 = xmax * 1.04

    for t in tree.get_terminals():
        grp, sra = groups.get(t.name, ("?", ""))
        c = COLOURS.get(grp, "#000000")
        ax.plot(strip1, y[t], marker="s", markersize=2.2, color=c,
                markeredgewidth=0)
        if sra in HIGHLIGHTS:
            hc, _ = HIGHLIGHTS[sra]
            ax.plot(strip2, y[t], marker="o", markersize=4.5, color=hc,
                    markeredgecolor="black", markeredgewidth=0.3)

    # Bootstrap labels on internal nodes: show only if >=80 AND has >=5
    # descendants
    for clade in tree.get_nonterminals():
        if clade == tree.root or clade.confidence is None:
            continue
        n_desc = len(clade.get_terminals())
        if clade.confidence < 80 or n_desc < 5:
            continue
        ax.text(x[clade], y[clade], f"{int(clade.confidence)}",
                fontsize=5.5, va="center", ha="right",
                color="#333333",
                bbox=dict(boxstyle="round,pad=0.12",
                          facecolor="white", edgecolor="none", alpha=0.85))

    # Cosmetics
    ax.set_yticks([])
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.set_xlabel("Substitutions per site", fontsize=9)
    ax.set_title("L4.13 ML phylogeny (268 taxa, BIN+G, 100 BS)",
                 fontsize=10, pad=8)
    ax.tick_params(axis="x", labelsize=8)
    ax.set_xlim(-0.05, xmax * 1.10)

    # Legend
    handles = [
        Patch(color=COLOURS["L4.13.1"], label=f"L4.13.1 (n={n_by_grp['L4.13.1']})"),
        Patch(color=COLOURS["L4.13.2"], label=f"L4.13.2 (n={n_by_grp['L4.13.2']})"),
        Patch(color=COLOURS["L4.4.2"],
              label=f"L4.4.2 outgroup (n={n_by_grp['L4.4.2']})"),
        Patch(color=COLOURS["L4.2.1"],
              label=f"L4.2.1 outgroup (n={n_by_grp['L4.2.1']})"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#ffbf00",
               markeredgecolor="black", markersize=6,
               label="Basal candidate (SRR33894822)"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#ff7f0e",
               markeredgecolor="black", markersize=6,
               label="Multi-system outlier"),
    ]
    ax.legend(handles=handles, loc="lower right", fontsize=7,
              frameon=True, framealpha=0.92)

    out_pdf = FIG_DIR / "phase4_tree_rect.pdf"
    out_png = FIG_DIR / "phase4_tree_rect.png"
    plt.tight_layout()
    plt.savefig(out_pdf, dpi=600, bbox_inches="tight")
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    print(f"[write] {out_pdf}")
    print(f"[write] {out_png}")
    plt.close()

    # Summary of key supports printed to text
    print("\n[key supports]")
    # L4.13 MRCA = non-trivial ancestor
    l413_leaves = [t for t in tree.get_terminals()
                   if groups.get(t.name, (None,))[0] in ("L4.13.1", "L4.13.2")]
    # Helper: find MRCA of a set
    def mrca(leaves):
        return tree.common_ancestor(leaves)
    m413 = mrca(l413_leaves)
    m131 = mrca([t for t in tree.get_terminals()
                 if groups.get(t.name, (None,))[0] == "L4.13.1"])
    m132 = mrca([t for t in tree.get_terminals()
                 if groups.get(t.name, (None,))[0] == "L4.13.2"])
    print(f"  L4.13 clade bs = {m413.confidence}")
    print(f"  L4.13.1 clade bs = {m131.confidence}")
    print(f"  L4.13.2 MRCA bs = {m132.confidence}"
          f"  (includes {len(m132.get_terminals())} leaves)")

    print("\n" + "=" * 70)
    print("Phase 4b complete.")
    print("=" * 70)


if __name__ == "__main__":
    main()
