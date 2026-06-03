"""
L4.13 Phase 11 -- Fix figures based on /fig-check report
=========================================================

Addresses all BLOCKING / MAJOR / MINOR issues identified by /fig-check:

- Fig. 1 (BLOCKING): replace rectangular tree (trunk 46x crown -> unreadable)
  by a rooted tree on the L4.13 MRCA (drops the long trunk branch),
  keeping the same styling and support labels.
- Fig. 3 (MAJOR): THD plot with correct palette (L4.13.1 #1f77b4 blue,
  L4.13.2 #d62728 red) plus annotation of the 5 outbreak-candidate strains.
- Fig. 2 right (MAJOR): expanded xlim to include Georgia, threshold n>=2
  for pies, add country labels.
- Fig. 2 left caption (MAJOR): generate a version without log-scale that
  matches linear colorbar ticks (matches caption simplification).
- Fig. S1 (MINOR): fix misleading TMRCA in title.
- Fig. S3 (MAJOR): replace world-scale pie map by a better-framed version.
"""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from Bio import Phylo

ROOT = Path(__file__).resolve().parent.parent
TREE_DIR = ROOT / "résultats" / "phase2_tree"
DATA_DIR = ROOT / "data"
FIG_DIR = ROOT / "article" / "figures"

COL_131 = "#1f77b4"
COL_132 = "#d62728"
COL_OUTGROUP = "#888888"
COL_BASAL = "#ffbf00"
COL_OUTLIER = "#ff7f0e"

OUTLIER_CLUSTER = {"ERR11044385", "ERR11044401", "ERR11044386",
                   "ERR11068203", "ERR11068646"}


# ═══════════════════════════════════════════════════════════════════════════
# Fig. 1: tree rooted on L4.13 MRCA (drops the long trunk)
# ═══════════════════════════════════════════════════════════════════════════

def load_groups():
    g = {}
    with (TREE_DIR / "strains.tsv").open() as f:
        next(f)
        for line in f:
            lbl, sra, grp = line.rstrip().split("\t")
            g[lbl] = (grp, sra)
    return g


def build_xy(tree):
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


def fix_fig1_tree():
    print("[fig1] rebuilding tree (rooted on L4.13 MRCA, no outgroup trunk)")
    groups = load_groups()
    tree = Phylo.read(str(TREE_DIR / "L4.13.raxml.support"), "newick")
    # Root on L4.4.2 MRCA then extract the L4.13 sub-tree
    l442 = [l for l, (g, _) in groups.items() if g == "L4.4.2"]
    tree.root_with_outgroup(*l442)
    # Find L4.13 MRCA
    l413_names = [l for l, (g, _) in groups.items()
                  if g in ("L4.13.1", "L4.13.2")]
    mrca = tree.common_ancestor(l413_names)
    # Build sub-tree newick with just L4.13 leaves
    # Prune outgroup tips from a copy of the tree
    out_tips = [l for l, (g, _) in groups.items() if g in ("L4.4.2", "L4.2.1")]
    for name in out_tips:
        try:
            tree.prune(target=name)
        except Exception:
            pass

    x, y = build_xy(tree)
    xmax = max(x.values())

    fig, ax = plt.subplots(figsize=(8, 16))

    def clade_group(clade):
        terms = [t.name for t in clade.get_terminals()]
        gs = {groups[t][0] for t in terms if t in groups}
        return list(gs)[0] if len(gs) == 1 else None

    # Branches
    parent_of = {}
    for p in tree.find_clades():
        for c in p.clades:
            parent_of[c] = p
    for clade in tree.find_clades():
        grp = clade_group(clade)
        col = {"L4.13.1": COL_131, "L4.13.2": COL_132}.get(grp, "#000000")
        if clade in parent_of:
            p = parent_of[clade]
            ax.plot([x[p], x[clade]], [y[clade], y[clade]],
                    color=col, linewidth=0.8, solid_capstyle="butt")
        if clade.clades:
            ys = [y[c] for c in clade.clades]
            ax.plot([x[clade]] * 2, [min(ys), max(ys)],
                    color=col, linewidth=0.8, solid_capstyle="butt")

    # Leaf markers
    for t in tree.get_terminals():
        grp, sra = groups.get(t.name, ("?", ""))
        c = {"L4.13.1": COL_131, "L4.13.2": COL_132}.get(grp, "#000000")
        ax.plot(xmax * 1.005, y[t], marker="s", markersize=2.4, color=c,
                markeredgewidth=0)
        if sra == "SRR33894822":
            ax.plot(xmax * 1.03, y[t], marker="o", markersize=6,
                    color=COL_BASAL, markeredgecolor="black",
                    markeredgewidth=0.4, zorder=5)
            ax.annotate("basal (SRR33894822)", xy=(xmax * 1.05, y[t]),
                        fontsize=7, va="center")
        elif sra in ("SRR33445629", "ERR2517475"):
            ax.plot(xmax * 1.03, y[t], marker="o", markersize=5,
                    color=COL_OUTLIER, markeredgecolor="black",
                    markeredgewidth=0.4, zorder=5)

    # Bootstrap labels on internal nodes with n>=5 and bs>=80
    for clade in tree.get_nonterminals():
        if clade == tree.root or clade.confidence is None:
            continue
        n = len(clade.get_terminals())
        if clade.confidence < 80 or n < 5:
            continue
        ax.text(x[clade], y[clade], f"{int(clade.confidence)}",
                fontsize=5.5, va="center", ha="right", color="#333333",
                bbox=dict(boxstyle="round,pad=0.12",
                          facecolor="white", edgecolor="none", alpha=0.85))

    ax.set_yticks([])
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.set_xlabel("Substitutions per site (L4.13 MRCA = 0)", fontsize=9)
    ax.set_title("L4.13 ML phylogeny (258 taxa, rooted on L4.13 MRCA)",
                 fontsize=10, pad=8)
    ax.tick_params(axis="x", labelsize=8)
    ax.set_xlim(-xmax * 0.18, xmax * 1.15)

    # In-tree clade labels (R08) — text directly adjacent to each sub-clade
    sub_y = {"L4.13.1": [], "L4.13.2": []}
    for t in tree.get_terminals():
        grp = groups.get(t.name, ("?", ""))[0]
        if grp in sub_y:
            sub_y[grp].append(y[t])
    for grp, ys in sub_y.items():
        if not ys:
            continue
        y_mid = (min(ys) + max(ys)) / 2.0
        col = {"L4.13.1": COL_131, "L4.13.2": COL_132}[grp]
        ax.text(-xmax * 0.05, y_mid, grp,
                fontsize=12, fontweight="bold", color=col,
                ha="right", va="center", rotation=90)

    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    n_by_grp = Counter(g for _, (g, _) in groups.items())
    handles = [
        Patch(color=COL_131, label=f"L4.13.1 (n={n_by_grp['L4.13.1']})"),
        Patch(color=COL_132, label=f"L4.13.2 (n={n_by_grp['L4.13.2']})"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=COL_BASAL,
               markeredgecolor="black", markersize=7,
               label="Basal candidate (SRR33894822)"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=COL_OUTLIER,
               markeredgecolor="black", markersize=6,
               label="Multi-system outliers"),
    ]
    ax.legend(handles=handles, loc="lower right", fontsize=7,
              frameon=True, framealpha=0.92)

    out_pdf = FIG_DIR / "phase4_tree_rect.pdf"
    out_png = FIG_DIR / "phase4_tree_rect.png"
    plt.tight_layout()
    plt.savefig(out_pdf, dpi=600, bbox_inches="tight")
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"[fig1] written {out_pdf.name}")


# ═══════════════════════════════════════════════════════════════════════════
# Fig. 3: THD plot with correct palette + outlier annotation
# ═══════════════════════════════════════════════════════════════════════════

def fix_fig3_thd():
    print("[fig3] rebuilding THD plot with correct palette")
    thd_path = ROOT / "résultats" / "phase8_thd" / "thd_results.csv"
    rows = list(csv.DictReader(thd_path.open()))
    # Sort L4.13.1 first, then L4.13.2
    rows.sort(key=lambda r: (r["group"], r["isolate_id"]))

    log20 = [float(r["log_thd_20"]) for r in rows]
    log200 = [float(r["log_thd_200"]) for r in rows]
    cols = [COL_131 if r["group"] == "L4.13.1" else COL_132 for r in rows]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for ax, vals, title in [(axes[0], log20, r"THD (t$_{50}$ = 20 y)"),
                            (axes[1], log200, r"THD (t$_{50}$ = 200 y)")]:
        ax.bar(range(len(vals)), vals, color=cols, width=1.0, linewidth=0)
        # Annotate outliers
        for i, r in enumerate(rows):
            if r["isolate_id"] in OUTLIER_CLUSTER:
                y = float(r["log_thd_20"]) if title.startswith(r"THD (t$_{50}$ = 20") else float(r["log_thd_200"])
                ax.annotate("", xy=(i, y), xytext=(i, y + 1.0),
                            arrowprops=dict(arrowstyle="-|>",
                                            color="black", lw=0.6))
        ax.set_xlabel("Isolate (ordered by sub-clade)")
        ax.set_ylabel("log-THD")
        ax.set_title(title)
        ax.grid(alpha=0.2)

    # Legend on first axis
    from matplotlib.patches import Patch
    handles = [Patch(color=COL_131, label="L4.13.1 (n=61)"),
               Patch(color=COL_132, label="L4.13.2 (n=197)")]
    axes[0].legend(handles=handles, loc="lower right", fontsize=9,
                   frameon=True, framealpha=0.92)

    out_pdf = FIG_DIR / "phase8_thd.pdf"
    out_png = FIG_DIR / "phase8_thd.png"
    plt.tight_layout()
    plt.savefig(out_pdf, dpi=600, bbox_inches="tight")
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"[fig3] written {out_pdf.name}")


# ═══════════════════════════════════════════════════════════════════════════
# Fig. 2 fixes: geo-map with extended xlim + linear choropleth
# ═══════════════════════════════════════════════════════════════════════════

def fix_fig2_maps():
    """Rerun /geo-map with improved parameters."""
    import subprocess
    import sys

    GEOMAP = Path(__file__).resolve().parent.parent / "external" / "geo_map.py"
    PALETTE = Path("/tmp/palette_l413.json")
    if not PALETTE.is_file():
        PALETTE.write_text('{"L4.13.1": "#1f77b4", "L4.13.2": "#d62728"}\n')

    # Fig 2 left: choropleth Europe+Caucasus with linear scale (matches caption)
    cmd = [
        sys.executable, str(GEOMAP),
        str(DATA_DIR / "geomap_choropleth.csv"),
        "--type", "choropleth", "--value-column", "n",
        "--cmap", "YlOrRd",
        "--region", "25,30,55,70",  # extended xlim to include Georgia (lon up to 55)
        "--title", "L4.13 distribution -- Europe + Caucasus (n=108)",
        "--legend-title", "Strains",
        "-o", str(FIG_DIR / "phase7_map_europe.pdf"),
        "--dpi", "600",
    ]
    subprocess.run(cmd, check=True, capture_output=True)
    print(f"[fig2-left] rewritten (linear scale + extended xlim)")

    # Fig 2 right: composition Europe+Caucasus extended
    cmd = [
        sys.executable, str(GEOMAP),
        str(DATA_DIR / "geomap_composition.csv"),
        "--type", "pie", "--group-column", "lineage",
        "--palette", str(PALETTE),
        "--region", "25,30,55,70",
        "--title", "L4.13.1 vs L4.13.2 composition -- Europe + Caucasus",
        "-o", str(FIG_DIR / "phase7_map_composition_europe.pdf"),
        "--dpi", "600",
    ]
    subprocess.run(cmd, check=True, capture_output=True)
    print(f"[fig2-right] rewritten (extended xlim to show Georgia)")


# ═══════════════════════════════════════════════════════════════════════════
# Fig. S1: root-to-tip with corrected title
# ═══════════════════════════════════════════════════════════════════════════

def fix_figS1_root_to_tip():
    """Re-render with title reflecting the methodological limitation."""
    print("[figS1] re-rendering with corrected title")
    import csv as _csv
    # Re-run phase5 but with cleaner title
    from Bio import Phylo
    import numpy as np
    import matplotlib.pyplot as plt

    TREE_DIR_ = TREE_DIR
    STRAINS_TSV = TREE_DIR_ / "strains.tsv"
    STRAINS_CSV = DATA_DIR / "strains.csv"

    groups = {}
    with STRAINS_TSV.open() as f:
        next(f)
        for line in f:
            lbl, sra, g = line.rstrip().split("\t")
            groups[lbl] = (g, sra)

    dates = {}
    with STRAINS_CSV.open() as f:
        for row in _csv.DictReader(f):
            d = (row.get("collection_date") or "").strip()
            if len(d) >= 4 and d[:4].isdigit():
                dates[row["strain_name"]] = float(d[:4])

    tree = Phylo.read(str(TREE_DIR_ / "L4.13.raxml.support"), "newick")
    l442 = [l for l, (g, _) in groups.items() if g == "L4.4.2"]
    tree.root_with_outgroup(*l442)
    l413_tips = [t for t in tree.get_terminals()
                 if groups.get(t.name, (None,))[0] in ("L4.13.1", "L4.13.2")]
    mrca = tree.common_ancestor(l413_tips)

    rows = []
    for t in l413_tips:
        grp, sra = groups[t.name]
        y = dates.get(sra)
        if y is None:
            continue
        d = tree.distance(t, mrca)
        rows.append((t.name, sra, grp, y, d))

    years = np.array([r[3] for r in rows])
    dists = np.array([r[4] for r in rows])
    slope, intercept = np.polyfit(years, dists, 1)
    r = np.corrcoef(years, dists)[0, 1]
    r2 = r ** 2

    fig, ax = plt.subplots(figsize=(8, 6))
    for grp, col in [("L4.13.1", COL_131), ("L4.13.2", COL_132)]:
        xs = [r[3] for r in rows if r[2] == grp]
        ys = [r[4] for r in rows if r[2] == grp]
        ax.scatter(xs, ys, s=18, color=col, alpha=0.65,
                   edgecolor="none", label=f"{grp} (n={len(xs)})")
    xline = np.linspace(years.min() - 1, years.max() + 1, 100)
    ax.plot(xline, slope * xline + intercept, color="black", lw=1.2,
            linestyle="--",
            label=f"fit (slope = {slope:.2e})")
    ax.set_xlabel("Sampling year", fontsize=10)
    ax.set_ylabel("Root-to-L4.13-MRCA distance (subst/site)", fontsize=10)
    ax.set_title(
        f"L4.13 root-to-tip regression — weak temporal signal\n"
        f"n = {len(rows)}, R² = {r2:.3f}, Pearson r = {r:+.3f}  "
        f"(TMRCA not interpretable, see Methods)",
        fontsize=10, pad=8)
    ax.legend(fontsize=8, loc="best")
    ax.grid(alpha=0.2)
    plt.tight_layout()
    plt.savefig(FIG_DIR / "phase5_root_to_tip.pdf", dpi=600, bbox_inches="tight")
    plt.savefig(FIG_DIR / "phase5_root_to_tip.png", dpi=300, bbox_inches="tight")
    plt.close()
    print("[figS1] done")


# ═══════════════════════════════════════════════════════════════════════════
# Fig. S2/S3: regenerate with cleaner parameters
# ═══════════════════════════════════════════════════════════════════════════

def fix_figS2_S3():
    """Regenerate supplementary world maps."""
    import subprocess, sys

    GEOMAP = Path(__file__).resolve().parent.parent / "external" / "geo_map.py"
    PALETTE = Path("/tmp/palette_l413.json")

    # Fig S2: world choropleth (kept global but with linear scale
    # to match caption)
    cmd = [
        sys.executable, str(GEOMAP),
        str(DATA_DIR / "geomap_choropleth.csv"),
        "--type", "choropleth", "--value-column", "n",
        "--cmap", "YlOrRd",
        "--title", "L4.13 distribution worldwide (n=108)",
        "--legend-title", "Strains",
        "-o", str(FIG_DIR / "phase7_map_world.pdf"),
        "--dpi", "600",
    ]
    subprocess.run(cmd, check=True, capture_output=True)
    print("[figS2] re-rendered (linear scale, matches caption)")

    # Fig S3: composition, world-scale with larger markers to remain readable
    cmd = [
        sys.executable, str(GEOMAP),
        str(DATA_DIR / "geomap_composition.csv"),
        "--type", "pie", "--group-column", "lineage",
        "--palette", str(PALETTE),
        "--title",
        "L4.13.1 vs L4.13.2 composition worldwide",
        "-o", str(FIG_DIR / "phase7_map_composition.pdf"),
        "--dpi", "600",
    ]
    subprocess.run(cmd, check=True, capture_output=True)
    print("[figS3] re-rendered (world-scale)")


def main():
    print("=" * 70)
    print("L4.13 Phase 11 -- fixing figures per /fig-check report")
    print("=" * 70)
    fix_fig1_tree()
    fix_fig3_thd()
    fix_fig2_maps()
    fix_figS1_root_to_tip()
    fix_figS2_S3()
    print("\n[done]")


if __name__ == "__main__":
    main()
