"""
L4.13 Phase 4 -- Publication-ready tree figure
===============================================

Builds two versions of the L4.13 RAxML tree with colour-coded sub-clades
and ancillary annotations:

    1. Rectangular cladogram with bootstrap support (>=70) at internal
       nodes, L4.13.1/L4.13.2 coloured strips, country annotation if
       available, outlier/basal strain marks.
    2. Circular tree (same annotations) for compact display.

Also emits iTOL-compatible dataset files (COLORSTRIP + TEXT) so the
tree can be refined interactively on the iTOL server.

Dependencies: ete3, PyQt5/QtPDF (for PDF export). Falls back to SVG if
PDF export fails.

Outputs:
    article/figures/phase4_tree_rect.pdf
    article/figures/phase4_tree_circ.pdf
    article/figures/itol/newick.nwk
    article/figures/itol/colorstrip_subclade.txt
    article/figures/itol/label_country.txt
    article/figures/itol/annotations_outliers.txt

Usage:
    python phase4_tree_figure.py
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

from ete3 import Tree, TreeStyle, NodeStyle, TextFace, faces

ROOT = Path(__file__).resolve().parent.parent
TREE_DIR = ROOT / "résultats" / "phase2_tree"
FIG_DIR = ROOT / "article" / "figures"
ITOL_DIR = FIG_DIR / "itol"
FIG_DIR.mkdir(parents=True, exist_ok=True)
ITOL_DIR.mkdir(parents=True, exist_ok=True)

SUPPORT_TREE = TREE_DIR / "L4.13.raxml.support"
STRAINS_TSV = TREE_DIR / "strains.tsv"
STRAINS_CSV = ROOT / "data" / "strains.csv"

COLOURS = {
    "L4.13.1": "#1f77b4",  # blue
    "L4.13.2": "#d62728",  # red
    "L4.4.2":  "#7f7f7f",  # dark grey outgroup
    "L4.2.1":  "#bcbcbc",  # light grey outgroup
}

HIGHLIGHTS = {
    "SRR33894822": "#ffbf00",  # basal candidate (orange)
    "SRR33445629": "#ff7f0e",  # multi-system outlier (darker orange)
    "ERR2517475":  "#ff7f0e",  # multi-system outlier
}


def load_groups():
    groups = {}
    with STRAINS_TSV.open() as f:
        next(f)
        for line in f:
            lbl, sra, g = line.rstrip().split("\t")
            groups[lbl] = (g, sra)
    return groups


def load_country_map():
    """Return sra -> country short label."""
    country = {}
    if not STRAINS_CSV.is_file():
        return country
    with STRAINS_CSV.open() as f:
        for row in csv.DictReader(f):
            c = (row.get("country") or "").split(":")[0].strip()
            if c:
                country[row["strain_name"]] = c
    return country


def layout_leaf(node):
    """Custom leaf layout adding coloured strip and optional country tag."""
    if not node.is_leaf():
        return
    grp, sra = node.grp, node.sra
    colour = COLOURS.get(grp, "#000000")

    # Thin coloured rectangle next to leaf name
    strip = faces.RectFace(4, 10, colour, colour)
    strip.margin_left = 4
    faces.add_face_to_node(strip, node, column=1, position="aligned")

    # Highlight the outliers
    if sra in HIGHLIGHTS:
        mark = faces.CircleFace(5, HIGHLIGHTS[sra])
        mark.margin_left = 2
        faces.add_face_to_node(mark, node, column=2, position="aligned")

    # Country tag (truncated)
    c = node.country
    if c:
        tf = TextFace(c[:3].upper(), fsize=7, fgcolor="#404040")
        tf.margin_left = 4
        faces.add_face_to_node(tf, node, column=3, position="aligned")


def layout_internal_bs(node):
    """Show bootstrap support >=70 on internal nodes."""
    if node.is_leaf() or node.support is None:
        return
    if node.support >= 70:
        tf = TextFace(f"{int(node.support)}", fsize=6, fgcolor="#333333")
        tf.margin_right = 2
        faces.add_face_to_node(tf, node, column=0, position="branch-top")


def main():
    print("=" * 70)
    print("L4.13 Phase 4 -- tree figure")
    print("=" * 70)

    groups = load_groups()
    country_map = load_country_map()

    t = Tree(str(SUPPORT_TREE), format=0)
    # Reroot on L4.4.2 MRCA
    l442 = [l for l, (g, _) in groups.items() if g == "L4.4.2"]
    t.set_outgroup(t.get_common_ancestor(l442))

    # Attach metadata to each leaf
    for leaf in t.get_leaves():
        grp, sra = groups.get(leaf.name, ("?", leaf.name))
        leaf.grp = grp
        leaf.sra = sra
        leaf.country = country_map.get(sra, "")

    # --- iTOL export ---
    (ITOL_DIR / "newick.nwk").write_text(t.write(format=0))

    # COLORSTRIP dataset for sub-clade
    colorstrip = [
        "DATASET_COLORSTRIP",
        "SEPARATOR TAB",
        "DATASET_LABEL\tSubclade",
        "COLOR\t#1f77b4",
        "STRIP_WIDTH\t30",
        "LEGEND_TITLE\tSub-clade",
        "LEGEND_SHAPES\t1\t1\t1\t1",
        "LEGEND_COLORS\t{}\t{}\t{}\t{}".format(
            COLOURS["L4.13.1"], COLOURS["L4.13.2"],
            COLOURS["L4.4.2"], COLOURS["L4.2.1"]),
        "LEGEND_LABELS\tL4.13.1\tL4.13.2\tL4.4.2 (outgroup)\tL4.2.1 (outgroup)",
        "DATA",
    ]
    for leaf in t.get_leaves():
        colorstrip.append(f"{leaf.name}\t{COLOURS.get(leaf.grp, '#000000')}\t{leaf.grp}")
    (ITOL_DIR / "colorstrip_subclade.txt").write_text("\n".join(colorstrip) + "\n")

    # TEXT dataset for country
    if country_map:
        labels = [
            "DATASET_TEXT",
            "SEPARATOR TAB",
            "DATASET_LABEL\tCountry",
            "COLOR\t#606060",
            "DATA",
        ]
        for leaf in t.get_leaves():
            if leaf.country:
                labels.append(f"{leaf.name}\t{leaf.country}\t-1\t#404040\tnormal\t1")
        (ITOL_DIR / "label_country.txt").write_text("\n".join(labels) + "\n")

    # Outlier highlights
    notes = [
        "DATASET_SYMBOL",
        "SEPARATOR TAB",
        "DATASET_LABEL\tOutliers",
        "COLOR\t#ff7f0e",
        "DATA",
    ]
    for leaf in t.get_leaves():
        if leaf.sra in HIGHLIGHTS:
            notes.append(
                f"{leaf.name}\t2\t6\t{HIGHLIGHTS[leaf.sra]}\t1\t1")
    (ITOL_DIR / "annotations_outliers.txt").write_text("\n".join(notes) + "\n")
    print(f"[iTOL] exported datasets to {ITOL_DIR}/")

    # --- Render rectangular ---
    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.mode = "r"
    ts.scale = 60
    ts.branch_vertical_margin = 1.2
    ts.show_scale = True
    ts.layout_fn = [layout_internal_bs, layout_leaf]

    for node in t.traverse():
        ns = NodeStyle()
        ns["size"] = 0
        ns["hz_line_width"] = 0.7
        ns["vt_line_width"] = 0.7
        # Colour by clade
        if node.is_leaf():
            ns["hz_line_color"] = COLOURS.get(node.grp, "#000000")
        else:
            leaves = node.get_leaves()
            grps = {l.grp for l in leaves}
            if len(grps) == 1:
                ns["hz_line_color"] = COLOURS.get(list(grps)[0], "#000000")
                ns["vt_line_color"] = COLOURS.get(list(grps)[0], "#000000")
        node.set_style(ns)

    rect_pdf = FIG_DIR / "phase4_tree_rect.pdf"
    rect_svg = FIG_DIR / "phase4_tree_rect.svg"
    ok = False
    for target in (rect_pdf, rect_svg):
        try:
            t.render(str(target), tree_style=ts, w=1200, units="px")
            print(f"[render] {target}")
            ok = True
            break
        except Exception as e:
            print(f"[fail] {target.name}: {e}")
    if not ok:
        # Last resort: plain ASCII
        (FIG_DIR / "phase4_tree_rect.txt").write_text(t.get_ascii(show_internal=False))
        print(f"[fallback] wrote ASCII tree")

    # --- Circular version ---
    ts_c = TreeStyle()
    ts_c.mode = "c"
    ts_c.show_leaf_name = False
    ts_c.layout_fn = [layout_leaf]
    ts_c.scale = 60
    ts_c.arc_start = -90
    ts_c.arc_span = 360

    for node in t.traverse():
        ns = NodeStyle()
        ns["size"] = 0
        ns["hz_line_width"] = 0.5
        ns["vt_line_width"] = 0.5
        if node.is_leaf():
            ns["hz_line_color"] = COLOURS.get(node.grp, "#000000")
        node.set_style(ns)

    circ_pdf = FIG_DIR / "phase4_tree_circ.pdf"
    circ_svg = FIG_DIR / "phase4_tree_circ.svg"
    ok = False
    for target in (circ_pdf, circ_svg):
        try:
            t.render(str(target), tree_style=ts_c, w=1500, units="px")
            print(f"[render] {target}")
            ok = True
            break
        except Exception as e:
            print(f"[fail] {target.name}: {e}")

    print("\n" + "=" * 70)
    print("Phase 4 complete.")
    print("=" * 70)


if __name__ == "__main__":
    main()
