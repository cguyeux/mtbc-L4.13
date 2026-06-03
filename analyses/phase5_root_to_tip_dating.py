"""
L4.13 Phase 5 -- Root-to-tip regression dating (TempEst-style)
===============================================================

Non-Bayesian point-estimate dating:
  1. Load the RAxML support tree, re-rooted on L4.4.2 MRCA.
  2. For each L4.13 leaf with a known collection_date, compute distance
     to the L4.13 MRCA (in substitutions per site).
  3. Linear regression: distance = rate * (date - t_root).
  4. Report clock rate, TMRCA, and bifurcation date L4.13.1 / L4.13.2.
  5. Quality checks: R², Pearson r, residuals, outlier tips.

Also estimates independently the TMRCA of each sub-clade.

Input:
    résultats/phase2_tree/L4.13.raxml.support
    résultats/phase2_tree/strains.tsv
    data/strains.csv (with collection_date column)

Output:
    résultats/phase5_root_to_tip_regression.txt
    article/figures/phase5_root_to_tip.pdf

Usage:
    python phase5_root_to_tip_dating.py
"""

from __future__ import annotations

import csv
import re
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from Bio import Phylo

ROOT = Path(__file__).resolve().parent.parent
TREE_DIR = ROOT / "résultats" / "phase2_tree"
DATA_DIR = ROOT / "data"
RES_DIR = ROOT / "résultats"
FIG_DIR = ROOT / "article" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

SUPPORT_TREE = TREE_DIR / "L4.13.raxml.support"
STRAINS_TSV = TREE_DIR / "strains.tsv"
STRAINS_CSV = DATA_DIR / "strains.csv"


def parse_date(s):
    """Return float year (e.g. 2018.45) or None."""
    if not s:
        return None
    s = s.strip()
    # YYYY
    m = re.match(r"^(\d{4})$", s)
    if m:
        return float(m.group(1))
    # YYYY-MM
    m = re.match(r"^(\d{4})-(\d{1,2})$", s)
    if m:
        y, mo = int(m.group(1)), int(m.group(2))
        return y + (mo - 0.5) / 12.0
    # YYYY-MM-DD
    m = re.match(r"^(\d{4})-(\d{1,2})-(\d{1,2})$", s)
    if m:
        y, mo, d = int(m.group(1)), int(m.group(2)), int(m.group(3))
        return y + (mo - 1 + (d - 1) / 30.0) / 12.0
    return None


def load_groups():
    g = {}
    with STRAINS_TSV.open() as f:
        next(f)
        for line in f:
            lbl, sra, grp = line.rstrip().split("\t")
            g[lbl] = (grp, sra)
    return g


def load_dates():
    """Return sra -> float year."""
    d = {}
    if not STRAINS_CSV.is_file():
        return d
    with STRAINS_CSV.open() as f:
        for row in csv.DictReader(f):
            y = parse_date(row.get("collection_date", ""))
            if y is not None:
                d[row["strain_name"]] = y
    return d


def root_tree(tree, groups):
    l442 = [l for l, (g, _) in groups.items() if g == "L4.4.2"]
    tree.root_with_outgroup(*l442)


def compute_root_distances(tree, target_mrca):
    """Return dict leaf -> branch distance to target_mrca."""
    dist = {}
    for leaf in tree.get_terminals():
        # Biopython Phylo .distance(a, b): branch length between two clades
        try:
            dist[leaf.name] = tree.distance(target_mrca, leaf)
        except Exception:
            pass
    return dist


def main():
    print("=" * 70)
    print("L4.13 Phase 5 -- root-to-tip regression dating")
    print("=" * 70)

    groups = load_groups()
    dates = load_dates()
    tree = Phylo.read(str(SUPPORT_TREE), "newick")
    root_tree(tree, groups)

    # L4.13 ingroup leaves
    l413_labels = [l for l, (g, _) in groups.items()
                   if g in ("L4.13.1", "L4.13.2")]
    l413_tips = [t for t in tree.get_terminals() if t.name in l413_labels]
    l413_mrca = tree.common_ancestor(l413_tips)

    # Sub-clade MRCAs
    l131_labels = [l for l, (g, _) in groups.items() if g == "L4.13.1"]
    l132_labels = [l for l, (g, _) in groups.items() if g == "L4.13.2"]
    l131_mrca = tree.common_ancestor([t for t in tree.get_terminals()
                                      if t.name in l131_labels])
    l132_mrca = tree.common_ancestor([t for t in tree.get_terminals()
                                      if t.name in l132_labels])

    dist_to_l413 = compute_root_distances(tree, l413_mrca)

    # Build regression dataset: leaves with date
    rows = []
    for leaf in l413_tips:
        grp, sra = groups[leaf.name]
        y = dates.get(sra)
        d = dist_to_l413.get(leaf.name)
        if y is None or d is None:
            continue
        rows.append((leaf.name, sra, grp, y, d))

    print(f"\n[data] {len(rows)}/{len(l413_tips)} L4.13 tips have "
          f"collection_date")
    if len(rows) < 20:
        print("[warn] too few dated tips for a reliable regression")

    years = np.array([r[3] for r in rows])
    dists = np.array([r[4] for r in rows])

    # Linear regression
    slope, intercept = np.polyfit(years, dists, 1)
    # R^2
    pred = slope * years + intercept
    ss_res = np.sum((dists - pred) ** 2)
    ss_tot = np.sum((dists - dists.mean()) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    # Pearson r
    r = np.corrcoef(years, dists)[0, 1]

    # Estimate TMRCA (year where distance would be 0 at the MRCA)
    # pred = slope * year + intercept -> year0 = -intercept/slope
    tmrca_l413 = -intercept / slope if slope != 0 else float("nan")

    # Clock rate (substitutions/site/year)
    rate = slope

    # Sub-clade MRCAs: use the same rate to back-date
    # Distance from l413_mrca to each sub-clade MRCA
    d_to_131 = tree.distance(l413_mrca, l131_mrca)
    d_to_132 = tree.distance(l413_mrca, l132_mrca)

    # If rate > 0 and TMRCA estimated, sub-clade TMRCAs:
    tmrca_l131 = tmrca_l413 + (d_to_131 / rate) if rate > 0 else float("nan")
    tmrca_l132 = tmrca_l413 + (d_to_132 / rate) if rate > 0 else float("nan")

    # Also mean sampling year
    mean_year = years.mean()

    # Bootstrap CI on rate and TMRCA (resample dated tips)
    rng = np.random.default_rng(42)
    n_boot = 1000
    rates_boot = []
    tmrca_boot = []
    for _ in range(n_boot):
        idx = rng.integers(0, len(rows), len(rows))
        yb, db = years[idx], dists[idx]
        try:
            s_b, i_b = np.polyfit(yb, db, 1)
            if s_b > 0:
                rates_boot.append(s_b)
                tmrca_boot.append(-i_b / s_b)
        except Exception:
            pass
    if rates_boot:
        rate_lo, rate_hi = np.percentile(rates_boot, [2.5, 97.5])
        tmrca_lo, tmrca_hi = np.percentile(tmrca_boot, [2.5, 97.5])
    else:
        rate_lo = rate_hi = tmrca_lo = tmrca_hi = float("nan")

    # --- Report ---
    lines = []
    lines.append("L4.13 -- Root-to-tip regression dating\n")
    lines.append("=" * 50 + "\n\n")
    lines.append(f"Dated L4.13 tips used: {len(rows)} / {len(l413_tips)} "
                 f"({100*len(rows)/len(l413_tips):.0f}%)\n")
    lines.append(f"Year range: {years.min():.0f}  -  {years.max():.0f}  "
                 f"(mean {mean_year:.1f})\n\n")

    lines.append("Regression statistics\n")
    lines.append("-" * 40 + "\n")
    lines.append(f"  Slope (clock rate): {slope:.3e} subst/site/year\n")
    lines.append(f"  Intercept: {intercept:.5f}\n")
    lines.append(f"  R^2 = {r2:.4f},   Pearson r = {r:.4f}\n")
    lines.append(f"  Bootstrap 95% CI on rate: "
                 f"[{rate_lo:.3e}, {rate_hi:.3e}]\n\n")

    lines.append("TMRCA estimates\n")
    lines.append("-" * 40 + "\n")
    lines.append(f"  L4.13 MRCA: {tmrca_l413:.1f} "
                 f"(95% CI [{tmrca_lo:.1f}, {tmrca_hi:.1f}])\n")
    lines.append(f"  L4.13.1 MRCA (branch from L4.13 = {d_to_131:.5f}): "
                 f"{tmrca_l131:.1f}\n")
    lines.append(f"  L4.13.2 MRCA (branch from L4.13 = {d_to_132:.5f}): "
                 f"{tmrca_l132:.1f}\n\n")

    # Reference MTBC rate for comparison
    lines.append("Reference MTBC clock rates (literature)\n")
    lines.append("-" * 40 + "\n")
    lines.append("  Duchene et al. 2016 (shallow MTBC tree): ~5e-8\n")
    lines.append("  Menardo et al. 2019 (range):             0.3-0.7 x 1e-7\n")
    lines.append("  Gagneux lab typical:                     5e-8 -- 1e-7\n")
    lines.append(f"  -- observed here:                       {slope:.2e}\n\n")

    if slope < 1e-8 or slope > 2e-7:
        lines.append("WARN: estimated rate is outside typical MTBC range.\n")
        lines.append("  Possible causes: weak temporal signal, dating\n")
        lines.append("  uncertainty, recent-biased sampling. Report with\n")
        lines.append("  caution, preferably cross-check with LSD2 or BEAST2.\n\n")

    # Outlier tips (residuals > 3 MAD)
    residuals = dists - pred
    mad = np.median(np.abs(residuals - np.median(residuals))) * 1.4826
    thresh = 3 * mad if mad > 0 else 3 * residuals.std()
    outliers = [rows[i] for i in range(len(rows))
                if abs(residuals[i]) > thresh]
    lines.append(f"Outliers (|residual| > {thresh:.5f} = 3*MAD)\n")
    lines.append("-" * 40 + "\n")
    lines.append(f"  n = {len(outliers)}\n")
    for label, sra, grp, y, d in outliers[:10]:
        lines.append(f"    {sra:<14s} [{grp}] year={y:.1f} "
                     f"dist={d:.5f}\n")
    if len(outliers) > 10:
        lines.append(f"    ... (+{len(outliers)-10} more)\n")

    txt = "".join(lines)
    print(txt)
    (RES_DIR / "phase5_root_to_tip_regression.txt").write_text(txt)

    # --- Figure ---
    fig, ax = plt.subplots(figsize=(8, 6))
    colours = {"L4.13.1": "#1f77b4", "L4.13.2": "#d62728"}
    for grp in ("L4.13.1", "L4.13.2"):
        xs = [r[3] for r in rows if r[2] == grp]
        ys = [r[4] for r in rows if r[2] == grp]
        ax.scatter(xs, ys, s=18, color=colours[grp], alpha=0.65,
                   edgecolor="none", label=f"{grp} (n={len(xs)})")

    xline = np.linspace(min(years) - 1, max(years) + 1, 100)
    ax.plot(xline, slope * xline + intercept,
            color="black", lw=1.2,
            label=f"y = {slope:.2e}·year + {intercept:.4f}")
    ax.set_xlabel("Sampling year", fontsize=10)
    ax.set_ylabel("Root-to-L4.13-MRCA distance\n(subst/site)", fontsize=10)
    ax.set_title(
        f"L4.13 root-to-tip regression (n={len(rows)}, R²={r2:.3f})\n"
        f"Clock rate = {slope:.2e} /site/year, TMRCA ≈ {tmrca_l413:.0f} "
        f"[{tmrca_lo:.0f}-{tmrca_hi:.0f}]",
        fontsize=10, pad=8)
    ax.legend(fontsize=8, loc="best")
    ax.grid(alpha=0.2)
    plt.tight_layout()
    out_pdf = FIG_DIR / "phase5_root_to_tip.pdf"
    out_png = FIG_DIR / "phase5_root_to_tip.png"
    plt.savefig(out_pdf, dpi=600, bbox_inches="tight")
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    print(f"\n[write] {out_pdf}")
    print(f"[write] {out_png}")

    print("\n" + "=" * 70)
    print("Phase 5 complete.")
    print("=" * 70)


if __name__ == "__main__":
    main()
