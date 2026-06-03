#!/usr/bin/env python3
"""
L4.13 Analysis -- Phase 1c: Classification via all SNP barcodings
==================================================================

For each of the 258 L4.13 strains, classify the sample against every
published MTBC taxonomy scheme defined in the reference barcode marker
set (Coll, Freschi, Napier, Stucki, Coscolla, Shitikov23,
Palittapongarnpim, Shuaib, Lipworth, Netikul, Merker, Ates, Gisch,
Shitikov, Thawornwattana, Zwyer).

The "best match" per system is chosen as the deepest lineage code for
which all defining markers are satisfied (100% match, positive and
negative markers). If no full 100% match exists, we keep the deepest
code whose percentage is >= a threshold (default 50%).

Outputs:
    data/classify_barcodings.csv   -- wide table (1 row per strain)
    résultats/phase1c_discordances.txt -- discordance summary
    résultats/phase1c_summary.txt  -- system-by-system distribution

Usage:
    python phase1c_classify_barcodings.py
"""

import csv
import importlib.util
import sys
from collections import Counter, defaultdict
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent
REPO_ROOT = PROJECT_DIR.parent

BDD = Path(__file__).resolve().parent.parent / "external" / "spdi_db"  # per-strain SPDI profiles; see external/README.md
SUB_CLADES = ["L4.13.1", "L4.13.2"]
DATA_DIR = PROJECT_DIR / "data"
RESULTS_DIR = PROJECT_DIR / "résultats"

REFERENCE = "NC_000962.3"

# Import the classify helpers from the mtbc-lineages skill
SKILL_SCRIPT = Path(__file__).resolve().parent.parent / "external" / "lineages.py"
spec = importlib.util.spec_from_file_location("mtbc_lineages", SKILL_SCRIPT)
_mod = importlib.util.module_from_spec(spec)
sys.modules["mtbc_lineages"] = _mod
spec.loader.exec_module(_mod)
load_lignees = _mod.load_lignees
read_spdi_set = _mod.read_spdi_set


def classify_sample(sample_spdis, lignees, min_pct=50.0):
    """Return {system: best_code} for the sample.

    Best = deepest (most dotted components) code among lineages whose
    markers are >= min_pct satisfied. If a lineage has multiple markers,
    the percentage is the fraction satisfied (positive markers require
    PRESENCE, markers prefixed with '-' require ABSENCE).
    """
    result = {}
    for sys_name, entries in lignees.items():
        by_lineage = defaultdict(list)
        for code, spdi in entries:
            by_lineage[code.strip()].append(spdi)

        matches = []  # (code, pct, depth)
        for code, spdi_list in by_lineage.items():
            total = len(spdi_list)
            matched = 0
            for spdi in spdi_list:
                if spdi.startswith("-"):
                    if spdi[1:] not in sample_spdis:
                        matched += 1
                else:
                    if spdi in sample_spdis:
                        matched += 1
            if total == 0:
                continue
            pct = matched / total * 100.0
            if pct >= min_pct and matched > 0:
                depth = code.count(".") + 1
                matches.append((code, pct, depth, matched, total))

        if not matches:
            continue
        # Prefer 100% matches; among those, deepest code; fallback: highest pct.
        full = [m for m in matches if m[1] == 100.0]
        pool = full if full else matches
        pool.sort(key=lambda x: (-x[2], -x[1], x[0]))
        best_code, pct, _, matched, total = pool[0]
        if total == 1:
            result[sys_name] = best_code
        else:
            result[sys_name] = f"{best_code} ({matched}/{total})"
    return result


def list_valid_strains(directory):
    if not directory.is_dir():
        return []
    return sorted(
        entry.name
        for entry in directory.iterdir()
        if entry.is_dir() and (entry / REFERENCE / "spdi.txt").is_file()
    )


def main():
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("L4.13 Analysis -- Phase 1c: Barcoding classifications")
    print("=" * 70)

    print("\n[load] lignees.py")
    lignees, lignees_path, mtime, degraded = load_lignees()
    if not lignees:
        print("ERROR: could not load lignees.py")
        sys.exit(1)
    print(f"  path: {lignees_path}")
    print(f"  loaded: {mtime:%Y-%m-%d %H:%M}"
          + (" -- DEGRADED" if degraded else ""))
    print(f"  systems: {len(lignees)}")

    systems = ["moi"] + sorted(k for k in lignees.keys() if k != "moi")

    rows = []
    for sub in SUB_CLADES:
        sub_dir = BDD / sub
        strains = list_valid_strains(sub_dir)
        print(f"\n[scan] {sub}: {len(strains)} strains")
        for sra in strains:
            spdi_file = sub_dir / sra / REFERENCE / "spdi.txt"
            sample = read_spdi_set(str(spdi_file))
            cls = classify_sample(sample, lignees, min_pct=50.0)
            row = {"strain_name": sra, "sub_lineage_dir": sub, "n_spdis": len(sample)}
            for s in systems:
                row[s] = cls.get(s, "")
            rows.append(row)

    # Write CSV
    fieldnames = ["strain_name", "sub_lineage_dir", "n_spdis"] + systems
    csv_path = DATA_DIR / "classify_barcodings.csv"
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        for row in sorted(rows,
                          key=lambda r: (r["sub_lineage_dir"], r["strain_name"])):
            w.writerow(row)
    print(f"\n[write] {csv_path} ({len(rows)} rows, {len(systems)} systems)")

    # System-by-system distribution
    summary_lines = ["L4.13 -- Phase 1c: barcoding classifications\n", "=" * 50 + "\n\n"]
    summary_lines.append(f"Total strains: {len(rows)}\n")
    summary_lines.append(f"Systems: {len(systems)}\n")
    summary_lines.append(f"lignees.py loaded: {mtime:%Y-%m-%d %H:%M}\n\n")
    for s in systems:
        summary_lines.append(f"--- {s} ---\n")
        c = Counter(r[s] or "(no match)" for r in rows)
        for label, n in c.most_common():
            summary_lines.append(f"  {label or '(no match)'}: {n}\n")
        summary_lines.append("\n")

    # Discordance analysis vs physical ranging
    print("\n[discordance analysis]")
    # For each physical sub_lineage, count the "moi" assignments
    moi_vs_dir = defaultdict(Counter)
    for r in rows:
        moi_vs_dir[r["sub_lineage_dir"]][r["moi"] or "(no match)"] += 1
    disc_lines = ["L4.13 -- Discordance analysis (barcode vs placement)\n",
                  "=" * 50 + "\n\n"]
    for sub in SUB_CLADES:
        disc_lines.append(f"Physical directory {sub} ({sum(moi_vs_dir[sub].values())} strains):\n")
        for label, n in moi_vs_dir[sub].most_common():
            disc_lines.append(f"  moi = {label:<25s}  {n}\n")
        disc_lines.append("\n")
        print(f"  {sub} ({sum(moi_vs_dir[sub].values())} strains):")
        for label, n in moi_vs_dir[sub].most_common():
            print(f"    moi = {label:<25s}  {n}")

    # Strains where moi != expected
    disc_lines.append("Strains with 'moi' NOT matching the directory:\n")
    n_conflicts = 0
    for r in rows:
        expected = r["sub_lineage_dir"].lstrip("L")  # "4.13.1" or "4.13.2"
        moi = (r["moi"] or "").strip()
        # consider conflict if moi is defined but doesn't start with "4.13"
        # or starts with 4.13 but differs on last component
        if not moi:
            continue
        moi_code = moi.split()[0]
        if not moi_code.startswith("4.13"):
            disc_lines.append(f"  {r['strain_name']:<15s} dir={r['sub_lineage_dir']:<8s} "
                              f"moi={moi_code}\n")
            n_conflicts += 1
        elif moi_code != expected and moi_code != "4.13":
            disc_lines.append(f"  {r['strain_name']:<15s} dir={r['sub_lineage_dir']:<8s} "
                              f"moi={moi_code}\n")
            n_conflicts += 1
    disc_lines.insert(0, f"Conflicts (strain in dir X but moi assigns Y != X): {n_conflicts}\n\n")
    print(f"\n  Conflicts (dir vs moi): {n_conflicts}")

    with open(RESULTS_DIR / "phase1c_discordances.txt", "w") as f:
        f.writelines(disc_lines)
    with open(RESULTS_DIR / "phase1c_summary.txt", "w") as f:
        f.writelines(summary_lines)
    print(f"[write] {RESULTS_DIR / 'phase1c_discordances.txt'}")
    print(f"[write] {RESULTS_DIR / 'phase1c_summary.txt'}")

    print("\n" + "=" * 70)
    print("Phase 1c complete.")
    print("=" * 70)


if __name__ == "__main__":
    main()
