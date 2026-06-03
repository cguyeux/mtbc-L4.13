"""
L4.13 Phase 7 -- Prepare CSVs for /geo-map
==========================================

Builds two CSVs:

    data/geomap_choropleth.csv -- country, n  (total L4.13 strains)
    data/geomap_composition.csv -- country, lineage, n  (L4.13.1 vs L4.13.2)

These feed the /geo-map skill for two figures:
    figures/phase7_map_choropleth.pdf  -- total density per country
    figures/phase7_map_composition.pdf -- L4.13.1/L4.13.2 pies by country
"""

import csv
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
STRAINS_CSV = DATA / "strains.csv"


def clean_country(c):
    if not c: return None
    c = c.split(":")[0].strip()
    if not c or c.lower() in ("missing", "not provided", "not applicable",
                               "unknown", "none"):
        return None
    # Natural Earth normalisation
    fixes = {
        "USA": "United States of America",
        "United States": "United States of America",
        "UK": "United Kingdom",
        "Czech Republic": "Czechia",
        "South Korea": "Korea",
    }
    return fixes.get(c, c)


def main():
    by_country = Counter()
    by_country_lineage = defaultdict(Counter)

    with STRAINS_CSV.open() as f:
        for row in csv.DictReader(f):
            c = clean_country(row.get("country", ""))
            if c is None: continue
            sub = row.get("sub_lineage", "")
            by_country[c] += 1
            by_country_lineage[c][sub] += 1

    # Choropleth CSV
    out1 = DATA / "geomap_choropleth.csv"
    with out1.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["country", "n"])
        for country, n in sorted(by_country.items(), key=lambda x: -x[1]):
            w.writerow([country, n])
    print(f"[write] {out1}  ({len(by_country)} countries, "
          f"{sum(by_country.values())} strains)")

    # Composition CSV
    out2 = DATA / "geomap_composition.csv"
    with out2.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["country", "lineage", "n"])
        for country in sorted(by_country_lineage, key=lambda x: -by_country[x]):
            for sub, n in sorted(by_country_lineage[country].items()):
                w.writerow([country, sub, n])
    print(f"[write] {out2}  "
          f"({sum(len(v) for v in by_country_lineage.values())} rows)")

    print("\nTop 10 countries:")
    for c, n in sorted(by_country.items(), key=lambda x: -x[1])[:10]:
        comp = by_country_lineage[c]
        print(f"  {c:<30s} {n:>3d}  (L4.13.1={comp.get('L4.13.1',0)}, "
              f"L4.13.2={comp.get('L4.13.2',0)})")


if __name__ == "__main__":
    main()
