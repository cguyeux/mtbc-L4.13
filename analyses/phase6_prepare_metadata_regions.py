"""
L4.13 Phase 6 -- Prepare metadata for ancestral reconstruction
===============================================================

Builds a CSV compatible with the /ancestral-reconstruction script:
tree labels are '<subclade>_<sra>' (as stored in strains.tsv), each
paired with a country and a broader region.

Region mapping (11 regions) chosen to balance resolution vs power:
    North America, Latin America, Western Europe, Eastern Europe,
    Caucasus, Middle East, East Asia, South Asia, Southeast Asia,
    Sub-Saharan Africa, Oceania, Unknown

Output:
    data/ancestral_metadata.csv (strain_id, country, region)
"""

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
TREE_DIR = ROOT / "résultats" / "phase2_tree"

STRAINS_CSV = DATA / "strains.csv"
STRAINS_TSV = TREE_DIR / "strains.tsv"
OUT = DATA / "ancestral_metadata.csv"

# Country -> region mapping
COUNTRY_REGION = {
    "USA": "North America", "Canada": "North America", "Mexico": "North America",
    "Brazil": "Latin America", "Peru": "Latin America", "Colombia": "Latin America",
    "Argentina": "Latin America", "Chile": "Latin America",
    "United Kingdom": "Western Europe", "UK": "Western Europe",
    "France": "Western Europe", "Germany": "Western Europe",
    "Netherlands": "Western Europe", "Belgium": "Western Europe",
    "Switzerland": "Western Europe", "Italy": "Western Europe",
    "Spain": "Western Europe", "Portugal": "Western Europe",
    "Ireland": "Western Europe", "Sweden": "Western Europe",
    "Norway": "Western Europe", "Denmark": "Western Europe",
    "Finland": "Western Europe", "Austria": "Western Europe",
    "Romania": "Eastern Europe", "Bulgaria": "Eastern Europe",
    "Poland": "Eastern Europe", "Czech": "Eastern Europe",
    "Slovakia": "Eastern Europe", "Hungary": "Eastern Europe",
    "Serbia": "Eastern Europe", "Kosovo": "Eastern Europe",
    "Croatia": "Eastern Europe", "Russia": "Eastern Europe",
    "Ukraine": "Eastern Europe", "Belarus": "Eastern Europe",
    "Moldova": "Eastern Europe", "Albania": "Eastern Europe",
    "Macedonia": "Eastern Europe",
    "Georgia": "Caucasus", "Armenia": "Caucasus", "Azerbaijan": "Caucasus",
    "Turkey": "Middle East", "Lebanon": "Middle East", "Israel": "Middle East",
    "Iran": "Middle East", "Iraq": "Middle East", "Jordan": "Middle East",
    "Saudi Arabia": "Middle East", "Egypt": "Middle East", "Syria": "Middle East",
    "China": "East Asia", "Japan": "East Asia", "Korea": "East Asia",
    "South Korea": "East Asia", "Taiwan": "East Asia", "Mongolia": "East Asia",
    "India": "South Asia", "Pakistan": "South Asia", "Bangladesh": "South Asia",
    "Sri Lanka": "South Asia", "Nepal": "South Asia",
    "Vietnam": "Southeast Asia", "Thailand": "Southeast Asia",
    "Indonesia": "Southeast Asia", "Philippines": "Southeast Asia",
    "Malaysia": "Southeast Asia", "Myanmar": "Southeast Asia",
    "Cambodia": "Southeast Asia", "Laos": "Southeast Asia",
    "Singapore": "Southeast Asia",
    "Tanzania": "Sub-Saharan Africa", "Uganda": "Sub-Saharan Africa",
    "Kenya": "Sub-Saharan Africa", "Ethiopia": "Sub-Saharan Africa",
    "Ghana": "Sub-Saharan Africa", "Nigeria": "Sub-Saharan Africa",
    "South Africa": "Sub-Saharan Africa", "Malawi": "Sub-Saharan Africa",
    "Zambia": "Sub-Saharan Africa", "Zimbabwe": "Sub-Saharan Africa",
    "Cameroon": "Sub-Saharan Africa", "Ivory Coast": "Sub-Saharan Africa",
    "Cote d'Ivoire": "Sub-Saharan Africa",
    "Australia": "Oceania", "New Zealand": "Oceania",
}


def to_region(country_raw):
    if not country_raw or country_raw.strip().lower() in ("", "missing",
                                                           "not provided",
                                                           "not applicable",
                                                           "unknown"):
        return "Unknown"
    # Strip after colon (e.g. 'United Kingdom: Midlands' -> 'United Kingdom')
    c = country_raw.split(":")[0].strip()
    return COUNTRY_REGION.get(c, "Unknown_" + c)


def main():
    # Load tree labels
    labels = {}  # sra -> label
    with STRAINS_TSV.open() as f:
        next(f)
        for line in f:
            lbl, sra, grp = line.rstrip().split("\t")
            labels[sra] = lbl

    # Load country per SRA
    sra_country = {}
    with STRAINS_CSV.open() as f:
        for row in csv.DictReader(f):
            c = row.get("country", "")
            sra_country[row["strain_name"]] = c

    rows = []
    region_counts = {}
    for sra, lbl in labels.items():
        c = sra_country.get(sra, "")
        region = to_region(c)
        # Clean country (same as for region)
        clean_country = c.split(":")[0].strip() if c else ""
        rows.append({
            "strain_id": lbl,
            "country": clean_country or "Unknown",
            "region": region if not region.startswith("Unknown_") else "Unknown",
        })
        region_counts[region] = region_counts.get(region, 0) + 1

    # Write
    with OUT.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["strain_id", "country", "region"])
        w.writeheader()
        for row in rows:
            w.writerow(row)

    print(f"[write] {OUT} ({len(rows)} rows)")
    print(f"\nRegion distribution:")
    for r, n in sorted(region_counts.items(), key=lambda x: -x[1]):
        print(f"  {r:<25s} {n}")

    # Report unknown granular countries
    unknown_countries = set()
    for row in rows:
        r = to_region(row["country"])
        if r.startswith("Unknown_"):
            unknown_countries.add(r[8:])
    if unknown_countries:
        print(f"\n[warn] countries not in mapping: {sorted(unknown_countries)}")


if __name__ == "__main__":
    main()
