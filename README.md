# mtbc-L4.13

Code, derived data and figures accompanying:

> **Characterisation of L4.13, a cryptic European sub-lineage of *Mycobacterium tuberculosis* with candidate positive selection on accessory ESX-2/4 type VII secretion modules**
> Christophe Guyeux, Bassam AlKindy, Huda AlNayyef.

This repository lets readers inspect the analysis pipeline and re-run the
downstream statistics and figures of the paper from the provided derived data.

## Summary

L4.13 is a lineage-4 sub-lineage of the *M. tuberculosis* complex (MTBC),
defined by diagnostic markers in Freschi *et al.* (2021) but absent from the
Coll, Stucki and Napier barcodes and never individually characterised in the
post-2021 literature. From 258 whole-genome sequenced strains (61 L4.13.1,
197 L4.13.2) we propose three SNP markers validating the bipartite structure,
document an orphan spoligotype signature, reconstruct a Western-European
ancestral origin, date the L4 root, and detect a candidate focal signal of
positive selection on the accessory ESX-2/4 type VII secretion modules.

The three proposed SNP markers (SPDI, against H37Rv / NC_000962.3):

| Lineage | Marker | Inclusion | Leakage |
|---------|--------|-----------|---------|
| L4.13   | `NC_000962.3:3222055:G:A` | 258/258 | 0/1485 |
| L4.13.1 | `NC_000962.3:187910:T:C`  | 61/61   | 0/1485 |
| L4.13.2 | `NC_000962.3:1129237:G:A` | 197/197 | 0/1485 |

## Repository structure

```
mtbc-L4.13/
├── data/              # Derived data tables consumed by the analysis scripts
├── analyses/          # Phased analysis scripts (phase1*..phase11)
├── figures/           # Publication figures (PDF + PNG)
├── supplementary/     # Supplementary tables S1-S5, supplementary.tex/pdf
├── external/          # Placeholders for inputs not redistributed here (see external/README.md)
├── requirements.txt
├── CITATION.cff
└── LICENSE
```

## Reproducibility: three data tiers

The pipeline runs over three tiers. Only the last is redistributed here, which
is sufficient to reproduce every statistic and figure in the paper.

1. **Raw reads (public).** All 258 strains are public whole-genome sequencing
   datasets in the European Nucleotide Archive (ENA). Accessions, lineage
   assignment and per-strain quality metrics are in `data/strains.csv` (and
   `supplementary/table_S1_strains.tsv`).
2. **Per-strain SPDI profiles (not redistributed).** Reads are mapped against
   the H37Rv reference (NC_000962.3) to call variants in SPDI format
   (`NC_000962.3:pos:ref:alt`). The upstream scripts (`phase1*`, `phase2`,
   `phase3`, `phase3b`) consume one `spdi.txt` per strain; the directory holding
   them is configurable (see `external/README.md`). These profiles are not
   shipped because they are a large reprocessing of the public ENA reads.
3. **Derived data (provided in `data/`).** The downstream scripts (`phase3c`,
   `phase3e`, `phase5`, `phase6`, `phase7`, `phase8`, `phase4*`, `phase2b`,
   `phase10`) read only the CSV/TSV tables in `data/` and reproduce the
   statistics, tables and figures without needing tier 2.

## `data/` contents

| File | Description |
|------|-------------|
| `strains.csv` | The 258 strains: accession, sub-lineage, QC metrics, country, collection date, BioProject/BioSample, platform |
| `ena_metadata.csv` | Raw ENA metadata as retrieved |
| `classify_barcodings.csv` | Per-strain classification across 17 published MTBC schemes |
| `core_exclusive_all.csv`, `core_exclusive_annotated.csv`, `core_exclusive_L4.13*.txt` | Core-exclusive SPDI sets (L4.13 and the two sub-clades) and their functional annotation |
| `mk_input_raw.csv`, `mk_input_annotated.csv` | Tier 1 / Tier 4 input for the McDonald-Kreitman tests |
| `spoligotypes.csv` | In-silico 43-spacer spoligotypes and octal codes |
| `ancestral_metadata*.csv` | Region assignment for ancestral-state reconstruction |
| `geomap_choropleth.csv`, `geomap_composition.csv` | Inputs for the geographic maps |
| `thd_haplotypes.csv` | 258 x 4117 polymorphic-site matrix for Time-scaled Haplotypic Density |

## Running the downstream analyses

```bash
pip install -r requirements.txt
cd analyses
python phase3c_mk_test.py        # global McDonald-Kreitman test
python phase3e_mk_stratified.py  # pathway-stratified MK (ESX-2/4, PDIM, PE/PPE)
python phase5_root_to_tip_dating.py
python phase8_prepare_thd_input.py
# ...
```

These scripts resolve their inputs relative to the repository root
(`data/`) and write to a `résultats/` directory. The upstream scripts and the
convergence scripts additionally require the inputs documented in
`external/README.md`.

## Software

Python >= 3.10 with the packages in `requirements.txt`
(biopython, ete3, matplotlib, numpy, pandas, scipy).
External phylogenetic tools used by the upstream steps:
RAxML-NG v1.1.0 (https://github.com/amkozlov/raxml-ng) and
IQ-TREE 2.4.0 with the bundled LSD2 (http://www.iqtree.org/).

## Data availability

All sequencing data are public in the ENA; accessions are listed in
`data/strains.csv`. The spoligotype queries used SITVIT2
(http://www.pasteur-guadeloupe.fr:8081/SITVIT2/).

## Citation

See `CITATION.cff`. Please cite the paper (in preparation) when using this code or data.

## License

Code is released under the MIT License (`LICENSE`).
Derived data tables are released under CC-BY-4.0.
