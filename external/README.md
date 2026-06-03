# external/ inputs

This directory holds inputs that the analysis scripts reference but that are
**not redistributed** in this repository, either because they are large
reprocessings of public data, public reference files better obtained from their
authoritative source, or derived data from companion studies. Place the
corresponding files here (or set the path in the relevant script) to run the
upstream and convergence steps. None of these are needed to reproduce the
downstream statistics and figures from `data/`.

| Expected path | Used by | What it is / where to get it |
|---------------|---------|------------------------------|
| `external/spdi_db/L4.13.1/<SRA>/NC_000962.3/spdi.txt`, `.../L4.13.2/...` | `phase1*`, `phase2`, `phase3`, `phase3b` | Per-strain SPDI variant calls, one `spdi.txt` per strain, obtained by mapping the ENA reads (accessions in `data/strains.csv`) against H37Rv (NC_000962.3). |
| `external/NC_000962.3.gff3`, `external/NC_000962.3.gb` | `phase3b` | H37Rv reference annotation (GFF3 and GenBank). Public; download from NCBI assembly `GCF_000195955.2` / RefSeq `NC_000962.3`. |
| `raxml-ng` (on `PATH`) | `phase2` | RAxML-NG v1.1.0 phylogenetic inference. https://github.com/amkozlov/raxml-ng |
| `external/lineages.py` | `phase1c` | Helper exposing `load_lignees()`, returning the reference SNP-barcode marker definitions for the published MTBC schemes used in multi-system classification. |
| `external/annotate_spdis.py` | `phase3b` | Helper annotating SPDI variants (gene, effect, amino-acid change) against the H37Rv GFF3/GenBank. |
| `external/geo_map.py` | `phase11` | Plotting helper for the geographic maps (final maps are already provided in `../figures/`). |
| `external/convergence_reference/*.tsv` | `phase3d` | List of ~1986 genes convergently mutated in the four animal-adapted MTBC ecotypes, from a companion study (animal vs human comparative genomics). |
| `external/animal_ecotypes_annotated/annotated_<Species>.tsv` | `phase3f` | Per-species annotated SPDI tables for the four animal ecotypes (Suricattae, Dassie bacillus, M. mungi, chimpanzee bacillus), from the same companion study. |
| `external/L4.15_core_exclusive_annotated.csv` | `phase3d` | Core-exclusive annotated SPDI set of the sister human sub-lineage L4.15, from a companion study. |

The convergence-reference and animal-ecotype tables come from companion work
that is being published separately; they are referenced in the paper as a
curated convergence list rather than redistributed here.
