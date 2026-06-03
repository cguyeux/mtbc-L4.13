# Supplementary materials — L4.13

This directory is **self-contained**: every file referenced by the
supplementary document (`supplementary.tex` / `supplementary.pdf`)
and by the main manuscript is available here.

## Contents

| File | Purpose |
|------|---------|
| `supplementary.tex` | LaTeX source of the supplementary document |
| `supplementary.pdf` | Compiled supplementary document (3 pages) |
| `table_S1_strains.tsv` | 258 L4.13 strains + 10 outgroup — metadata (QC + ENA) |
| `table_S2_core_exclusive.tsv` | 249 core-exclusive SPDIs (187 L4.13 + 23 L4.13.1 + 39 L4.13.2) with functional annotation |
| `table_S3_marker_validation.tsv` | Inclusion + exclusion tests of the three L4.13 SNP markers across 17 sister sub-lineages (n = 1 485) |
| `table_S4_spoligotypes.tsv` | In silico 43-spacer spoligotypes of the 188 strains with full coverage data |
| `table_S5_mk_stratified.tsv` | Stratified McDonald-Kreitman table (7 strata, Dn/Ds/Pn/Ps/DoS/p) |
| `figures/phase5_root_to_tip.pdf` | Supplementary Figure S1 — root-to-tip regression |
| `figures/phase7_map_world.pdf` | Supplementary Figure S2 — global choropleth |
| `figures/phase7_map_composition.pdf` | Supplementary Figure S3 — global sub-clade composition |

## Compilation

```bash
pdflatex supplementary.tex
```

No BibTeX pass is required (no citations in the supplementary
document).

## Citation

If used, please cite the main manuscript (Guyeux C., *Characterisation
of L4.13*, in preparation).
