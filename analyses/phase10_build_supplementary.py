"""
L4.13 Phase 10 -- Build supplementary materials
================================================

Generates:
    supplementary_materials/table_S1_strains.tsv       -- 258 strains + metadata
    supplementary_materials/table_S2_core_exclusive.tsv -- 249 SPDI + annotation
    supplementary_materials/table_S3_marker_validation.tsv
    supplementary_materials/table_S4_spoligotypes.tsv
    supplementary_materials/table_S5_mk_stratified.tsv
    supplementary_materials/supplementary.tex          -- compilable LaTeX
"""

import csv
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
RES = ROOT / "résultats"
SUPP = ROOT / "article" / "supplementary_materials"
SUPP.mkdir(parents=True, exist_ok=True)


def copy_as_tsv(src_csv, dst_tsv, description=""):
    """Convert CSV to TSV, preserving contents."""
    with src_csv.open() as f, dst_tsv.open("w", newline="") as g:
        reader = csv.reader(f)
        writer = csv.writer(g, delimiter="\t")
        for row in reader:
            writer.writerow(row)
    print(f"[write] {dst_tsv.name} ({description})")


def main():
    # Table S1: 258 strains + 10 outgroups with metadata
    src = DATA / "strains.csv"
    dst = SUPP / "table_S1_strains.tsv"
    # Read strain metadata
    with src.open() as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        fieldnames = reader.fieldnames[:]
    # Add 10 outgroup rows from strains.tsv (phase2_tree)
    strains_tsv = ROOT / "résultats" / "phase2_tree" / "strains.tsv"
    if strains_tsv.is_file():
        existing_ids = {r["strain_name"] for r in rows}
        with strains_tsv.open() as f:
            next(f)  # header
            for line in f:
                _lbl, sra, grp = line.rstrip().split("\t")
                if grp not in ("L4.4.2", "L4.2.1"):
                    continue
                if sra in existing_ids:
                    continue
                # Empty row with sub_lineage = outgroup clade
                new_row = {k: "" for k in fieldnames}
                new_row["strain_name"] = sra
                if "sub_lineage" in fieldnames:
                    new_row["sub_lineage"] = grp
                rows.append(new_row)
    with dst.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"[write] {dst.name} ({len(rows)} strains = 258 L4.13 + 10 outgroup)")

    # Table S2: core-exclusive SPDIs (249 total, annotated)
    src = DATA / "core_exclusive_annotated.csv"
    dst = SUPP / "table_S2_core_exclusive.tsv"
    copy_as_tsv(src, dst, "249 core-exclusive SPDIs + annotation")

    # Table S3: marker validation (inclusion + exclusion)
    src = RES / "phase1d_marker_validation.tsv"
    dst = SUPP / "table_S3_marker_validation.tsv"
    copy_as_tsv(src, dst, "3 marker validation (inclusion + exclusion)")

    # Table S4: spoligotypes
    src = DATA / "spoligotypes.csv"
    dst = SUPP / "table_S4_spoligotypes.tsv"
    copy_as_tsv(src, dst, "188 strains spoligotype octals")

    # Table S5: MK stratified (parse phase3e text into TSV)
    mk_text = (RES / "phase3e_mk_stratified.txt").read_text()
    dst = SUPP / "table_S5_mk_stratified.tsv"
    with dst.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["Stratum", "Dn", "Ds", "Pn", "Ps",
                    "Dn_Ds", "Pn_Ps", "Fisher_p_two",
                    "Fisher_p_greater", "DoS", "CI95_low", "CI95_high"])
        in_table = False
        for line in mk_text.split("\n"):
            if line.startswith("Stratum"):
                in_table = True
                continue
            if not in_table:
                continue
            # Skip separator dashes but don't break on them
            if line.startswith("---"):
                continue
            if line.startswith("Notes:"):
                break
            if not line.strip():
                continue
            parts = line.split()
            if len(parts) < 11:
                continue
            # Parse: Stratum Dn Ds Pn Ps Dn_Ds Pn_Ps pFisher pGreater DoS [CI95]
            try:
                stratum = parts[0]
                Dn, Ds, Pn, Ps = parts[1], parts[2], parts[3], parts[4]
                dn_ds = parts[5]
                pn_ps = parts[6]
                p_two = parts[7]
                p_gt = parts[8]
                dos = parts[9]
                ci = parts[10].strip("[]")
                ci_lo, ci_hi = ci.split(",")
                w.writerow([stratum, Dn, Ds, Pn, Ps, dn_ds, pn_ps,
                            p_two, p_gt, dos, ci_lo, ci_hi])
            except (ValueError, IndexError):
                continue
    print(f"[write] {dst.name} (MK stratified 6 strata)")

    # LaTeX supplementary document
    supp_tex = SUPP / "supplementary.tex"
    content = r"""\documentclass[11pt,a4paper]{article}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage[english]{babel}
\usepackage[margin=2cm]{geometry}
\usepackage{amsmath}
\usepackage{graphicx}
\usepackage{booktabs}
\usepackage{longtable}
\usepackage{xcolor}
\usepackage{hyperref}
\hypersetup{colorlinks=true, linkcolor=blue!60!black, urlcolor=blue!70!black}

\title{Supplementary Materials --\\
  Characterisation of L4.13, a cryptic European sub-lineage of
  \textit{Mycobacterium tuberculosis} with focal signatures of
  positive selection on ESX-2/4 secretion and PDIM biosynthesis}
\author{Christophe Guyeux}
\date{}

\begin{document}
\maketitle

\section*{Supplementary tables}

The following tables are provided as tab-separated files in the
\texttt{supplementary\_materials/} directory of the project
repository. Each is self-describing through its header row.

\begin{description}
\item[Table S1 --- Strains with metadata.]
File: \texttt{table\_S1\_strains.tsv}.
One row per strain (258 L4.13 + 10 outgroup), columns: strain
accession, sub-lineage assignment, SPDI count, QC metrics (missing
genes, missing RD, mean depth, covered bases percent, mean MAPQ,
mean BASEQ, GC content), and ENA metadata (country, collection
date, host, BioProject, BioSample, platform).

\item[Table S2 --- Core-exclusive SPDIs with annotation.]
File: \texttt{table\_S2\_core\_exclusive.tsv}. One row per SPDI
(249 total: 187 L4.13-global, 23 L4.13.1-specific, 39
L4.13.2-specific). Columns: SPDI identifier, clade (L4.13 /
L4.13.1 / L4.13.2), position, ref/alt, gene name, locus tag,
strand, predicted effect (\texttt{synonymous\_variant},
\texttt{missense\_variant}, \texttt{stop\_gained},
\texttt{frameshift\_variant}, \texttt{deletion},
\texttt{intergenic\_region}), impact severity, amino-acid change
where applicable, gene product, and Mycobrowser functional
category.

\item[Table S3 --- SNP marker validation.]
File: \texttt{table\_S3\_marker\_validation.tsv}. For each of the
three proposed markers (L4.13 \texttt{3222055:G:A}, L4.13.1
\texttt{187910:T:C}, L4.13.2 \texttt{1129237:G:A}), rows
reporting inclusion penetrance in the target clade
(L4.13.1 n=61, L4.13.2 n=197) and exclusion leakage in each of
17 tested sister sub-lineages (up to 50 strains per sister, seed
42, total n=1485).

\item[Table S4 --- In silico spoligotypes.]
File: \texttt{table\_S4\_spoligotypes.tsv}. Per-strain 43-spacer
binary pattern (columns) and 15-digit SITVIT2-style octal code,
for the 188 strains with complete \texttt{known\_coverage} data.
SITVIT2 returned zero matching isolates for all five
most-frequent octals of our dataset (see main text).

\item[Table S5 --- Stratified McDonald-Kreitman test.]
File: \texttt{table\_S5\_mk\_stratified.tsv}. Six functional
strata (global, PDIM strict, PDIM extended, ESX strict, ESX
extended, PE/PPE, rest-of-genome), with $D_{n}$, $D_{s}$,
$P_{n}$, $P_{s}$, $D_{n}/D_{s}$ : $P_{n}/P_{s}$, two-sided and
one-sided Fisher exact p-values, Direction of Selection (DoS)
with 10\,000-bootstrap 95\% confidence interval.
\end{description}

\section*{Supplementary figures}

\begin{figure}[!htb]
\centering
\includegraphics[width=0.85\textwidth]{../figures/phase5_root_to_tip.pdf}
\caption{\textbf{Supplementary Figure~S1 --- Root-to-tip regression.}
Scatter of root-to-L4.13-MRCA distance (in BIN+G
substitutions/site) vs sampling year for the 104 L4.13 tips with
known collection date. Blue: L4.13.1 (n=28 dated tips). Red:
L4.13.2 (n=76 dated tips). Linear regression yields $R^{2}=0.09$
and a negative slope, demonstrating that the BIN+G binary tree
carries insufficient temporal signal for reliable
tip-date-only molecular clock inference. This motivates the
ancient-sample-anchored LSD2 approach used in the main text.}
\end{figure}

\begin{figure}[!htb]
\centering
\includegraphics[width=\textwidth]{../figures/phase7_map_world.pdf}
\caption{\textbf{Supplementary Figure~S2 --- Global choropleth
distribution of L4.13.} Colour scale (\texttt{YlOrRd},
log-transformed) encodes the number of L4.13 strains per country
among the 108 strains with ENA-declared country of origin.
Countries in grey are either absent from our sample or have
\texttt{Unknown} / \texttt{not provided} metadata.}
\end{figure}

\begin{figure}[!htb]
\centering
\includegraphics[width=\textwidth]{../figures/phase7_map_composition.pdf}
\caption{\textbf{Supplementary Figure~S3 --- Global sub-clade
composition.} Pie charts on each country with $n\geq3$ showing
the proportion of L4.13.1 (blue) vs L4.13.2 (red). The United
States is the only region with substantial co-circulation of
both sub-clades; Georgia is exclusively L4.13.2.}
\end{figure}

\section*{Supplementary methods}

\paragraph{Ancient-sample anchors for LSD2 dating.}
The IQ-TREE~2.4.0 + LSD2 pipeline (main text, Molecular dating
section) used the following nine ancient calibrators retrieved
from the SPAAM ancient metagenome directory: three pinnipedii
genomes from Bos et al.\ 2014 (\texttt{SRR1238557},
\texttt{SRR1238558}, \texttt{SRR1238559}, $^{14}$C date
\texttt{b(1028,1280)}); three pinnipedii genomes from
V{\aa}gene~2022 (\texttt{S82} \texttt{b(1250,1470)},
\texttt{S281} \texttt{b(1265,1380)}, \texttt{S386}
\texttt{b(1260,1400)}); Winstrup's Lund~1 (Sabin~2020,
historical date 1679); and two Kay~2015 V{\'a}c body samples
(Body~68 cleaned for mixed infection via inter-L4 sharing,
Body~92, \texttt{b(1731,1838)} historical interval). The H37Rv
MRCA was constrained to \texttt{b(1900,1910)} as an additional
nodal anchor. The multi-seed workflow (5~seeds, objective
function selection) selected seed 335004 as the best run
(LSD2 objective = 0.432), yielding the point estimates quoted
in the main text.

\paragraph{Binary alignment and outgroup rationale.}
The 268-taxon binary alignment used 5\,049 parsimony-informative
sites (SPDIs present in $\geq$\,2 and absent in $\geq$\,2
strains). The outgroup was chosen on the basis of the global L4
topology (\texttt{lignees.py[moi]}) which places L4.13 as sister
to the L4.4-containing clade: 5 \texttt{L4.4.2} strains (direct
sister) plus 5 \texttt{L4.2.1} strains (distal L4 anchor) were
selected at random with seed 42.

\paragraph{Reproducibility.}
All analysis scripts (phases 1 through 9) are archived in the
\texttt{analyses/} directory of the project repository. Each
phase script is self-contained and reproducibly executes its
corresponding analysis from the raw TBannotator SPDI profiles
and ENA metadata.

\end{document}
"""
    supp_tex.write_text(content)
    print(f"[write] {supp_tex.name}")


if __name__ == "__main__":
    main()
