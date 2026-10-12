# CRC-Ferroptosis-MR-scRNAseq

## Overview

This repository contains R scripts supporting the study **“Integrative Mendelian Randomization and Single-Cell Transcriptomic Analyses Reveal Ferroptosis-Related Genes Associated with Colorectal Cancer Risk and Tumor Microenvironmental Remodeling.

The analyses integrate bulk transcriptomics, Mendelian randomization, Bayesian colocalization, proteomics, immune-signature analysis, and single-cell RNA sequencing to investigate ferroptosis-annotated candidate genes in colorectal cancer.

## Data Sources

- GSE44861: Discovery cohort and diagnostic model development
- GSE44076 and GSE25070: Independent validation cohorts
- GSE294300: Single-cell RNA-sequencing analysis
- eQTLGen: cis-eQTL exposure data
- IEU OpenGWAS: Colorectal cancer GWAS (ebi-a-GCST90018808)
- CPTAC/PDC000116: Proteomic assessment
- FerrDb V2: Ferroptosis-related gene annotation

## Analysis Scripts

The `scripts/` directory contains 12 R scripts covering:

1. Expression data preprocessing
2. Differential expression analysis
3. Mendelian randomization
4. Bayesian colocalization
5. Diagnostic model construction
6. Expression validation
7. Hallmark GSEA
8. Immune-signature analysis
9. Predicted ceRNA network construction
10. Single-cell RNA-seq analysis
11. Patient-level pseudobulk analysis
12. Epithelial-state and pseudotime analysis

## Software

Analyses were performed using R version 4.5.2. Individual scripts specify their required R packages.

## Reproducibility

The scripts are organized according to the main analytical workflow. Users should obtain the publicly available input datasets from the repositories listed above and configure the input and output paths before execution.

## Data Availability

Key analytical results are reported in the manuscript and supplementary tables. Additional processed data may be requested from the corresponding author.

## License

This repository is distributed under the MIT License.
