# Transcriptomic Stability During Long-Term Cell Culture

[![R](https://img.shields.io/badge/R-4.5+-276DC3?logo=r)](https://www.r-project.org/)
[![Bioconductor](https://img.shields.io/badge/Bioconductor-DESeq2-green)](https://bioconductor.org/packages/DESeq2/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## Overview

Long-term in vitro culture can alter gene expression and affect the reproducibility of biological experiments. This project reproduces the downstream RNA-seq differential expression analysis of the publicly available **GSE268437** dataset, investigating transcriptomic changes associated with prolonged culture across multiple human cancer cell lines.

The workflow was implemented in **R/Bioconductor** using **DESeq2**, including differential expression analysis, principal component analysis (PCA), volcano plot visualization, and clustered heatmap generation.

This project was completed as part of my bioinformatics portfolio while preparing for PhD applications in computational genomics.

---

## Objectives

- Reproduce a published RNA-seq downstream analysis.
- Explore transcriptomic alterations after long-term cell culture.
- Perform differential expression analysis using DESeq2.
- Visualize global transcriptomic variation using PCA, Volcano Plot and Heatmap.
- Build a reproducible RNA-seq analysis workflow.

---

## Dataset

**Accession:** GSE268437

**Source:** NCBI Gene Expression Omnibus (GEO)

**Species:** Homo sapiens

**Sequencing platform:** Illumina NextSeq 500

**Input data:** Processed count matrix provided by the authors

---

## Analysis Workflow

```
GEO Dataset
      │
      ▼
Processed Count Matrix
      │
      ▼
Metadata Construction
      │
      ▼
DESeq2
      │
      ▼
Differential Expression Analysis
      │
      ├────────► PCA
      │
      ├────────► Volcano Plot
      │
      └────────► Heatmap
```

---

## Repository Structure

```
transcriptomics-stability-rnaseq/

├── data/
│   └── processed/
│
├── figures/
│   ├── volcano_m12_vs_m0.png
│   ├── pca_cellline.png
│   ├── pca_time.png
│   └── heatmap_top50_genes.png
│
├── results/
│   └── DEGs_m12_vs_m0.csv
│
├── scripts/
│
├── notebooks/
│
├── README.md
└── LICENSE
```

---

## Results

### Principal Component Analysis

![PCA Cell Line](figures/pca_cellline.png)

Cell line identity explains the largest proportion of transcriptomic variance across samples.

---

### Temporal Transcriptomic Variation

![PCA Time](figures/pca_time.png)

Samples show progressive transcriptomic changes following prolonged culture.

---

### Differential Gene Expression

![Volcano Plot](figures/volcano_m12_vs_m0.png)

Differential expression analysis identified multiple genes significantly altered after long-term culture, including genes involved in DNA repair, RNA processing and cell-cycle regulation.

---

### Heatmap

![Heatmap](figures/heatmap_top50_genes.png)

Hierarchical clustering demonstrates consistent transcriptional differences between culture time points.

---

## Software

- R
- Bioconductor
- DESeq2
- EnhancedVolcano
- pheatmap
- readxl
- ggplot2
- WSL (Linux)
- Conda
- STAR (environment setup)
- GENCODE reference genome

---

## Reproducibility

The complete downstream analysis can be reproduced directly from the processed count matrix included in this repository.

---

## Limitations

This repository focuses on **downstream RNA-seq analysis** using the processed count matrix supplied by the original study. Raw FASTQ preprocessing, alignment, and read quantification were not reproduced as part of this project.

---

## Future Improvements

- Functional enrichment analysis (GO / KEGG)
- Gene Set Enrichment Analysis (GSEA)
- Sample distance analysis
- Batch effect exploration
- Reproduction of the complete RNA-seq pipeline from raw FASTQ files

---

## Citation

Dataset:

> GSE268437 — NCBI Gene Expression Omnibus

Analysis:

> Love MI, Huber W, Anders S. *Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.* Genome Biology (2014).