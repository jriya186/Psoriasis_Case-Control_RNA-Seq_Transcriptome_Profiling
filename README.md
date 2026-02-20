# Psoriasis RNA-seq Analysis: Oxidative Stress Transcriptional Signature

A bulk RNA-seq analysis pipeline examining differential gene expression between lesional psoriatic and normal skin, with a focused investigation into the oxidative stress transcriptional landscape in psoriasis.


## Background

Psoriasis is a chronic inflammatory skin disease affecting approximately 2-3% of the world's population. While the inflammatory transcriptional signature of psoriasis is well established, the role of oxidative stress, an imbalance between reactive oxygen species (ROS) production and antioxidant defense, remains an important and underexplored contributor to disease pathogenesis.

This project uses a large publicly available RNA-seq dataset to perform differential expression analysis and pathway enrichment, then focuses specifically on the oxidative stress gene expression landscape in psoriatic lesional skin compared to healthy normal skin.

---

## Dataset

**GEO Accession:** [GSE54456](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE54456)

**Samples:** 92 lesional psoriatic skin biopsies and 82 normal skin biopsies (171 samples after QC filtering)

**Raw counts** were obtained via NCBI's standardized recount pipeline from the Sequence Read Archive (SRA).

**Original Publication:**

> Li B, Tsoi LC, Swindell WR, Gudjonsson JE, Tejasvi T, Johnston A, Ding J, Stuart PE, Xing X, Kochkodan JJ, Voorhees JJ, Kang HM, Nair RP, Abecasis GR, Elder JT. Transcriptome analysis of psoriasis in a large case-control sample: RNA-seq provides insights into disease mechanisms. *Journal of Investigative Dermatology*. 2014 Jul;134(7):1828-1838. doi: [10.1038/jid.2014.28](https://doi.org/10.1038/jid.2014.28)

---

## Analysis Overview

- **Quality Control** — PCA, sample distance heatmap, dispersion plot
- **Differential Expression** — DESeq2 with negative binomial model (padj < 0.05, |log2FC| > 1)
- **Visualization** — Volcano plot, heatmaps, dot plots
- **Oxidative Stress Gene Panel** — Targeted analysis of antioxidant, ROS-generating, and stress response genes
- **Pathway Enrichment (ORA)** — Gene Ontology enrichment (BP, MF, CC) using clusterProfiler
- **Pathway Enrichment (GSEA)** — Gene Set Enrichment Analysis on full ranked gene list

---

## Key Findings

### Completed
- Over 21,000 genes differentially expressed between psoriatic and normal skin (padj < 0.05)
- Applying |log2FC| > 1 threshold yields 1,290 upregulated and 1,419 downregulated genes
- PCA shows clean separation of conditions on PC1 (55% variance explained)
- From a curated oxidative stress gene panel, 5 genes were significantly differentially expressed:
  - **NOS2** — strongly upregulated (log2FC ~5.6), driving nitric oxide and ROS production
  - **SOD2** — upregulated, mitochondrial antioxidant response
  - **GPX2** — upregulated, glutathione peroxidase stress response
  - **HMOX1** — upregulated, classic oxidative stress response gene
  - **CYCS** — upregulated, consistent with mitochondrial apoptotic signaling
- GO enrichment confirms dominant inflammatory, skin barrier, and cell cycle signatures

### In Progress
- GSEA on full ranked gene list
- GSEA focused on oxidative stress pathways

## Repository Structure

```
├── main.R                        # Main analysis script
├── functions.R                   # Reusable functions
├── data/
│   ├── counts.tsv                # Raw count matrix (171 samples x 39,376 genes)
│   └── metadata.csv              # Sample metadata from SRA Run Selector
├── plots/
│   ├── PCA.png
│   ├── dispersion_plot.png
│   ├── sampledist_heatmap.png
│   ├── volcano_plot.png
│   ├── top50_heatmap.png
│   ├── oxidative_stress_heatmap.png
│   └── ora_barplots/
└── results/
    └── oxidative_stress_results.csv
```

---

## Requirements

### R Packages

```r
# CRAN
install.packages(c("tidyverse", "ggplot2", "pheatmap"))

# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
  "DESeq2",
  "org.Hs.eg.db",
  "AnnotationDbi",
  "EnhancedVolcano",
  "clusterProfiler"
))
```

---

## Usage

```r
# 1. Source functions
source("functions.R")

# 2. Run main analysis
source("main.R")
```

---

## Acknowledgements

Raw sequencing data were obtained from the NCBI Gene Expression Omnibus (GEO) database (accession GSE54456). All credit for data generation goes to Li et al. (2014) and the University of Michigan Department of Dermatology. This project is an independent re-analysis for educational purposes.

