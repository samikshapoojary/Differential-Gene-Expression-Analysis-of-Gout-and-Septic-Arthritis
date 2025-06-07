# Differential-Gene-Expression-Analysis-of-Gout-and-Septic-Arthritis
This script analyzes RNA-seq data from Gout, Septic Arthritis, and Healthy patients to identify key differentially expressed genes, visualize expression patterns, and assess clinical effects like sex and neutrophil count on gene expression. Developed as coursework for MSc Bioinformatics 2024/25, University of Glasgow.
---

## Overview

* Loaded and processed RNA-seq expression data, annotation files, and clinical sample information.
* Grouped samples by condition: Gout, SA, and Healthy Controls.
* Calculated sex ratios and compared neutrophil counts across groups.
* Identified significant genes differentially expressed between Gout vs HC and SA vs HC.
* Performed PCA for clustering analysis on significant genes.
* Investigated correlation of gene expression with sex and neutrophil count.
* Identified significantly different genes between Gout and SA.
* Created visualizations including histograms, violin plots, PCA plots, and volcano plots to summarize findings.

---

## How to Use

1. Install required R packages:

```r
install.packages("ggplot2")
install.packages("ggfortify")
```

2. Load data files and run the provided R script to reproduce analyses and plots.

3. Customize file paths in the script to your data locations.

---

## File Descriptions

* `Annotations.csv` — Gene annotation data.
* `DE_GOUT_vs_HC.csv` — Differential expression results for Gout vs Healthy Controls.
* `DE_SA_vs_HC.csv` — Differential expression results for Septic Arthritis vs Healthy Controls.
* `Expression_Table.csv` — Gene expression matrix.
* `Sample_Information.csv` — Clinical metadata for each sample.
