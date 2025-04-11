# Vulvar Microbiome Analysis Pipeline
This repository contains computational workflows for analyzing vulvar microbiome data through RNA sequencing and metagenomic approaches, as described in our publication. The analysis is structured into two primary components: microbiome characterization and methodological comparison.

## Repository Structure
```
project/
├── data/
├── results/
├── scripts/
│ ├── MainAnalysis.Rmd # Characterization & clustering
│ └── ComparativeDNA-RNA.Rmd
└── README.md
```
<img src="https://github.com/user-attachments/assets/01ca6b41-6e4a-446d-84f4-dc33419d6ad0" width="600" height="500">

## Analysis 1: Vulvar Microbiome Characterization
**Code:** `MainAnalysis.Rmd`  
Processes RNA-seq data to:
- Import taxonomic abundances from Kraken/Bracken
- Perform **class discovery** through:
  - Dirichlet Multinomial Mixtures (CST identification)
  - PAM clustering with silhouette width validation
- Calculate alpha diversity metrics:
  - Shannon Index
  - Chao1 Richness
- Assess dysbiosis using:
  - Community State Type (CST) prevalence analysis
  - Differential abundance testing (DESeq2)

## Analysis 2: Methodological Comparison
**Code:** `MethodComparison.Rmd`  
Evaluates RNA-seq against shotgun metagenomics for:
- Taxonomic concordance (Genus-level Bray-Curtis similarity)
- Sensitivity analysis for low-abundance taxa
- Computational resource benchmarking
- Correlation of diversity metrics between platforms

## Installation
