# Vulvar Microbiome Analysis Pipeline
This repository contains computational workflows for analyzing vulvar microbiome data through RNA sequencing and metagenomic approaches, as described in our publication. The analysis is structured into two primary components: microbiome characterization and methodological comparison.

## Repository Structure
```
project/
├── local/
│ ├── data/
├── results/
├── scripts/
│ ├── MainAnalysis.Rmd
│ └── ComparativeDNA-RNA.Rmd
└── README.md
```
<img src="https://github.com/user-attachments/assets/01ca6b41-6e4a-446d-84f4-dc33419d6ad0" width="600" height="500">

## Analysis 0: Decontamination and Filtering
**Code:** `DataPreparation.Rmd`
- Decontamination and Filtering of Abundance Data
**Input:** `local/data/bracken_merged_abbundances.num.txt`
**Results:**
  -`countsToT_decont.txt` - Abundances decontaminated
  -`countsToT_decont_normtss.txt` Abundances decontaminated and normalized (TSS)
  -`counts_decont_filtered.txt` Abundances decontaminated and filtered ( 0.1% relative abudance in at least 1 sample )
  
## Analysis 1: Methodological Comparison (RNA - DNA)
**Code:** `ComparisonDNA-RNA.Rmd`  
Evaluates RNA-seq against shotgun metagenomics for:
- Taxonomic concordance
- Sensitivity analysis for low-abundance taxa
- Correlation of diversity metrics between platforms

## Analysis 2: Vulvar Microbiome Characterization
**Code:** `MainAnalysis.Rmd`  
Processes RNA-seq data to:
- Import taxonomic abundances from Kraken/Bracken
- Perform **class discovery** through:
  - CST identification
  - Hierarchical Clustering 
- Calculate alpha diversity metrics:
  - Shannon Index
- Assess dysbiosis using
- Differential abundance testing

## Analysis 3: Cross Studies 

