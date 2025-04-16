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
Decontamination and Filtering of Abundance Data 

**Code:** `DataPreparation.Rmd`

**Input:** 
- `bracken_merged_abbundances.num.txt` Raw abundances from Kraken2/Braken pipeline
- `metadata.xlsx` Metadata

**Results:**
 - `countsToT_decont.txt` - Abundances decontaminated
 - `countsToT_decont_normtss.txt` Abundances decontaminated and normalized (TSS)
 - `counts_decont_filtered.txt` Abundances decontaminated and filtered ( 0.1% relative abudance in at least 1 sample )
  
## Analysis 1: Methodological Comparison (RNA - DNA)
Evaluate the concordance and differences between DNA (metagenomics) and RNA sequencing approaches for the same set of samples.

**Code:** `ComparisonDNA-RNA.Rmd`  
- Merges datasets based on common samples.
- Assesses the overall similarity/dissimilarity between DNA and RNA samples based on their abundance profiles.
    - Distance matrix heatmap.
    - Dendrogram.
    - PCoA to plots the samples in a 2D space (based on Bray-Curtis). The plot helps visualize how samples cluster based on type (DNA/RNA) and individual ID.
- Shannon diversity for each sample
    - Boxplot and paired t-test between the two techniques.
- Taxonomic Comparison
    - Visualizes relative abundance differences using stacked bar plots.
- Information Gain : By discretizing abundances and calculating entropies, it determines the average Information Gain, indicating how much knowing the DNA abundance reduces uncertainty about the RNA abundance across the dataset.

**Input:** 
- `countsToT_decont_normtss.txt` Abundances decontaminated and normalized (TSS)
- `bracken_merged_abbundances_dna.num.txt` Dna Raw Abundances from Kraken2/Braken pipeline
- `metadata.xlsx` Metadata

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

