# Vulvar Microbiome Analysis Pipeline

This repository contains computational workflows for analyzing vulvar microbiome data through RNA sequencing and metagenomic approaches, as described in our publication. The analysis is structured into two primary components: microbiome characterization and methodological comparison.   
![Fig 1_Microbiome (1)](https://github.com/user-attachments/assets/e809c287-ea82-436c-ad1c-6e0a8473ad8c)


---


- **[Analysis 0: Decontamination and Filtering](#analysis-0-decontamination-and-filtering)**, 
- **[Analysis 1: Methodological Comparison (RNA - DNA)](#analysis-1-methodological-comparison-rna---dna)** 
- **[Analysis 2: Vulvar Microbiome Characterization](#analysis-2-vulvar-microbiome-characterization)** 
- **[Analysis 3: Cross Studies](#analysis-3-cross-studies)** 
- **[Analysis 4: Host Differential Gene Expression](#analysis-4-host-differential-gene-expression)** 

---


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
   
--- 

## Analysis 1: Methodological Comparison (RNA - DNA)
Evaluate the concordance and differences between DNA (metagenomics) and RNA sequencing approaches for the same set of samples.

**Code:** `ComparisonDNA-RNA.Rmd`  
**Input:** 
- `countsToT_decont_normtss.txt` Abundances decontaminated and normalized (TSS)
- `bracken_merged_abbundances_dna.num.txt` Dna Raw Abundances from Kraken2/Braken pipeline
- `metadata.xlsx` Metadata

**Analysis:**
- Merges datasets based on common samples.
- Assesses the overall similarity/dissimilarity between DNA and RNA samples based on their abundance profiles.
    - Distance matrix heatmap.
    - Dendrogram.
    - PCoA to plots the samples in a 2D space (based on Bray-Curtis). The plot helps visualize how samples cluster based on type (DNA/RNA) and individual ID.
      
- Shannon diversity for each sample
    - Boxplot and paired t-test between the two techniques.

<img src="https://github.com/user-attachments/assets/927d9632-9877-4604-846f-b2faa053f9cf" width="300" height="200">



- Taxonomic Comparison
    - Visualizes relative abundance differences using stacked bar plots.
- Information Gain : By discretizing abundances and calculating entropies, it determines the average Information Gain, indicating how much knowing the DNA abundance reduces uncertainty about the RNA abundance across the dataset.

--- 

## Analysis 2: Vulvar Microbiome Characterization
Processes microbiome data to characterize community structure, identify community types (CSTs), assess diversity and dysbiosis, and perform differential abundance analysis between identified community types.

**Code:** `VulvarCharacterization.Rmd`  
**Input:**
-  `../results/counts_decont_filtered.txt`   Abundances decontaminated and filtered ( 0.1% relative abudance in at least 1 sample )
-  `../local/data/MicroPhenoDBAssociationScores.csv` External database file containing scores linking microbes to Bacterial Vaginosis.
  
**Analysis:**
- Normalization: Total Sum Scaling (TSS) normalization on input count data.
  
- Visualization of community structure in samples:
    - Principal Component Analysis (PCA) on scaled TSS data and plots PC1 vs PC2.
    - Non-Metric Multidimensional Scaling (NMDS) using Bray-Curtis distance on TSS data and plots NMDS1 vs NMDS2.
      
- Clustering & Community Typing (CST Identification):
    - Bray-Curtis distance matrix.
    - Hierarchical Clustering (Euclidean distance, ward.D2).
    - Dominant species identification in each sample as a proxy for Community State Types (CSTs).
    - Stacked bar plot showing dominant species/CST composition across samples.  ( <img src="https://github.com/user-attachments/assets/013ffe81-eb08-4ad4-a548-d711fe882b95" width="90" height="30"> )

 - **Heatmap of row-scaled TSS data**, annotated by cluster and CSTs.
   
       - Differential Abundance Testing (Between Clusters):
          Normalizes counts using edgeR TMM method and calculates logCPM.
          Performs Wilcoxon rank-sum tests comparing taxon abundance (logCPM) between each cluster and all others.
          Identifies Differentially Abundant Microorganisms (DAMs) based on FDR and median change thresholds.
          Generates a heatmap visualizing only the significant DAMs.
      
       - Cluster Composition Visualization:
          Calculates the average taxonomic composition for each dominant species group (CST).
          Generates pie charts visualizing the average composition for each identified CST.

- Alpha Diversity Analysis:
    - Rarefies count data to the minimum library size.
    - Shannon diversity index calculation.
       - Box plot comparing Shannon diversity between CSTs, including statistical comparisons.
       - Dot plot visualizing Shannon diversity per sample.
 
<img src="https://github.com/user-attachments/assets/e4afc8f8-c618-474a-b771-92726ea1bc0d" width="900" height="70">



- Dysbiosis Score Assessment:
    - Loading of external microbe-phenotype association scores.
    - Sample-specific dysbiosis score calculation based on weighted microbial abundances.
       - Dot plot visualizing dysbiosis scores per sample.

<img src="https://github.com/user-attachments/assets/117bffec-cc94-4566-aa79-f71f82e4b6ad" width="900" height="60">

- Genus-Level Analysis
   - Aggregates TSS data to the genus level.
   - Filters for 5% most abundant genera.
   - Generates genus-level stacked bar plots and bubble plots.
  
<img src="https://github.com/user-attachments/assets/33a1ce86-3b1b-44ee-b582-c9fd6a7ad899" width="500" height="300">

  Creates faceted box plots comparing genus abundances across CSTs.

- Correlation Analysis:
    - Calculates and visualizes the correlation between Shannon diversity and the dysbiosis score using a scatter plot.


---

## Analysis 3: Cross Studies 
Integrates microbiome data from different body sites (vulvar, vaginal, anal) and different studies/batches. 
Corrects for batch effects using ComBat and compares community structure and composition.

Code:`ComparativeOtherStudy.Rmd` 

- Data Loading & Preprocessing
- Normalization : Total Sum Scaling (TSS) normalization.
- Batch Effect Correction using ComBat function (sva package) to the TSS data to adjust for batch effects defined in the metadata.        
- Performs NMDS on the batch-corrected abundance data (combdata).
        
<img src="https://github.com/user-attachments/assets/75c8cfe5-1c4d-4897-b3ac-8636dd78fea6" width="500" height="350">

- Compositional Comparison (Top Taxa):
   - Euler diagram visualizing the overlap of the Top 100 abundant taxa sets between the vaginal, anal, and vulvar sites.
     
   <img src="https://github.com/user-attachments/assets/4884853b-c514-4019-a1e2-060002993905" width="300" height="230">

- Alpha Diversity Analysis & Dysbiosis Score Assessment

**Input:**

  - `../local/data/df_metadataPRJEB61325.txt` Metadata file for samples, must contain site and batch information.
  - `../local/data/bracken_merged_abbundances.num_PRJEB61325.txt` Raw taxonomic abundance counts (Kraken/Bracken output) for PRJEB61325 samples.
  - `../results/countsToT_decont.txt` Pre-processed (e.g., decontaminated) count data for a second set of samples (potentially from Analysis 2).
  - `../local/data/Vulvodinia_Samples.xlsx` Spreadsheet with sample information (loaded but not explicitly used in the main analysis steps shown).
  - `../local/data/MicroPhenoDBAssociationScores.csv` External database file linking microbes to phenotype/disease association scores.

---


## Analysis 4: Host Differential Gene Expression


