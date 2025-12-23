# Vulvar Microbiome Analysis Pipeline

This repository contains computational workflows for analyzing vulvar microbiome data through RNA sequencing and metagenomic approaches, as described in our publication. The analysis is structured into two primary components: microbiome characterization and methodological comparison.   
![Fig 1_Microbiome (1)](https://github.com/user-attachments/assets/e809c287-ea82-436c-ad1c-6e0a8473ad8c)


---


- **[Analysis 0: Decontamination and Filtering](#analysis-0-decontamination-and-filtering)**, 
- **[Analysis 1: Methodological Comparison (RNA - DNA)](#analysis-1-methodological-comparison-rna---dna)** 
- **[Analysis 2: Vulvar Microbiome Characterization](#analysis-2-vulvar-microbiome-characterization)** 
- **[Analysis 3: Cross Studies](#analysis-3-cross-studies)** 
- **[Analysis 4: Host Differential Gene Expression](#analysis-4-host-differential-gene-expression)** 
- **[Analysis 5: Microbiome Data Correlation Analysis for Cytoscape](#analysis-5-Microbiome-Data-Correlation-Analysis-for-Cytoscape)**
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
Identifies differentially expressed genes (DEGs) between vaginal community state types (CSTs) and correlates them with epithelial keratinization pathways and dysbiosis scores.

**Code:** Differentially_Expressed_Genes.Rmd

**Output:** 
  - `../results/DEG_dysbio.xlsx`

**Input:**
  - `../local/data/GEP.count.xlsx` Raw RNA-seq counts

  - `../results/cluster_col.csv`  CST Classifications metadata

  - `../local/data/zscore_reactome_all.xlsx` Precomputed pathway scores

1.  **Data Loading and Preprocessing**:
    * Loads gene count data (`GEP.count.xlsx`) and sample cluster information (`cluster_col.csv`).
    * Maps microbial dominance categories to CST types (CSTI, CSTII, CSTIII, CSTIV).
    * Filters and normalizes gene counts using `edgeR` (TMM normalization).
    * Filters genes based on expression presence across clusters.

2.  **Differential Gene Expression (DGE) Analysis - CSTs**:
    * Performs DGE analysis using `edgeR`'s GLM QL F-test to compare different CSTs.
    * Identifies DEGs based on P-value and log-fold change thresholds.
    * Generates heatmaps (`pheatmap`) of DEGs for each contrast, scaled by row.
    * Saves lists of DEGs to an Excel file (`DEG_clusters.xlsx`).

3.  **Principal Component Analysis (PCA)**:
    * Selects the top 1000 most variable genes based on row variances.
    * Performs PCA on the normalized counts of these genes.
    * Generates PCA plots (PC1 vs. PC2) colored by:
        * Microbial `Dominance` (CST type).
        * `Keratinization` score (loaded from `zscore_reactome_all.xlsx`).
    * Combines plots using `patchwork`.

4.  **Keratinization Score Analysis**:
    * Visualizes the "Formation.of.the.Cornified.Envelope" pathway score (referred to as Keratinization Score) across different CSTs using boxplots.
    * Performs statistical comparisons between CSTs using `ggpubr::stat_compare_means`.

5.  **Pathway Enrichment Analysis - CST DEGs**:
    * For DEGs identified in each CST comparison:
        * Converts gene symbols to Entrez IDs using `org.Hs.eg.db`.
        * Performs Reactome pathway enrichment analysis using `ReactomePA::enrichPathway`.
        * Generates and saves dot plots (`plotCount` custom function) for significant pathways.

6.  **Differential Gene Expression (DGE) Analysis - Dysbiosis Score**:
    * Categorizes samples into "lower," "medium," and "higher" groups based on a `dysbio` score (derived from quantiles).
    * Performs DGE analysis between "lower" and "higher" dysbiosis groups.
    * Generates a heatmap of DEGs.
    * Saves the list of DEGs related to dysbiosis (`DEG_dysbio.xlsx`).

7.  **Pathway Enrichment Analysis - Dysbiosis DEGs**:
    * Performs Reactome pathway enrichment analysis for DEGs identified in the dysbiosis comparison, similar to the CST DEG enrichment.
    * Generates and saves dot plots.

---
   
## Analysis 5: Microbiome Data Correlation Analysis for Cytoscape

**Code:** `correlation_host_microbes.Rmd`

**Input:**
- `GEP.count.xlsx` Raw read counts for host genes
- `meta.xlsx` metadata file listing samples used in the analysis
- `bracken_merged_abundances.num.filtered.nocontam.txt` Filtered and decontaminated abundances from Kraken2/Bracken pipeline (output from Analysis 0)

**Results:**
- `results/edge_list.csv` - Edge list containing microbe-microbe, host-microbe and host-host correlations, p-values, and source/target for network visualization (e.g., in Cytoscape).
- `results/node_list.csv` - Node list that can be decorated with metadata column for Cytoscape analysis

---

This script documents a sophisticated **integrative multi-omics pipeline** designed to explore the interactions between the **host transcriptome** (RNA-seq) and the **vulvar microbiome** (taxonomic abundances).

By combining gene co-expression networks with cross-domain correlations, the analysis identifies how specific microbial taxa influence epithelial biological processes such as inflammation, keratinization, and immune response.

---

### **Analysis Workflow Overview**

The code is organized into several distinct functional blocks:

#### **1. Data Preprocessing & Normalization**

* **Host Data:** Utilizes `DESeq2` for filtering and size-factor normalization of gene expression counts.
* **Microbial Data:** Applies **CLR (Centered Log-Ratio)** transformation to Bracken-derived abundances to handle the compositional nature of microbiome data.
* **Variance Stabilization:** Uses `vst` to prepare host data specifically for network construction.

#### **2. WGCNA (Weighted Gene Co-expression Network Analysis)**

* Identifies modules of highly co-regulated genes.
* Determines the optimal soft-thresholding power (selected as 7 in this script) to achieve a scale-free topology (R^2 > 0.9).
* **Module-Trait Relationship:** Correlates these gene modules with microbial "dominance" clusters (e.g., *Lactobacillus*-dominated vs. others) to find biological signatures linked to specific community states.

#### **3. Cross-Domain Correlation & Network Construction**

* Calculates **Spearman correlations** between host genes and microbial taxa.
* Filters for high-confidence interactions (|R| > 0.6 and p < 0.05).
* Generates **Edge and Node lists** formatted for **Cytoscape**, allowing for the visualization of Host-Host, Microbe-Microbe, and Host-Microbe interaction networks.

#### **4. Functional Enrichment & Annotation**

* **EnrichR:** Performs pathway analysis (GO Biological Process, Reactome, etc.) for each WGCNA module.
* **Custom Annotation:** Maps host genes to specific epithelial functions like *apoptotic process*, *keratinization*, and *innate immune response* using `org.Hs.eg.db`.

#### **5. Microbial Functional Impact & "Bridge" Analysis**

* **Diverging Barplots:** Visualizes the percentage of positive vs. negative correlations between specific microbes (like *Gardnerella* or *Lactobacillus*) and host biological processes.
* **Bridge Analysis:** Identifies "Bridge Genes"—critical host nodes that interact with multiple distinct microbial taxa—visualized through a clustered heatmap to show patterns of co-regulation.

---

### **Key Biological Output**

The final outputs of this script are intended to answer:

1. Which host pathways are activated or suppressed by the presence of specific bacteria?
2. Which genes act as the primary "sensors" or "responders" to shifts in the microbiome?
3. How do different microbes differ in their molecular crosstalk with the host epithelium?


