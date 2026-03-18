# Big_data_analysis_in_genomics

## RNA-seq Analysis of ESCC vs Normal Tissue
Project description

This project presents an RNA-seq analysis comparing gene expression profiles between normal esophageal tissue and esophageal squamous cell carcinoma (ESCC) samples. The goal of the analysis was to identify differentially expressed genes and investigate biological processes and pathways associated with tumor development.

The analysis includes quality control, mapping evaluation, gene quantification, exploratory data analysis, differential expression analysis, and functional interpretation using enrichment methods.

### Data description

The dataset consists of RNA-seq samples obtained from publicly available sequencing data:
3 normal samples and
3 tumor samples (ESCC)

Sample IDs:
- Normal: SRR11647692, SRR11647693, SRR11647694
- Tumor: SRR11647702, SRR11647703, SRR11647704

The data includes:
- BAM alignment files
- Gene count matrix generated using featureCounts
- Sample metadata (sampleInfo.txt)
- Software and packages used

Software:
- R (version 4.3.3)
- Bash (Linux environment)
- RSeQC (for QC)
- featureCounts (for gene quantification)

### R packages

- **DESeq2** – differential expression analysis  
- **dplyr, tidyr** – data manipulation  
- **ggplot2, ggrepel** – data visualization  
- **pheatmap** – heatmaps  
- **RColorBrewer** – color palettes  
- **clusterProfiler** – enrichment analysis  
- **org.Hs.eg.db** – gene annotation  
- **reshape2** – data reshaping  
- **msigdbr** – MSigDB gene sets  
- **enrichplot** – enrichment visualization  
- **aPEAR** – pathway network visualization  
- **arules** – dependency for aPEAR  

### How to reproduce the analysis
1. **Preprocessing** (Bash)

Run preprocessing script:

    bash HW1_preprocess.sh

These steps include:

- Quality control (RSeQC)
- Alignment evaluation
- Generation of gene counts

2. **Run R analysis**

Run the full analysis script:

    Rscript HW1_R_dalis.R

or inside R:

    source("HW1_R_dalis.R")

3. **Outputs**

The analysis generates:

- QC plots (mapping rate, feature assignment)
- PCA plot
- Correlation heatmap
- Volcano plot
- MA plot
- Normalization comparison plots
- GO enrichment plots
- GSEA plots (GO and MSigDB)
- Pathway network visualization (aPEAR)


The analysis identifies significant transcriptional differences between normal and tumor samples. Functional analysis highlights enrichment of immune-related pathways, inflammatory signaling, and tumor-associated biological processes, supporting the role of gene expression dysregulation in ESCC.
