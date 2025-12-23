# Gene-Expression-Analysis-and-Interpretation
This repository contains the code and analysis for exploring HER2 amplification in breast cancer using the BRCA TCGA Pan-Cancer Atlas 2018 dataset. The project involves differential expression analysis, pathway enrichment, and survival modeling.
## Code Structure

The analysis is divided into sequential R scripts, numbered for execution order:

## 1. Data Preparation
*   **01_load_and_prepare_data.R**:
    *   Loads raw RNA-seq counts, clinical data, and Copy Number Alteration (CNA) data.
    *   Determines HER2 status (Amplified vs. Not Amplified) based on *ERBB2* CNVs.
    *   Merges clinical and expression data into processed files.

### 2. Differential Expression
*   **02_deseq2_DE_analysis.R**:
    *   Perform differential expression analysis using **DESeq2**.
    *   Identifies significantly differentially expressed genes (DEGs) between HER2-amplified and non-amplified tumors.
    *   Exports results and identifying the top 10 genes by fold change.

### 3. Pathway Analysis
*   **03_pathway_enrichment.R**:
    *   Uses **clusterProfiler** to perform Gene Ontology (GO) and KEGG pathway enrichment analysis on significant DEGs.
    *   Saves enrichment tables for biological interpretation.

### 4. Visualization
*   **04_vst_PCA_heatmap.R**:
    *   Applies Variance Stabilizing Transformation (VST) to count data.
    *   Generates a Principal Component Analysis (PCA) plot to visualize sample variance.
    *   Creates a heatmap of the top 50 DEGs to show expression patterns across samples.

### 5. Survival Analysis
*   **05_lasso_cox_glmnet.R**:
    *   Performs Lasso-regularized Cox Proportional Hazards regression using **glmnet**.
    *   Selects a subset of prognostic genes from the top 100 DEGs that best predict overall survival.
    *   Outputs the cross-validation curve and the list of selected prognostic genes with their coefficients.

## Results
All output files, including gene lists, plots (PDF/PNG), and enrichment tables, are saved in the 'results/' directory. .
