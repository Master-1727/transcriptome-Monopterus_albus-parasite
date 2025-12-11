# transcriptome-Monopterus_albus-parasite


Here’s a concise English “Methods / Pipeline” description you can drop into your GitHub project (e.g. in `README.md`). I’ll write it in a neutral, paper-friendly style and clearly mention the path / server setting and the two scripts: `work.r` and `drawmulltiwgcnanetwork.r`.

---

## Overview

This repository contains an end-to-end RNA-seq analysis pipeline for comparing parasite-infected and control samples across multiple tissues (intestine, kidney, liver, spleen). The workflow is implemented in R and covers:

* data import and QC
* PCA and sample clustering
* expression distribution and correlation analysis
* differential expression analysis
* GO and KEGG enrichment analysis
* multi-condition WGCNA (multiWGCNA)
* network visualization and module characterization (GSVA, enrichment)

The main analysis is performed in **`work.r`**, and a custom helper function for network visualization is defined in **`drawmulltiwgcnanetwork.r`**.

> ⚠️ **Important:** All file paths in the scripts (e.g. `setwd()`, input/output directories) are **hard-coded for our lab’s server** and must be adapted to your own directory structure before running.

---

## Scripts

### 1. `work.r` (main pipeline)

This script runs the full analysis pipeline in the following order:

1. **Environment setup and data import**

   * Clears the workspace, sets a random seed, and defines the working directory via `setwd()`.
   * Loads required R packages (e.g. `tidyverse`, `data.table`, `PCAtools`, `pheatmap`, `clusterProfiler`, `multiWGCNA`, `GSVA`).
   * Reads in:

     * gene-level TPM matrices and annotation tables
     * PCA input data
     * differential expression results for each tissue
     * precomputed sample correlation matrix
     * GO/KEGG enrichment results where applicable

2. **Principal Component Analysis (PCA)**

   * Reshapes the expression matrix to a gene × sample format.
   * Builds a sample metadata table (treatment vs control; tissue type) based on sample names.
   * Runs PCA using **PCAtools** and visualizes:

     * PCA colored by tissue (organ)
     * PCA colored by treatment condition
   * Highlights sample clustering patterns and potential outliers.

3. **Expression distribution and sample correlation**

   * Converts TPM values to log10 scale (with a small offset to avoid `log(0)`).
   * Groups samples by tissue and condition and plots:

     * **Boxplots** of log10(TPM) per group
     * **Heatmap** of sample–sample correlations using `pheatmap`

4. **Differential expression summary and visualization**

   * Standardizes DEG classification labels (`Up`, `Down`, `Nonsig`).
   * Summarizes DEG counts per tissue and direction (up/down/total) and visualizes them using a bar plot.
   * Generates:

     * **MA plots** highlighting the top up- and down-regulated genes
     * **Volcano plots** with significance thresholds (|logFC| > 1, FDR < 0.05) and gene labels for top hits

5. **GO enrichment analysis**

   * Converts gene identifiers (e.g. SYMBOL/GENENAME to ENTREZID) using a custom `OrgDb` object.
   * Performs GO enrichment (BP/CC/MF) with **clusterProfiler**.
   * For each tissue, imports GO enrichment results, standardizes ontology labels, and:

     * computes GO term similarity
     * draws GO tree plots (semantic similarity clusters)
     * generates simplified GO term sets and network-style visualizations of enriched GO terms

6. **KEGG enrichment and ridge plots**

   * Imports KEGG enrichment results for each tissue and filters by significance.
   * For each pathway, maps log2 fold changes of member genes onto **ridge density plots**, showing the distribution of logFC for genes within each KEGG pathway.
   * Produces one KEGG ridge plot per tissue and saves them as PDF files.

7. **Shared DEG analysis (Venn / UpSet / intersection GO)**

   * For each tissue, extracts significant DEGs (excluding non-significant genes) and removes duplicate gene IDs.
   * Constructs:

     * a **Venn diagram** of DEGs across tissues
     * an **UpSet plot** showing intersection sizes between tissue-specific DEG sets
   * Identifies the intersection of DEGs shared by all four tissues, maps them to annotation, and performs GO enrichment on the shared gene set, followed by network-style visualization of enriched GO terms.

8. **multiWGCNA network construction and module analysis**

   * Calls the custom function `drawMultiWGCNAnetwork2()` (from `drawmulltiwgcnanetwork.r`) to generate publication-ready network plots highlighting the module of interest and its overlaps across conditions.

   * Preprocesses the expression matrix:

     * removes genes with zero variance
     * keeps the top 75% most variable genes
     * applies log1p transformation
   * Builds a sample metadata table including **Status** (Treatment vs Controls) and **Tissue**.
   * Uses **multiWGCNA** to:

     * construct networks across conditions (combined, tissue, status)
     * detect modules and calculate module eigengenes
     * assess module–trait relationships
   * Computes:

     * **overlap** between modules across conditions
     * **differential module expression (DME)** between treatment and control
     * **preservation statistics** across disease status and tissues
   * Extracts genes belonging to a specific module of interest (e.g. `Treatment_015`), performs GO and KEGG enrichment for that module, and saves the results.

9. **Network visualization and GSVA**

   * Saves the gene IDs of the selected module as an `.rds` object.
   * Runs **GSVA** on the expression matrix using the module gene set as a “pathway” and plots:

     * GSVA enrichment scores across tissues and disease status
     * tissue-stratified boxplots for the module GSVA scores

---

### 2. `drawmulltiwgcnanetwork.r` (custom multiWGCNA network plot)

This script defines a modified version of the multiWGCNA plotting function:

* **`drawMultiWGCNAnetwork2()`**

Key features:

* Accepts:

  * the list of WGCNA networks (`WGCNAlist`)
  * the list of overlap comparison results (`comparisonList`)
  * the module ID of interest (e.g. `"Treatment_015"`)
  * the experimental design table (`design`)
* Filters out outlier modules if requested.
* Builds an **igraph** object where:

  * nodes represent modules from different networks/conditions
  * edges represent significant overlaps between modules (weighted by –log10 adjusted P-value)
* Colors nodes by condition (e.g. combined, Treatment, Controls, tissue-specific networks) and scales edge width by overlap significance.
* Produces a compact network layout suitable for inclusion in figures or supplementary materials.

This script is sourced by `work.r` and does **not** need to be run independently, as long as it is placed in the same project and sourced before calling `drawMultiWGCNAnetwork2()`.

---

## Dependencies

The pipeline uses (non-exhaustive list):

* **Core tidyverse / data handling**
  `tidyverse`, `data.table`, `stringr`, `readr`, `dplyr`, `tibble`
* **Visualization**
  `ggplot2`, `ggforce`, `ggrepel`, `pheatmap`, `ggridges`, `ComplexUpset`, `ggvenn`, `colorspace`
* **Dimension reduction / QC**
  `PCAtools`
* **Enrichment analysis**
  `clusterProfiler`, `DOSE`, `GOSemSim`, `aPEAR`, custom `OrgDb` objects, `enrichplot` (optional)
* **Network analysis**
  `multiWGCNA`, `WGCNA`, `igraph`
* **Pathway and gene set analysis**
  `GSVA`
* **Fonts and export**
  `showtext`, `extrafont`, `RColorBrewer`

Please ensure all required packages are installed before running the pipeline.

---

