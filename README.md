# scRNA-MSC-0steogenesis
# Identification of a surface marker of osteogenic tendency of mesenchymal stem cells

[![DOI](https://img.shields.io/badge/DOI-YOUR_PAPER_DOI-blue.svg)](https://doi.org/YOUR_PAPER_DOI)

This repository contains the code for the analysis in the manuscript titled: "**Identification of a surface marker of osteogenic tendency of mesenchymal stem cells**".

Our study identifies CD266 as a novel surface marker for isolating mesenchymal stem cells (MSCs) with enhanced osteogenic potential, leveraging single-cell RNA sequencing (scRNA-seq) analysis.

## Data Availability

The raw and processed sequencing data have been deposited in the Gene Expression Omnibus (GEO) database under accession number **GSEXXXXXX**.

## System Requirements

The analysis was performed on a Linux system and requires `conda` for environment management.

## Installation

To replicate the analysis, please follow these steps to create the computational environment:

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/your_username/MSC-Osteo-Marker.git
    cd MSC-Osteo-Marker
    ```

2.  **Create the conda environment:**
    All required packages are listed in the `environment.yml` file. Create and activate the environment using the following command:
    ```bash
    conda env create -f environment.yml
    conda activate hMSC
    ```

## Analysis Workflow

The analysis pipeline is organized into several scripts located in the `scripts/` directory. Please run them in the specified order.

*Before running, please modify the data, index, and result paths at the beginning of each script to match your local directory structure.*

1.  **`00_kallisto_quant.sh`**: Quantify transcript abundance from raw FASTQ files using Kallisto.

2.  **`01_tximport.R`**: Import Kallisto's quantification results into R and generate a gene-level count matrix.

3.  **`02_seurat_qc.R`**: Perform quality control (QC), filtering, and normalization on the count matrix using the Seurat package.

4.  **`03_pca_umap.R`**: Conduct Principal Component Analysis (PCA) and UMAP for dimensionality reduction.

5.  **`04_correlation_heatmap.R`**: Generate correlation heatmaps for genes of interest.

6.  **`05_go_analysis.R`**: Perform Gene Ontology (GO) enrichment analysis on identified gene clusters.

7.  **`06_expression_heatmap.R`**: Create heatmaps to visualize the expression of key genes across different cell conditions.

## Citation

If you use this code or data in your research, please cite our paper:

> Wengang Liu, Sihan Zhang, et al. (2024). "Identification of a surface marker of osteogenic tendency of mesenchymal stem cells". *Journal Name*.
