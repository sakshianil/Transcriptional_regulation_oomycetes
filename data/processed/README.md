# Processed Data Files

## Overview
This directory contains processed datasets derived from the raw experimental data of the transcriptional regulation in oomycetes project. The data in this directory is organized and ready for downstream analyses, providing an important foundation for gene expression and regulatory network studies.

## Directory Contents

### 1. **`final_counts_path.txt`**
- **Description**: Contains the processed gene expression counts across different time points. The columns include:
  - **Gene_id**: A unique identifier for each gene.
  - **Description**: Functional annotation for each gene, if available.
  - **Time Series Columns** (e.g., `T5min_1`, `T5min_2`, `T5min_3`): Expression counts at various time points, with triplicate measurements indicated by suffixes `_1`, `_2`, `_3`.
- **Purpose**: Essential for time-series expression analysis, allowing you to explore how gene expression changes across developmental stages or stress conditions.
- **Format Example**:
  ```
  Gene_id        Description           T5min_1 T5min_2 T5min_3 T15min_1 ...
  CEG35065       hypothetical_protein  0       0       0       0       ...
  ```

### 2. **`final_counts_path_with_description.txt`**
- **Description**: An extended version of `final_counts_path.txt` that includes additional metadata for each gene:
  - **Geneid**: Gene identifier.
  - **Description**: Functional annotation of the gene.
  - **Uniprot_id**: UniProt accession ID, where available.
  - **Chr, Start, End, Strand**: Genomic location information.
  - **Length**: Gene length in base pairs.
  - **Time Series Columns** (e.g., `R5min_1`, `R5min_2`): Similar expression count columns as above but labeled differently.
- **Purpose**: Facilitates integrative analyses of gene expression and genomic features by providing genomic context alongside the processed expression counts.
- **Format Example**:
  ```
  Geneid     Description         Uniprot_id            Chr         Start   End     Strand Length  R5min_1 R5min_2 ...
  CEG41710   hypothetical_protein Uniprot/SPTREMBL:A0A0P1AL91 Scaffold_614 4315    4632    -     318     0       0     ...
  ```

### 3. **`README.md`**
- **Description**: A guide for understanding the contents of the `/processed` folder, explaining the purpose and format of each file.

## How These Files Were Generated
The processed data files in this directory were generated using a series of bioinformatics tools for trimming, mapping, and quantifying gene expression. This was done using a shell script that performs the following tasks:

- **Script Location**: [`trimming_mapping_host_pathogen.sh`](../../scripts/shell/Pl_halstedii_study/trimming_mapping_host_pathogen.sh)
- **Tools Used**:
  - **Trimmomatic**: Used for trimming adapters and removing low-quality bases from the raw reads.
  - **STAR Mapper**: Used to align the cleaned reads to the reference genome.
  - **FeatureCounts** (Subread package): Used to quantify the expression levels of genes based on the mapped reads.

These steps were scripted to ensure consistent preprocessing for all samples, making them suitable for further analyses.

## Usage
These processed files serve as the basis for various downstream analyses:

- **Differential Gene Expression Analysis**: To identify genes with significant changes in expression levels between different conditions or time points.
- **Clustering and Time Series Analysis**: Grouping genes with similar expression patterns to discover co-regulated genes or shared biological functions.
- **Functional Annotation**: Using `Uniprot_id` and gene descriptions to assign functions and explore pathways involved in transcriptional regulation.


## Script for Data Generation
The processed files in this folder were generated using the shell script located here:
[`trimming_mapping_host_pathogen.sh`](../../scripts/shell/Pl_halstedii_study/trimming_mapping_host_pathogen.sh). This script performs trimming, mapping, and feature counting of RNA-seq reads, allowing you to reproduce the data processing steps described above.



