# Results

This document outlines the results generated from the scripts and workflows available in this GitHub repository for studying transcriptional regulation in oomycetes. The analyses covered include differential gene expression, regulatory motif analysis for transcription factors and binding sites studies and validation of effector genes across multiple genomes.

## 1. Differential Gene Expression (DGE) Analysis

### a. Differential Gene Expression Results
- **Tool Used**: DESeq2 R package.
- **Input**: The gene count matrix from **`data/processed/final_counts_path.txt`**.
- **Output**: A set of CSV files for each time-point comparison showing the significant genes with their Log2 Fold Changes.
  - **Example Output**:
    - `results_T15min_vs_T5min.csv`: Contains significant differentially expressed genes between 5 minutes and 15 minutes during zoospores release phase.
  
- **Visualization**:
  - **Before Normalization**: `counts_before_normalization.png` shows boxplots of gene counts for each sample before normalization.
  - **After Normalization**: `counts_after_normalization.png` shows boxplots after normalization.
  - **PCA Plot**: `PCA.png` visualizes the clustering of samples in reduced dimensions.

## 2. Clustering Analysis and Expression Profiles

### a. Gene Clustering
- **Tool Used**: DESeq2 and DEGreport.
- **Input**: Log2-transformed gene counts.
- **Output**: Clusters of genes with similar expression profiles.
  - **Files**:
    - `gene_clusters_minc15_reduceT.txt`: Contains the genes clustered together based on expression profiles.
    - `Phals_log2_transformed_counts.txt`: Log2-transformed count values for all genes.

### b. Cluster Visualizations
- **Gene Counts per Cluster**: `gene_counts_per_cluster.png` shows the number of genes assigned to each cluster.
- **Cluster Profiles**: Expression profiles for each cluster are visualized in `expression_profiles.png`.

## 3. Motif Analysis

### a. Motif Identification and Enrichment
- **Tools Used**: MEME Suite.
- **Input**: Promoter sequences from clustered genes.
- **Output**: Identified motifs enriched in promoters of each cluster.
  - **Files**:
    - `/results/meme/`: Contains the motif analysis results for each cluster.
    - `combined_fimo.tsv`: Summary of the motifs found using STREME and MEME.
  
### b. Conserved Motifs and Categories
- **Category Assignment**:
  - **Tool Used**: Custom shell scripts.
  - **Output**: `categorized_combined_fimo.tsv` - Contains motifs categorized based on conservation across species.

## 4. Effector Gene Analysis

### a. Effector Sequence Identification
- **Tool Used**: BLASTp.
- **Input**: Validated effector sequences from **`validated_effectors/effector_sequences.fa`**.
- **Output**: `effectors_blast_results.txt` lists the significant hits for effector proteins across the dataset. Furtheer processing generates dataset for blastn database.


## 5. Transcription Factor (TF) Analysis

### a. TF Binding Motif Analysis
- **Tools Used**: JASPAR, MEME Suite.
- **Input**: Promoters of TF-targeted genes.
- **Output**:
  - `combined_tomtom.tsv`: Contains the results from the motif similarity search against the JASPAR database.

### b. TF Conservation Analysis
- **FTFD Analysis**:
  - **Tool Used**: BLASTp against Fungal Transcription Factor Database (FTFD).
  - **Output**: Results stored in `blastp_output.txt` show homologous TFs across different species.

## 6. Summary of Key Findings

- **Differential Gene Expression**:
  - Significant changes in gene expression were observed at various time points of four asexual reproduction phases of pathogen. The PCA analysis showed distinct clusters of pathogencity-related gene groups enriched, indicating a clear temporal regulation of transcription.
  
- **Motif Discovery**:
  - Multiple clusters were enriched for known regulatory motifs, such as pathogencitiy related RxLR-like and CRN effector genes, which are indicative of host-pathogen interactions.

- **Effector and TF Analysis**:
  - A set of effectors conserved across multiple species were identified, along with their specific motifs, supporting their potential role in pathogenicity.
  - The TF analysis provided insight into the potential regulators of these conserved effectors.

## 7. Accessing the Results

All generated results are available within the `results/` directory, categorized by study:

- **Pl_halstedii_study**: Results for the Plasmopara halstedii gene expression and clustering analysis.
- **orthologs_motifs_study**: Results from ortholog analysis, conserved motif identification, and effector gene studies.
- **orthologs_TF_study**: Results for transcription factor motif discovery and similarity searches.

## 8. Interpretation Guide

- **Significance Levels**:
  - Differentially expressed genes were filtered using an adjusted p-value (padj) cutoff of 0.01.
  - Log2 Fold Change (LFC) was used to indicate the magnitude of differential expression.
  
- **Visualization**:
  - The various visual outputs provide an easy way to understand the results, such as the clustering of samples in PCA, heatmaps showing gene expression patterns, and plots of conserved motifs.

## 9. Additional Information

For further insights into the experimental design, refer to the [Methodology](../01_Methodology.md) section. For a list of all input data used in the analyses, refer to the [Data](../02_Data.md) section.


