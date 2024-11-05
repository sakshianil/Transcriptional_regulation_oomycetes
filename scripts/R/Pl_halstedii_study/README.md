# README for Differential Gene Expression Analysis Script (1.2_DGE_analysis.R)

## Overview
This script performs differential gene expression (DGE) analysis for the transcriptional regulation project in oomycetes. It utilizes the R paackages such as `DESeq2` to normalize gene expression data, visualize results, and identify differentially expressed genes (DEGs) across multiple time points. The workflow includes normalization, PCA visualization, clustering, and exporting results.

## Requirements
The following R libraries are required to run the script:
- `BiocManager` (to install bioconductor packages)
- `DESeq2` (for differential expression analysis)
- `pheatmap` (for heatmap visualization)
- `ggplot2` (for creating plots)
- `clusterProfiler` (for enrichment analysis)
- `DEGreport` (for clustering and visualizing expression patterns)

If these libraries are not installed, the script will automatically install them.

## Input Data
- **Counts Data**: The script requires a counts file, formatted as a tab-delimited text file, with gene IDs as row names and samples as columns.
  - **Path**: Specify the counts data file path at the beginning of the script.
- **Metadata**: The script creates metadata based on the sample names and experimental conditions, including different time points and replicates.

## Usage
1. **Set Input and Output Paths**
   - Update the `counts_file` variable to point to your counts data file.
   - Update `output_directory` to the folder where you want to save results.

2. **Run the Script**
   - Open RStudio or run via terminal using the command:
     ```R
     Rscript DGE_analysis.R
     ```
   - The script will normalize counts, perform differential expression analysis, and save various plots and output files.

## Outputs
- **Gene Counts Boxplots**
  - Before and after normalization to visualize the effects of size factor adjustments.
- **PCA Plot**
  - Displays principal component analysis (PCA) of the transformed data to visualize sample clustering.
- **Heatmaps**
  - Heatmaps of sample correlations and gene clusters.
- **Differential Expression Results**
  - CSV files for each comparison, containing differentially expressed genes (DEGs).
  - A combined CSV file containing log2 fold change values for all conditions.
- **Gene Clustering Results**
  - CSV files containing genes grouped by similar expression profiles across conditions.

## Script Details
1. **Normalization and Size Factor Estimation**
   - Normalizes raw counts using DESeq2 and estimates size factors to control for sequencing depth.
2. **Differential Expression Analysis**
   - Uses a likelihood ratio test (LRT) to identify genes with significant changes across time points.
   - Results are further processed using `lfcShrink` for stable effect size estimates.
3. **Clustering Analysis**
   - Clusters genes based on their expression profiles and removes outliers.
   - Genes are grouped into clusters with similar patterns of expression over time.
4. **Visualization**
   - Visualizes gene expression data using PCA plots, heatmaps, and bar plots to show the number of genes in each cluster.

## How to Interpret the Results
- **Boxplots**: Boxplots for raw and normalized counts show the effect of normalization on reducing technical biases.
- **PCA**: The PCA plot helps in understanding the variation across samples, revealing potential batch effects or sample outliers.
- **Heatmaps**: Correlation heatmaps provide insights into the relationships between samples, while the clustering heatmap reveals gene groupings with similar expression dynamics.
- **Differential Expression Results**: CSV files contain information about significantly regulated genes, including log fold change and adjusted p-values for each pairwise comparison.

## Troubleshooting
- If the script stops at the installation of packages, ensure you have the correct permissions to install R packages.
- Make sure the input counts file is formatted correctly with gene IDs as row names and sample counts in columns.
- If color palette errors occur, you may need to adjust the number of colors used or use an alternative color palette.

## Contact
For questions or issues related to this script, please contact [Sakshi Bharti](mailto:sakshi.bharti@senckenberg.de).

## License
This script is distributed under the MIT License. Please see the LICENSE file in the root directory for more details.


