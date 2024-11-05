
# Script Name: 3_DGE_analysis.R
# Author: Sakshi Bharti
# Description: This script performs differential gene expression analysis using DESeq2, generates visualizations, and clusters expression profiles for RNA-Seq data of oomycetes.
# Dependencies: R (>= 4.0), DESeq2, ggplot2, pheatmap, clusterProfiler, DEGreport,RColorBrewer, tidyr, dplyr, tibble, ggforce

# ---------------------------
# Load Necessary Libraries
# ---------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("DESeq2")) BiocManager::install("DESeq2")
if (!require("DEGreport")) BiocManager::install("DEGreport")
if (!require("pheatmap")) install.packages("pheatmap")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("clusterProfiler")) install.packages("clusterProfiler")
if (!require("RColorBrewer")) install.packages("RColorBrewer")
if (!require("tidyr")) install.packages("tidyr")
if (!require("dplyr")) install.packages("dplyr")
if (!require("tibble")) install.packages("tibble")
if (!require("ggforce")) install.packages("ggforce")
# Load the libraries
library(DESeq2)
library(pheatmap)
library(ggplot2)
library(clusterProfiler)
library(DEGreport)
library(RColorBrewer)
library(tidyr)
library(dplyr)
library(tibble)
library(ggforce)

# ---------------------------
# Define Input and Output Paths
# ---------------------------

counts_file <- "/path_to/Transcriptional_regulation_oomycetes/data/processed/final_counts_path.txt"
output_directory <- "/path_to/Transcriptional_regulation_oomycetes/results/R/Pl_halstedii_study"
dir.create(output_directory, showWarnings = FALSE)

# ---------------------------
# Load and Prepare Data
# ---------------------------
# Load counts data
counts_data_raw <- read.table(counts_file, header = TRUE, row.names = 1, sep = "\t")

# Remove "Gene description" column for DESeq2 analysis
counts_data <- counts_data_raw[, -1]

# Rename columns to sample names
colnames(counts_data) <- c("T5min_1", "T5min_2", "T5min_3", "T15min_1", "T15min_2", "T15min_3", "T4h_1", "T4h_2", "T4h_3", "T8h_1", "T8h_2", "T8h_3", "T12h_1", "T12h_2", "T12h_3", "T24h_1", "T24h_2", "T24h_3", "T48h_1", "T48h_2", "T48h_3", "T72h_1", "T72h_2", "T72h_3", "T120h_1", "T120h_2", "T120h_3", "T221h_1", "T221h_2", "T221h_3", "T288h_1", "T288h_2", "T288h_3", "T290h_1", "T290h_2", "T290h_3", "T292h_1", "T292h_2", "T292h_3", "T294h_1", "T294h_2", "T294h_3", "T296h_1", "T296h_2", "T296h_3", "T296h_prim_1", "T296h_prim_2", "T296h_prim_3")

# ---------------------------
# Create Metadata
# ---------------------------
metadata <- data.frame(
  row.names = colnames(counts_data),
  condition = factor(rep(c("T5min", "T15min", "T4h", "T8h", "T12h", "T24h", "T48h", "T72h", "T120h", "T221h", "T288h", "T290h", "T292h", "T294h", "T296h", "T296h_prim"), each = 3)),
  samples = colnames(counts_data),
  time = factor(rep(1:16, each = 3))
)

# ---------------------------
# Create DESeqDataSet
# ---------------------------
dds <- DESeqDataSetFromMatrix(countData = counts_data, colData = metadata, design = ~ condition)

# Filter out lowly expressed genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# ---------------------------
# Create Color palette
# ---------------------------

palette_colors <- brewer.pal(n = 12, name = "Set3") 
additional_colors <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3") # Additional colors
all_colors <- c(palette_colors, additional_colors)

# Ensure unique conditions from metadata
unique_conditions <- unique(metadata$condition)
num_conditions <- length(unique_conditions)

# Check if additional colors are needed
if (num_conditions > length(all_colors)) {
    stop("Not enough colors in the palette to match the number of conditions.")
}

# Define specific color palette for conditions
condition_colors <- setNames(all_colors[1:num_conditions], unique_conditions)


# ---------------------------
# Plot Counts Before Normalization
# ---------------------------

counts_melt <- counts_data %>%
  rownames_to_column(var = "Gene") %>%
  pivot_longer(cols = -Gene, names_to = "Sample", values_to = "Count")

counts_melt$Sample <- factor(counts_melt$Sample, levels = colnames(counts_data))
counts_melt$condition <- factor(gsub("_[^_]+$", "", counts_melt$Sample), levels = unique(metadata$condition))

# Generate boxplot before normalization
p1 <- ggplot(counts_melt, aes(x = Sample, y = Count, fill = condition)) +
  geom_boxplot() +
  scale_y_log10() +
  labs(title = "Gene counts before normalization", x = "Samples", y = "Counts (log scale)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = rep(condition_colors, length.out = nlevels(counts_melt$condition)))

# Save the plot
plot_filepath <- file.path(output_directory, "counts_before_normalization.png")
ggsave(filename = plot_filepath, plot = p1, width = 10, height = 8)

# ---------------------------
# Normalize Counts
# ---------------------------


dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized = TRUE)

# Plot Counts After Normalization
normalized_counts_melt <- normalized_counts %>%
  as.data.frame() %>%
  rownames_to_column(var = "Gene") %>%
  pivot_longer(cols = -Gene, names_to = "Sample", values_to = "Normalized_Count")

normalized_counts_melt$Sample <- factor(normalized_counts_melt$Sample, levels = colnames(counts_data))
normalized_counts_melt$condition <- factor(gsub("_[^_]+$", "", normalized_counts_melt$Sample), levels = unique(metadata$condition))

p2 <- ggplot(normalized_counts_melt, aes(x = Sample, y = Normalized_Count, fill = condition)) +
  geom_boxplot() +
  scale_y_log10() +
  labs(title = "Gene counts after normalization", x = "Samples", y = "Counts (log scale)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = rep(condition_colors, length.out = nlevels(normalized_counts_melt$condition)))

# Save the plot
plot_filepath <- file.path(output_directory, "counts_after_normalization.png")
ggsave(filename = plot_filepath, plot = p2, width = 10, height = 8)


# ---------------------------
# Pricipal Component Analysis
# ---------------------------


#variance stablizing transformation on the samples
vst_dds <- vst(dds, blind = FALSE)

# PCA plot
pca_dat <- plotPCA(vst_dds, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_dat, "percentVar"))


# Assign colors based on conditions
condition_colors_pca <- rep(condition_colors, each = 3)[1:nrow(metadata)]

# Create PCA plot with ggplot2
p3 <- ggplot(pca_dat, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  geom_mark_ellipse(aes(label = condition, color = condition), show.legend = FALSE) +
  scale_color_manual(values = condition_colors_pca) +
  theme_minimal() +
  labs(title = "PCA of Pathogen Samples", 
       x = paste0("PC1: ", percentVar[1], "% variance"), 
       y = paste0("PC2: ", percentVar[2], "% variance")) +
  theme(legend.position = "none")

#Save the plot
plot_filepath <- file.path(output_directory, "PCA.png")
ggsave(filename = plot_filepath, plot = p3, width = 10, height = 8)



# ---------------------------
# Differential Expression Analysis
# ---------------------------
# Run DESeq2 analysis

# Relevel the condition factor to ensure the reference level is appropriate
dds$condition <- relevel(dds$condition, ref = "T5min")

# Run the DESeq2 differential expression pipeline

# For a full time-series analysis, consider using a likelihood ratio test (LRT)
dds_lrt <- DESeq(dds, test = "LRT", reduced = ~ 1)

# Extract the results for the comparison of interest
# Define the specific condition pairs
condition_pairs <- data.frame(
 
  condition1 = c("T15min", "T4h", "T8h", "T12h", "T24h", "T48h", "T72h", "T120h", "T221h", "T288h", "T290h", "T292h", "T294h", "T296h", "T296h"),
  condition2 = c("T5min", "T15min", "T4h", "T8h", "T12h", "T24h", "T48h", "T72h", "T120h", "T221h", "T288h", "T290h", "T292h", "T294h", "T296h_prim")
)


# Initialize a list to store results
results_list <- list()

# Initialize an empty data frame to store combined results
combined_results <- data.frame(gene = character())  # Start with an empty data frame with a gene column

# Loop through each condition pair
for (i in 1:nrow(condition_pairs)) {
  condition1 <- condition_pairs$condition1[i]
  condition2 <- condition_pairs$condition2[i]
  
  # Define results name
  results_name <- paste0(condition1, "_vs_", condition2)
  print(paste("Processing:", results_name))
  
  # Perform differential expression analysis
  results_lfc <- results(dds_lrt, contrast = c("condition", condition1, condition2), alpha = 0.01, lfcThreshold = 0)
  results_lfc_shrink <- lfcShrink(dds_lrt, contrast = c("condition", condition1, condition2), alpha = 0.01, lfcThreshold = 0, type = "normal", rownames=1)
  
  # Process results
  results_tb <- results_lfc_shrink %>%
    as.data.frame() %>%
    rownames_to_column(var = "gene") %>%
    as_tibble()
  
  # Filter significant results
  padj_cutoff <- 0.01
  results_fc <- as.data.frame(results_tb %>%
    filter(padj < padj_cutoff))
  
  # Save results to CSV file
  csv_filename <- paste0("results_", condition2, "_", condition1, ".csv")
  csv_filepath <- file.path(output_directory, csv_filename)  
  write.csv(results_fc, file = csv_filepath, row.names = FALSE)
  
  # Store results in list (optional)
  results_list[[results_name]] <- list(
    results = results_lfc,
    results_shrink = results_lfc_shrink,
    results_tb = results_tb,
    results_fc = results_fc
  )
# Select relevant columns (gene and log2FoldChange) for the combined results
  combined_results_subset <- results_tb %>%
    select(gene, log2FoldChange) %>%
    rename(!!paste0("foldchange_", condition2, "_", condition1) := log2FoldChange)

  # Merge the results with the combined_results data frame
  if (nrow(combined_results) == 0) {
    combined_results <- combined_results_subset
  } else {
    combined_results <- full_join(combined_results, combined_results_subset, by = "gene")
  }
}

# Save the combined results to a CSV file
combined_output_filepath <- file.path(output_directory, "combined_DE_results.csv")
write.csv(combined_results, file = combined_output_filepath, row.names = FALSE)


# Example: Checking if files are saved
list.files(output_directory,pattern = "results_.*\\.csv")

# Assuming you have a results_fc data frame
# Extract rownames
sig_genes <- results_fc %>% pull(gene)


# ---------------------------
# Save and Visualize Clusters
# ---------------------------
# Clustering significant genes

clustering_sig_genes <-  results_fc %>% arrange(padj) %>% head(n=11705)

#For clustering the genes (vst, blind=T) are the parameters
vst_dds_T <- vst(dds, blind = T)
vst_mat_T <- assay(vst_dds_T)
cluster_rlog <- vst_mat_T[clustering_sig_genes$gene, ]

# Assuming cluster_rlog and metaData are defined
# Save log2 transformed counts

write.table(x = as.data.frame(cluster_rlog), file = file.path(output_directory, "Phals_log2_transformed_counts.txt"), sep = '\t', quote = FALSE)


# Clustering the profiles by cutting the tree (hierarchical clustering) and making groups/clusters with same expression profiles
# Reduce=T removes the outliers, minc=15 at least 15 genes in a group
clusters <- degPatterns(cluster_rlog, metadata = metadata, time = "time", col = NULL, plot = FALSE, reduce = TRUE, minc = 15)

# Save genes in clusters
cluster_groups <- clusters$df
write.table(x = as.data.frame(cluster_groups), file = file.path(output_directory, "gene_clusters_minc15_reduceT.txt"), sep = '\t', quote = FALSE)

# Save the list of survived clusters
survived_clusters <- clusters$pass
write.table(x = as.data.frame(survived_clusters), file = file.path(output_directory, "number_clusters_minc15_reduceT.txt"), sep = '\t', quote = FALSE)


# Visualize the number of genes per cluster
gene_counts_per_cluster <- clusters$df %>%
  group_by(cluster) %>%
  summarize(gene_count = n())

# Create bar plot with number of genes annotated on each bar
p4 <- ggplot(gene_counts_per_cluster, aes(x = factor(cluster), y = gene_count)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    geom_text(aes(label = gene_count), vjust = -0.3, color = "black", size = 2.5) +
    theme_minimal() +
    labs(title = "Number of Genes in Each Cluster",
         x = "Cluster",
         y = "Number of Genes") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p4)


# Save the plot
plot_filepath <- file.path(output_directory, "gene_counts_per_cluster.png")
ggsave(filename = plot_filepath, plot = p4, width = 10, height = 8)



# ---------------------------
# Heatmap of fold changes
# ---------------------------


# Rename the column in clusters$df from 'genes' to 'gene'
clusters$df <- clusters$df %>%
  rename(gene = genes)

# Combine combined_results with clusters$df based on the gene column
final_combined_results <- clusters$df %>%
    left_join(combined_results, by = c("gene" = "gene"))

# Save the final combined results to a CSV file
final_output_filepath <- file.path(output_directory, "clustered_and_DEG_genes.csv")
write.csv(final_combined_results, file = final_output_filepath, row.names = FALSE)


# Prepare the data for the heatmap
heatmap_data <- final_combined_results %>%
  select(gene, starts_with("foldchange"), cluster) %>%
  column_to_rownames("gene") %>%
  as.matrix()

# Define colors for the heatmap
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)

# Plot the heatmap
p5 <- pheatmap(heatmap_data[, -ncol(heatmap_data)],  # Exclude cluster column for heatmap
         color = color_palette, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE,
         show_rownames = FALSE, 
         show_colnames = TRUE, 
         main = "Heatmap of Fold Changes",
         annotation_row = data.frame(cluster = heatmap_data[, "cluster"]))

# Save the plot
plot_filepath <- file.path(output_directory, "heatmap.png")
ggsave(filename = plot_filepath, plot = p5, width = 10, height = 8)



# ---------------------------
# Additional Visualizations
# ---------------------------

# Create 'class' column in pathogen_counts
pathogen_counts <- counts_data_raw  %>% 
   rownames_to_column(var = "gene") %>% mutate(class = case_when(
        grepl("hypothetical\\s*[_]?\\s*protein", counts_data_raw$Description, ignore.case = TRUE) ~ "hypothetical protein",
        grepl("RxLR\\s*[-]?\\s*like", counts_data_raw$Description, ignore.case = TRUE) ~ "RxLRs",
        grepl("CRN\\s*[-]?\\s*like", counts_data_raw$Description, ignore.case = TRUE) ~ "CRNs",
        grepl("Phospholipase\\s*|\\s*Phosphatidyl\\s*inositol\\s*|\\s*Phosphatidyl\\s*inositol", counts_data_raw$Description, ignore.case = TRUE) ~ "Phospholipid signalling",
        grepl("NLP\\s*[-]?\\s*like", counts_data_raw$Description, ignore.case = TRUE) ~ "NRPs/NPPs",
        grepl("transport|channel|calcineurin|Calmodulin", counts_data_raw$Description, ignore.case = TRUE) ~ "Transporters/Channels",
        grepl("ase", counts_data_raw$Description, ignore.case = TRUE) ~ "Proteases",
        grepl("inhibitor", counts_data_raw$Description, ignore.case = TRUE) ~ "Inhibitor proteins",   
        grepl("ribosom|RNA",counts_data_raw$Description, ignore.case = TRUE) ~ "Ribosomal and non-coding",
        grepl("sec|rab|gtpase|exocyst|snare", counts_data_raw$Description, ignore.case = TRUE) ~ "Vesicle transport proteins",
        grepl("transcription\\s*[-]?\\s*factor|ccaat-binding|transcription\\s*[-]?\\s*initiation|myb\\s*[-]?\\s*like|heat\\s*[-]?\\s*shock|zinc finger|helix\\s*[-]?\\s*turn|bzip|mads box|helip\\s*[-]?\\s*loop|", counts_data_raw$Description, ignore.case = TRUE) ~ "TF proteins", 
        TRUE ~ "others"))
# If clusters$df already has a column named 'gene', you might not need to rename it
colnames(clusters$df)[which(colnames(clusters$df) == "V1")] <- "gene" 

# Merge the data
merged_data <- clusters$df %>%
  left_join(pathogen_counts, by = "gene")

# Summarize the gene counts per cluster and class
gene_counts_per_cluster_class <- merged_data %>%
  group_by(cluster, class) %>%
  summarize(gene_count = n(), .groups = 'drop')

# Create stacked bar plot
p6 <- ggplot(gene_counts_per_cluster_class, aes(x = factor(cluster), y = gene_count, fill = class)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = gene_count), position = position_stack(vjust = 0.5), color = "black", size = 2.5) +
    theme_minimal() +
    labs(title = "Number of Genes in Each Cluster",
         x = "Cluster",
         y = "Number of Genes",
         fill = "Gene Class") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the plot
plot_filepath <- file.path(output_directory, "Gene_groups.png")
ggsave(filename = plot_filepath, plot = p6, width = 10, height = 8)


# Create a phase mapping for foldchange columns
phase_mapping <- c(
  "foldchange_T5min_T15min" = "ZS","foldchange_T15min_T4h" = "IF", "foldchange_T4h_T8h" = "IF","foldchange_T8h_T12h" = "IF", "foldchange_T12h_T24h" = "IF",
  "foldchange_T24h_T48h" = "IF", "foldchange_T48h_T72h" = "IF","foldchange_T72h_T120h" = "CO", "foldchange_T120h_T221h" = "CO","foldchange_T221h_T288h" = "CO",
  "foldchange_T288h_T290h" = "SP", "foldchange_T290h_T292h" = "SP","foldchange_T292h_T294h" = "SP", "foldchange_T294h_T296h" = "SP", "foldchange_T296h_prim_T296h" = "SP"
)

# Reorganize the foldchange columns by phases
final_combined_results_long <- final_combined_results %>%
  pivot_longer(cols = starts_with("foldchange"), names_to = "condition", values_to = "foldchange") %>%
  mutate(phase = phase_mapping[condition]) %>%
  arrange(cluster, phase, gene)

# Prepare data for plotting
line_plot_data <- final_combined_results_long %>%
  group_by(cluster, phase) %>%
  summarise(mean_foldchange = mean(foldchange), .groups = 'drop')

# Plot the line plot
p7<- ggplot(line_plot_data, aes(x = phase, y = mean_foldchange, color = cluster, group = cluster)) +
  geom_line() +
  geom_point() +
  labs(title = "Fold Changes by Phase and Cluster", x = "Phase", y = "Mean Fold Change") +
  theme_minimal()

# Save the plot
plot_filepath <- file.path(output_directory, "Phasee_mapping.png")
ggsave(filename = plot_filepath, plot = p7, width = 10, height = 8)




#create expression profiles of all clusters
final_combined_results_fc <- clusters$df %>%
    left_join(results_fc, by = c("gene" = "gene"))


candidate_genes <- final_combined_results_fc  %>% 
    filter(padj < 0.01) %>%    # filter table
    
    pull(gene) %>%             # extract the gene column as a vector
    unique()


vsd_mat_tb <- vst_mat_T %>% data.frame() %>% rownames_to_column(var="rownames") %>% as_tibble()

# Check if the metadata samples match the column names in vsd_mat_tb
matching_samples <- intersect(metadata$samples, colnames(vsd_mat_tb)[-1])

# Filter metadata to include only matching samples
metadata_filtered <- metadata %>% filter(samples %in% matching_samples)

# Filter vsd_mat_tb to include only candidate genes
vsd_mat_tb_filtered <- vsd_mat_tb %>% filter(rownames %in% candidate_genes)


trans_cts_corr <- vsd_mat_tb_filtered %>% 
  # remove the column "gene", which we do not want to calculate correlation on
  select(-rownames) %>% 
  # we use Spearman's correlation, a non-parametric metric based on ranks
  cor(method = "spearman")


trans_cts_mean <- vsd_mat_tb_filtered %>%
    pivot_longer(cols = starts_with("T"), names_to = "samples", values_to = "cts") %>%
    full_join(metadata, by = c("rownames" = "samples")) %>%  # Adjust if necessary
    filter(rownames %in% candidate_genes) %>%
    group_by(rownames) %>%
    mutate(cts_scaled = (cts - mean(cts)) / sd(cts)) %>%
    group_by(rownames, condition) %>%
    summarise(Zscore_scaled_expression = mean(cts_scaled),
              pren = n()) %>%
    ungroup()

hclust_matrix <- vsd_mat_tb_filtered %>% 
  select(-rownames) %>% 
  as.matrix()
# assign rownames
rownames(hclust_matrix) <- vsd_mat_tb_filtered$rownames


hclust_matrix <- hclust_matrix %>% 
  # transpose the matrix so genes are as columns
  t() %>% 
  # apply scalling to each column of the matrix (genes)
  scale() %>% 
  # transpose back so genes are as rows again
  t()


#calculate the distance between each gene (row) 
gene_dist <- dist(hclust_matrix)

gene_hclust <- hclust(gene_dist, method = "complete")

trans_cts_cluster <- trans_cts_mean %>% 
    inner_join(clusters$df, by = c("rownames" = "gene")) %>%
    filter(cluster %in% survived_clusters)

# Ensure conditions are sorted
trans_cts_cluster$condition <- factor(trans_cts_cluster$condition, levels = c("T5min", "T15min", "T4h", "T8h", "T12h","T24h", "T48h", "T72h", "T120h", "T221h","T288h", "T290h", "T292h", "T294h", "T296h", "T296h_prim"))
# Count number of genes per cluster for labeling
cluster_sizes <- trans_cts_cluster %>%
    group_by(cluster) %>%
    summarise(num_genes = n_distinct(rownames))

# Create the plot
p8 <- ggplot(trans_cts_cluster, aes(x = condition, y = Zscore_scaled_expression, group = rownames)) +
    geom_line(alpha = 0.3, color = 'red') +
    geom_line(stat = "summary", fun = median, color = "blue", size = 1.5) +
    facet_wrap(~ cluster, nrow = 10, labeller = as_labeller(function(cluster) {
        num_genes <- cluster_sizes %>% filter(cluster == !!cluster) %>% pull(num_genes)
        paste("Cluster", cluster, "\n(", num_genes, ")")
    })) +
    theme_minimal() +
    theme(
        strip.background = element_rect(color = "black", size = 1.5, linetype = "solid"),
        strip.text.x = element_text(size = 8, color = "black", face = "bold"),
        legend.position = "bottom",
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title = element_text(size = 10),
        strip.text = element_text(size = 7)
    ) +
    labs(title = "Expression Profiles for Survived Clusters",
         x = "Condition",
         y = "Z-score of Scaled Expression")

plot_filepath <- file.path(output_directory, "Expreession profile of clusters.png")
ggsave(filename = plot_filepath, plot = p8, width = 10, height = 8)

