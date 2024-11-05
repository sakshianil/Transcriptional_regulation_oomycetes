# Methodology

## Overview

This section details the methodologies used in the project to explore transcriptional regulation in oomycetes. Each stage of the analysis—from data acquisition and processing to computational workflows—has been meticulously designed to unravel the complex interactions and regulatory mechanisms within oomycetes. The project follows a step-by-step pipeline that includes pre-processing, analysis, and validation.

## Experimental Design

The methodology is divided into three parts, each targeting a specific objective of the research project:

1. **Plasmopara halstedii Study**: Focusing on the downy mildew pathogen *Plasmopara halstedii*, this part of the project aimed at identifying the regulatory elements involved in infection mechanisms.
2. **Orthologs Motifs Study**: Analyzed effector protein motifs and regulatory motifs across multiple oomycetes to discover conserved motifs and elucidate the role of effectors in host-pathogen interactions.
3. **Orthologs Transcription Factor Study**: Investigated transcription factors conserved across five different oomycete genomes, emphasizing their role in shared pathogenic traits

## Data Collection

### Transcriptomic Data
High-throughput RNA sequencing data were collected from *Plasmopara halstedii* infected samples across different time points during infection. The RNA was extracted, purified, and sequenced using Illumina technology.

### Genome Data
Genomic data for oomycetes (including *Phytophthora infestans*, *Phytophthora sojae*, and others) were obtained from Ensembl Protists release 49. These genomes were used for motif discovery and effector protein analysis. Corresponding annotation files (GFF, GTF) were also retrieved for comprehensive analysis.

### Effector Sequences and Transcription Factor Information
- **Effector Sequences**: Conserved effector proteins were gathered from EffectorO GitHub repository.
- **Transcription Factor Information**: Oomycete-specific transcription factors were downloaded from the Fungal Transcription Factor Database (FTFD).

## Data Pre-processing

### Quality Control
The sequencing reads were first subjected to quality control using **Trimmomatic**. Poor quality reads, adapters, and low-quality bases were filtered out to ensure high-quality downstream analysis. Quality metrics were assessed using **FastQC**.

### Read Alignment
Reads were aligned to reference genomes using **STAR Mapper**. Indexing was performed beforehand to facilitate the alignment process, and the alignment outputs were saved in BAM format for further analysis.

## Feature Counting

Aligned reads were analyzed using **FeatureCounts** to quantify the expression levels of genes and features. The resulting counts files served as the basis for differential expression analysis.

### Count File Processing
The raw count files were processed for normalization and converted to formats compatible with downstream analysis, including R scripts for visualization and statistical tests.

## Differential Expression Analysis

The processed counts were analyzed for differential expression using **R** and specific Bioconductor packages. The analysis focused on finding genes with significantly different expression levels across various stages of the infection timeline.

## Motif Discovery and Analysis

### MEME Suite Analysis
To identify conserved motifs within the transcriptional regulatory regions, **MEME Suite** was employed. The motif discovery involved the following tools:

- **MEME**: To find de novo motifs across target sequences.
- **Tomtom**: To compare discovered motifs with known databases (e.g., JASPAR) to find matches and functional annotations.

### Motif Comparison Using JASPAR
The motifs identified were compared to motifs from the JASPAR core non-redundant database using Tomtom, helping in understanding functional similarities between oomycetes and other organisms.

### Effector Motif Analysis
Conserved effector motifs such as RxLRs and CRNs were analyzed using previously described motifs from Sharma et al., 2015. The effector motifs were searched against effector proteins to understand their conserved functionality in host-pathogen interaction.

## Transcription Factor Identification

### Fungal Transcription Factor Database (FTFD)
Transcription factors for each genome were identified using the **FTFD**. The identified TFs were categorized based on their family and domain information. Comparisons were made to highlight conserved transcription factors across the selected oomycetes.

## Functional Annotation and Analysis

### BlastP Search
**BlastP** was used to identify homologous sequences across species. Only hits with 100% identity were retained for further analysis to confirm orthologous relationships between transcription factors and effectors in different oomycetes.

## Effector Analysis

### Effector Protein Validation
Effector proteins were validated through **BlastP** against effector sequences from the EffectorO database. Only proteins with 100% identity were retained. Additionally, effector sequences were analyzed for conserved motifs relevant to their roles in pathogenicity.

## Data Visualization

The data were visualized using **R** packages for various plots, including heatmaps, volcano plots, and motif enrichment plots. Specific attention was given to:

- **Gene Expression Trends**: Plots to show differential expression across infection stages.
- **Motif Enrichment Plots**: Highlighting significantly enriched motifs and their conservation among species.
- **Effector Protein Distribution**: Visualizing how different effectors are conserved across different oomycete species.

## Reproducibility

All scripts and workflows used in this project are provided in the `scripts/` folder, organized by type (`python`, `R`, `shell`). Each script has detailed comments to guide users through the process, ensuring reproducibility. The documentation (`docs/` folder) provides further context and step-by-step instructions.

