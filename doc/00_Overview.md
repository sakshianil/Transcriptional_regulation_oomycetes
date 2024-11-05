# Overview

## Project Title
**Transcriptional Regulation in Oomycetes**

## Project Description
This project aims to decipher transcriptional regulation in oomycetes using a combination of high-throughput RNA sequencing, computational analysis, and motif discovery approaches. Oomycetes, commonly known as water molds, include many agriculturally significant pathogens, and understanding their regulatory mechanisms can help in developing novel disease control methods.

This repository provides computational workflows and analytical scripts developed during the course of a PhD research project conducted in the Thines group. The research was a collaborative effort involving various institutions, including:

- **Senckenberg Biodiversity and Climate Research Centre (SBiK-F), Frankfurt (Main), Germany.**
- **Department of Biological Sciences, Institute of Ecology, Evolution, and Diversity, Goethe University, Frankfurt (Main), Germany.**
- **Center for Integrative Fungal Research (IPF), Frankfurt (Main), Germany.**
- **LOEWE Centre for Translational Biodiversity Genomics, Frankfurt (Main), Germany.**

The repository includes data processing scripts, analysis workflows, and visualization tools necessary to identify and analyze transcription factors, effector proteins, and regulatory motifs in oomycetes. The project aims to shed light on gene expression control in pathogenic oomycetes to better understand how these organisms interact with their hosts.

## Research Objectives
The primary objectives of the research are:

1. **Identify and classify transcription factors in oomycetes.**
2. **Analyze gene regulatory motifs across different species of oomycetes.**
3. **Investigate transcriptional changes during infection and identify key regulatory elements involved in pathogenicity.**

The repository is divided into three main sections based on the three studies undertaken during this research:

1. **Plasmopara halstedii Study**: Exploring the transcriptional regulation in the downy mildew pathogen *Plasmopara halstedii*.
2. **Orthologs Motifs Study**: Analysis of conserved motifs and effector proteins across oomycetes.
3. **Orthologs Transcription Factor Study**: Identification of conserved transcription factors across five oomycete genomes.
## Why Oomycetes?
Oomycetes are unique eukaryotic microorganisms that resemble fungi morphologically but are phylogenetically distinct. They are notorious plant pathogens, causing diseases such as late blight in potatoes and downy mildew in numerous crops. Understanding how oomycetes control gene expression can reveal critical insights into their lifecycle, host adaptation, and infection mechanisms.

## Technical Highlights

- **High-Throughput RNA-Seq Analysis**: Leveraging Illumina RNA sequencing data to study gene expression profiles.
- **Motif Discovery**: Identification of conserved motifs using MEME Suite and analysis of regulatory elements across different species.
- **Transcription Factor Classification**: Identification of transcription factors using the Fungal Transcription Factor Database (FTFD) and clustering based on functional domains.
- **Effector Protein Analysis**: Identifying key effector proteins and their motifs that play an important role in host-pathogen interaction.

## Workflow Summary
The computational workflows used for this research include:

1. **Pre-processing**: Using Trimmomatic to clean RNA-Seq reads and STAR to align them to reference genomes.
2. **Gene Quantification**: Using FeatureCounts for transcript quantification.
3. **Motif Analysis**: Utilizing tools from MEME Suite to find significant motifs and BLAST to identify orthologous sequences.
4. **Data Analysis and Visualization**: Scripts for statistical analysis (R) and data handling (Python), including downstream visualization of results.

The project is organized into well-defined modules, with each study having its own data, scripts, and results.

## Repository Structure
The repository is structured into folders representing different stages of the analysis:

- **data/**: Contains raw, processed, and supporting data files.
- **scripts/**: Contains all scripts for data processing, motif analysis, and downstream analysis. Subdivided into `python`, `R`, `shell`, and `perl`.
- **results/**: Outputs from the analysis, including processed data files, tables, and figures.
- **figures/**: Includes pipeline figures depicting the workflow for each study.
- **docs/**: Documentation files for understanding the project’s workflow, tools, and methodology.

## Technologies Used

- **Programming Languages**: Python, R, Perl, Bash.
- **Software**: Trimmomatic, STAR, FeatureCounts, MEME Suite, BLAST, Clustal Omega.
- **Data Analysis Libraries**: Pandas, NumPy, and various Bioconductor packages.

## Target Audience
This repository is aimed at researchers interested in transcriptomics, plant pathology, and bioinformatics. The workflows and scripts are intended for advanced users with some familiarity with command-line tools and bioinformatics analysis.

## Getting Started
To get started, please refer to the `README.md` in the root directory for general instructions on how to use the repository. Detailed instructions for running the workflows are provided in the `docs/` folder under individual markdown files.

## Acknowledgments
This research was supported by the Deutscher Akademischer Austauschdienst (DAAD) doctoral program and LOEWE Centre for Translational Biodiversity Genomics (TBG). The research would not have been possible without the generous support from these institutions.

For further inquiries, feel free to contact [Sakshi Bharti](mailto:sakshi.bharti@senckenberg.de).

