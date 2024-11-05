# Data Documentation

## Overview

The data used in this project consists of several datasets, which include genomic, transcriptomic, and effector-related information from various oomycetes. The project aims to unravel transcriptional regulation in these organisms by leveraging data from multiple sources, including publicly available databases and experimentally generated sequences.

## Data Types and Sources

The datasets are divided into three main categories based on the objectives of the project:

1. **Plasmopara halstedii Study**
   - **Genomic and Annotation Data**: The genomic sequences and corresponding annotation files (FASTA, GTF, GFF) for *Plasmopara halstedii* were obtained from Ensembl Protists release 49.
   - **RNA-Seq Data**: High-throughput RNA sequencing data from *Plasmopara halstedii* were collected at multiple time points during infection. This transcriptomic data was used to identify gene expression trends throughout different stages of infection. All raw data used in this study are available from GenBank under the accession number [PRJEB49134](https://www.ncbi.nlm.nih.gov/bioproject/49134).

2. **Orthologs Motifs Study**
   - **Effector Sequences**: Conserved effector protein sequences across multiple oomycete species were sourced from the EffectorO GitHub repository. These effectors were further analyzed for known motifs, such as RxLRs and CRNs, as described by Sharma et al., 2015.

3. **Orthologs Transcription Factor Study**
   - **Oomycete Transcription Factors**: Transcription factor data for oomycetes were downloaded from the Fungal Transcription Factor Database (FTFD). The dataset includes information on TF families, TF names, functional annotations, and the species in which each TF was identified. This dataset provided a comparative foundation for identifying conserved transcription factors among the selected oomycete genomes.
   - **JASPAR Motif Database**: To perform motif discovery and motif comparison analysis, we used the non-redundant motif database provided by JASPAR (https://jaspar2020.genereg.net/download/data/2020/CORE/JASPAR2020_CORE_non-redundant_pfms_meme.zip).

## Data Organization

All datasets are organized in the `data/` directory, which is divided into three subdirectories:

- **raw/**: Contains the raw data files, including RNA sequencing files, genome data, and TF datasets, as originally collected from public databases.
- **processed/**: Contains the processed datasets that are ready for downstream analysis. This includes normalized count matrices, filtered transcription factor datasets, and curated effector sequences.
- **support/**: Holds supporting data files necessary for interpreting the analysis, including:
  - Annotation files (FASTA, GTF, GFF) for all species studied.
  - Validated effector sequences from the EffectorO repository.
  - Transcription factor information for oomycetes from the FTFD.

## Data Sources

1. **Ensembl Protists (release 49)**: Genomic and annotation files for *Plasmopara halstedii*, *Phytophthora infestans*, *Phytophthora sojae*, and others were obtained from Ensembl's public FTP.
2. **Fungal Transcription Factor Database (FTFD)**: The oomycete-specific transcription factor dataset was downloaded from [FTFD's official website](http://ftfd.snu.ac.kr/download.php?a=list&o=Oomycota).
3. **EffectorO GitHub Repository**: Effector protein sequences were gathered from the EffectorO repository, available at [EffectorO GitHub](https://github.com/mjnur/oomycete-effector-prediction).
4. **JASPAR Database**: Motifs from the JASPAR database were downloaded from [JASPAR's official site](https://jaspar2020.genereg.net/download/data/2020/CORE/JASPAR2020_CORE_non-redundant_pfms_meme.zip).
5. **Sharma et al., 2015**: Conserved effector motifs such as RxLRs and CRNs were referenced from Sharma et al.'s work on effector proteins in oomycetes.

## Processed Data Files

The `processed/` subdirectory contains the data files that have been curated and prepared for further analysis. The key files include:

- **Gene Counts Matrix**: Normalized count files generated from RNA-Seq data, capturing the gene expression levels across multiple time points of infection.
- **Effector Sequences**: Curated effector sequences that were retained based on high sequence identity with known effector proteins.
- **Motif Analysis Outputs**: Processed results from motif discovery, including MEME,STREME output files and motif matches from JASPAR.

## Data File Overview

- **final_counts_path.txt**: Contains gene expression counts for *Plasmopara halstedii*, organized by time points (e.g., T5min, T15min, T4h, etc.). Each row represents a unique gene, with corresponding expression levels.
- **final_counts_path_with_description.txt**: An extended version of the `final_counts_path.txt` that includes additional annotations, such as gene IDs, description, and genomic locations.
- **Effector Sequences (`effector_sequences.fa`)**: This file in the `validated_effectors` folder contains sequences of validated effector proteins used in subsequent motif analysis.
- **BlastP Results (`effectors_blast_results.txt`)**: Contains the results of BlastP comparisons of validated effectors, listing only sequences with 100% identity to known effectors in non-redundant ncbi database.

## Notes on Data Usage

1. **Privacy and Ethical Use**: The genomic data used in this project is publicly available and should be used according to the respective data provider's usage policies.
2. **File Formats**: Data is provided in standard bioinformatics formats like FASTA, GTF, GFF, and TSV for ease of use across a wide range of bioinformatics tools.
3. **Accessibility**: Users can find the raw, processed, and supporting data in the respective subdirectories within `data/` to perform custom analyses or reproduce results as required.

## Contribution

Users are encouraged to use, analyze, and further process the data provided in this repository. Please acknowledge the original sources and this repository in any publication or presentation based on these datasets.
