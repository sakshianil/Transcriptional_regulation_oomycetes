# Supporting Data Files

## Overview
The `support/` directory provides essential supporting data for the analyses carried out in this project. This includes reference genomes, annotations, effectors, transcription factor information, and conserved motif data used across three main studies. Each subfolder within `support/` corresponds to a specific study: `Pl_halstedii_study`, `orthologs_motifs_study`, and `orthologs_TF_study`.

## Directory Structure

### Pl_halstedii_study
This folder contains supporting documents, such as reference genomes and annotations, for *Plasmopara halstedii*. Specifically, it includes:
- `Han.fa`, `Han.gff`, `Han.gtf`: Genomic FASTA, GFF, and GTF files for *Helianthus annuus*.
- `Phals.fa`, `Phals.gff`, `Phals.gtf`: Genomic FASTA, GFF, and GTF files for *Plasmopara halstedii*.

These files were obtained from Ensembl Protists release 49. The genomic data for *Plasmopara halstedii* can be accessed [here](https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-49/protists/fasta/protists_stramenopiles1_collection/plasmopara_halstedii_gca_900000015/).

### orthologs_motifs_study
This folder contains supporting documents for the analysis of orthologs and motifs, including:
- Genomic FASTA, GFF, and GTF files for *Phytophthora infestans* and *Phytophthora sojae*, downloaded from Ensembl release 49.
- The `validated_effectors` folder, containing:
  - `effector_sequences.fa`: Effector sequences downloaded from the EffectorO GitHub repository, available [here](https://github.com/mjnur/oomycete-effector-prediction).
  - `effectors_blast_results.txt`: The output of BLASTp performed for effector sequences against ncbi non-redundant database, filtered for 100% identity. This file includes details such as query and subject IDs, percentage identity, alignment length, and more.
  - `validated_effectors_species.txt`: Contains validated effectors for species considered in this study.
- Conserved motif information, including RxLRs and CRNs motifs, used in the study. These motifs were obtained from *Sharma et al. 2015*.

### orthologs_TF_study
This folder contains supporting data related to transcription factor (TF) analysis in oomycetes, including:
- `FTFD_oomycetes.tsv`: A TSV file downloaded from the [Fungal Transcription Factor Database (FTFD)](http://ftfd.snu.ac.kr/download.php?a=list&o=Oomycota). The file provided is only small part of a large file.
- `JASPAR2020_CORE_non-redundant_pfms_meme.zip`: Transcription factor motifs in MEME format, downloaded from the [JASPAR 2020 database](https://jaspar2020.genereg.net/download/data/2020/CORE/JASPAR2020_CORE_non-redundant_pfms_meme.zip) or https://meme-suite.org/meme/db/motifs. The file provided is only small part of a large file.

## How to Use
The supporting data files provided here are used throughout the computational analyses and pipelines in this project. Each folder within `support/` contains specific data sets used for each of the three studies: 

1. **Pl_halstedii_study**: Contains reference genome and annotation files for *Plasmopara halstedii* and *Helianthus annuus*, essential for transcriptional analysis.
2. **orthologs_motifs_study**: Contains effector and conserved motif information, supporting the analysis of orthologs and conserved regions across *Phytophthora* species.
3. **orthologs_TF_study**: Provides TF-related data, including FTFD files and motifs from JASPAR, used for studying transcription factors across multiple oomycetes species.
## Notes
- The data in `validated_effectors` and motif information from *Sharma et al. 2015* are crucial for understanding effector biology in oomycetes.
- The processed effector sequences have been subjected to BLASTp analysis to identify identical sequences, with results saved in `effectors_blast_results.txt` for further interpretation.


