# Python Scripts for Transcriptional Regulation Analysis

## Overview
This folder contains various Python scripts designed for data processing tasks involved in the transcriptional regulation project. These scripts include functions for extracting motifs, manipulating datasets, integrating different biological data, and more. Each script focuses on a specific aspect of the analysis, allowing efficient and reproducible workflows. The structure below categorizes scripts according to different studies conducted within the project.

### Directory Structure
- **`Pl_halstedii_study/`**: Scripts used for analyzing **Plasmopara halstedii**.
  - **`1.3_Bg_promoters_search_phal.py`**:
    - **Purpose**: Searches for background promoters in *Plasmopara halstedii*. This script identifies the upstream promoter regions for the genome, extracting them for subsequent analysis.
    - **Output**: The output file is **`bg_promoter_100.fasta`**, which serves as the input for the Perl script **`fasta-get-markov.pl`** to generate background models for motif analysis.
  
  - **`1.4_clusters_Promoters_search_coord.py`**:
    - **Purpose**: Searches for promoter regions based on the coordinates of gene clusters identified by the R script **`1.2_DGE_analysis.R`**.
    - **Output**: The results are stored in **`results/python/Pl_halstedii_study/clusters_fasta/`** directory, where each output file represents a specific cluster (e.g., `output_1` for cluster 1).

- **`orthologs_motif_study/`**: Scripts for orthologous motif analysis.
  - **`2.1_clusters_gene_extractor.py`**:
    - **Purpose**: Extracts gene sequences based on cluster information from **`1.2_DGE_analysis.R`**. The output is used as the query for local BLAST searches against a gene database.
    
  - **`2.2_promoter_extractor_from_genome.py`**:
    - **Purpose**: Extracts promoter sequences (upstream 400 nucleotides from the start codon) for three genomes used in the study.
    - **Output**: The results are saved in **`results/python/orthologs_motifs_study/bg_promoter_full_400.fasta`**. This file is used as a subject database for local BLAST searches.
  
  - **`2.3_clusters_promoter_extractor.py`**:
    - **Purpose**: Extracts promoter sequences for gene clusters using the cluster information generated by **`1.2_DGE_analysis.R`**.
    - **Output**: These promoter sequences are used as queries in local BLAST analysis.
  
  - **`2.4_gene_extractor_from_genome.py`**:
    - **Purpose**: Extracts genes from the genome. This script is used for all three genomes, with respective file names and paths adjusted in the script.
    - **Output**: These extracted genes serve as the local BLAST gene database.
  
  - **`validated_effector_genes/`**: Scripts for analyzing validated effector genes.
    - **`2.6_fetch_Gene_from_proteinids.py`**:
      - **Purpose**: Fetches gene IDs and sequences based on protein IDs, provided in the file `effectors_blast_results.txt`.
      - **Output**: This output is merged with **`gene_extractor_from_genome.py`** results to create the local BLAST gene database.
    
    - **`2.7_fetch_upstream_from_geneid.py`**:
      - **Purpose**: Extracts upstream regions based on gene IDs generated by **`fetch_Gene_from_proteinids_editted_genesonly.py`**.
      - **Output**: The generated FASTA file, **`singleline_upstream_sequences_with_coordinates.fasta`**, is used in combination with **`promoter_extractor_from_genome.py`** for subsequent analysis.

- **`orthologs_TF_study/`**: Scripts for transcription factor analysis.
  - **`3.2_extract_jaspar_info.py`**:
    - **Purpose**: Extracts motif information from the JASPAR database based on results generated by the shell script **`3.1_motif_analysis_against_jaspar.sh`**.
  
  - **`3.3_extract_protein_info.py`**:
    - **Purpose**: Extracts protein sequence information from UniProt based on IDs in **`combined_with_profile_summaries.tsv`**.
    - **Output**: Generates the query protein sequences file **`output_query_sequences.fasta`**.
  
  - **`3.4_extract_protein_info_ftfd.py`**:
    - **Purpose**: Extracts protein sequences from the Fungal Transcription Factor Database (FTFD) based on alignment with **`output_query_sequences.fasta`**.
    - **Output**: The extracted sequences are saved in **`output_subject_sequences.fasta`**.

## How to Use These Scripts

1. **Dependencies**: Ensure you have the necessary dependencies installed for each script, such as `Biopython`, `pandas`, and other Python libraries.

2. **Input Files**: Each script has specific input files that need to be properly formatted. Details about the input files can be found within the comments at the beginning of each script.

3. **Execution**: To execute the scripts, navigate to the appropriate directory and run the scripts using Python 3. For example:
   ```bash
   cd python/Pl_halstedii_study/
   python3 1.3_Bg_promoters_search_phal.py
   ```

4. **Output**: All outputs are saved in the **`results/`** directory, with each subfolder corresponding to the specific analysis conducted. The naming convention for files should help identify the related input script and data type.

## Important Notes
- **Reverse Complements**: Please note that for scripts involving upstream regions, the reverse complements for sequences on the negative strand have already been handled by the Python scripts.

- **Integration with Other Tools**: Several Python scripts generate outputs that serve as inputs for other tools, such as the Perl script **`fasta-get-markov.pl`** for background modeling and various shell scripts for BLAST and motif analysis.

- **Folder Organization**: Scripts are categorized according to different studies for ease of navigation. Ensure that inputs and outputs are correctly specified and placed in the expected directories to maintain the workflow.

## Contact
For any questions or clarifications regarding the usage of these Python scripts, please contact [Sakshi Bharti](mailto:sakshi.bharti@senckenberg.de).


