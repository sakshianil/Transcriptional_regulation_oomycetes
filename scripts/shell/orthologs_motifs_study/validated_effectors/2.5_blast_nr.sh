#!/bin/bash

# Define input and output files
QUERY="/path_to/Transcriptional_regulation_oomycetes/data/support/orthologs_motifs_study/validated_effectors/effector_sequences.fa"  # Replace with the path to query FASTA file
DATABASE="nr_db"                  # Replace with your BLAST database path, e.g., nr (Non-redundant)
OUTFILE="/path_to/Transcriptional_regulation_oomycetes/data/support/orthologs_motifs_study/validated_effectors/effectors_blast_results.txt"     # Output file for BLAST results

# Run standalone BLAST with 100% identity filter
blastp \
-query $QUERY \
-db $DATABASE \
-out $OUTFILE \
-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" \
-evalue 1e-5 \
-max_target_seqs 100 \
-perc_identity 100 \
-gapopen 0 \
-gapextend 0

echo "BLAST search completed. Results stored in $OUTFILE."

