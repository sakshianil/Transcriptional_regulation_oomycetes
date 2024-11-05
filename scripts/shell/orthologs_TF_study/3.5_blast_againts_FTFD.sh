#!/bin/bash

# Paths to input files
QUERY_FILE="/path_to/Transcriptional_regulation_oomycetes/results/python/orthologs_TF_study/output_query_sequences.fasta"  # Replace with the correct path to your query file
SUBJECT_FILE="/path_to/Transcriptional_regulation_oomycetes/results/python/orthologs_TF_study/output_subject_sequences.fasta"  # Replace with the correct path to your subject file
OUTPUT_FILE="/path_to/Transcriptional_regulation_oomycetes/results/shell/orthologs_TF_study/blastp_output.txt"  # Replace with the desired output path

# Make sure BLAST+ is installed and added to the PATH

# Create a BLAST database from the subject sequences
makeblastdb -in "$SUBJECT_FILE" -dbtype prot -out subject_db

# Run BLASTP with inferred parameters
blastp -query "$QUERY_FILE" \
       -db subject_db \
       -evalue 1e-5 \
       -outfmt "6 sseqid sseq sstart send qseqid qstart qend qseq length pident evalue" \
       -num_threads 4 \
       -out "$OUTPUT_FILE"

# Notify the user
echo "BLASTP search completed. Results saved in $OUTPUT_FILE."

