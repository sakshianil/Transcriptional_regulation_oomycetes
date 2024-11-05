#!/bin/bash

# Set up email for NCBI (if needed)
export BLASTDB="/path_to/Transcriptional_regulation_oomycetes/results/python/orthologs_motifs_study" # Specify your BLAST database directory
export PATH="$PATH:/path_to/ncbi-blast-2.11.0+/bin" # Ensure BLAST is in your PATH. Replace path_to with local path

# Define the paths to your cluster files
GENE_CLUSTER_PATH="/path_to/Transcriptional_regulation_oomycetes/results/python/orthologs_motifs_study/clusters/gene"
PROMOTER_CLUSTER_PATH="path_to/Transcriptional_regulation_oomycetes/results/python/orthologs_motifs_study/clusters/promoter"

# Define output database names
GENE_DB_NAME="gene_sequences_db"
PROMOTER_DB_NAME="promoter_sequences_db"

# Create directories for output (if they don't already exist)
mkdir -p "$BLASTDB/gene"
mkdir -p "$BLASTDB/promoter"

# Combine all gene FASTA files into one file
echo "Combining gene cluster FASTA files..."
cat "$BLASTDB/gene"/*.fasta > "$BLASTDB/gene/combined_genomes.fasta"

# Combine all promoter FASTA files into one file
echo "Combining promoter cluster FASTA files..."
cat "$BLASTDB/promoter"/*.fasta > "$BLASTDB/promoter/combined_promoters.fasta"

# Create BLAST databases
echo "Creating gene sequences BLAST database..."
makeblastdb -in "$BLASTDB/gene/combined_genomes.fasta" -dbtype nucl -out "$BLASTDB/gene/$GENE_DB_NAME"

echo "Creating promoter sequences BLAST database..."
makeblastdb -in "$BLASTDB/promoter/combined_promoters.fasta" -dbtype nucl -out "$BLASTDB/promoter/$PROMOTER_DB_NAME"

# Run BLAST for each gene cluster file against gene database
echo "Running BLAST for gene sequences against gene sequences..."
for gene_file in "$GENE_CLUSTER_PATH"/*.fasta; do
    cluster_name=$(basename "$gene_file" .fasta)  # Get the base name of the cluster file
    blastn -query "$gene_file" -db "$BLASTDB/gene/$GENE_DB_NAME" -perc_identity 50 -evalue 0.01 -outfmt "6 qseqid sseqid qstart qend qseq sstart send sseq length pident evalue" -out "$BLASTDB/gene/${cluster_name}.txt"
done

# Run BLAST for each promoter cluster file against promoter database
echo "Running BLAST for promoter sequences against promoter sequences..."
for promoter_file in "$PROMOTER_CLUSTER_PATH"/*.fasta; do
    cluster_name=$(basename "$promoter_file" .fasta)  # Get the base name of the cluster file
    blastn -query "$promoter_file" -db "$BLASTDB/promoter/$PROMOTER_DB_NAME" -perc_identity 70 -word_size 11 -outfmt "6 qseqid sseqid qstart qend qseq sstart send sseq length pident evalue" -out "$BLASTDB/promoter/${cluster_name}.txt"
done

echo "BLAST searches completed. Results saved in $BLASTDB/gene and $BLASTDB/promoter directories."

