#!/bin/bash

# Input and output directories
input_folder="/path_to/Transcriptional_regulation_oomycetes/results/python/orthologs_motifs_study/clusters/common_upstream"
output_folder="/path_to/Transcriptional_regulation_oomycetes/scripts/python/orthologs_motifs_study/clusters/common_upstream/fasta_output"

# Create the output directory if it doesn't exist
mkdir -p "$output_folder"

# Loop through all text files in the input folder
for input_file in "$input_folder"/*.txt; do
    output_file="$output_folder/$(basename "$input_file" .txt).fasta"

    # Read through the input file line by line
    while read -r line; do
        # Extract fields using awk
        query=$(echo "$line" | awk '{print $1}')
        subject=$(echo "$line" | awk '{print $2}')
        start_query=$(echo "$line" | awk '{print $3}')
        end_query=$(echo "$line" | awk '{print $4}')
        seq_query=$(echo "$line" | awk '{print $5}')
        start_subject=$(echo "$line" | awk '{print $6}')
        end_subject=$(echo "$line" | awk '{print $7}')
        seq_subject=$(echo "$line" | awk '{print $8}')

	# Remove any characters from sequences that are not letters (A-Z or a-z)
        seq_query_clean=$(echo "$seq_query" | tr -d '-')
        seq_subject_clean=$(echo "$seq_subject" | tr -d '-')

        # Format the FASTA headers with gene name and coordinates only (desired format)
        query_header=">${query}(${start_query}-${end_query})"
        subject_header=">${subject}(${start_subject}-${end_subject})"

        # Output the FASTA formatted sequences to the output file
        echo "$query_header" >> "$output_file"
        echo "$seq_query" >> "$output_file"
        echo "$subject_header" >> "$output_file"
        echo "$seq_subject" >> "$output_file"
    done < "$input_file"
done

echo "FASTA files with desired format have been generated in $output_folder"

