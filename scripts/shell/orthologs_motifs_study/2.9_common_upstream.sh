#!/bin/bash

# Define paths to the gene, promoter, and output directories
gene_folder="/path_to/Transcriptional_regulation_oomycetes/results/python/orthologs_motifs_study/clusters/gene"
promoter_folder="/path_to/Transcriptional_regulation_oomycetes/results/python/orthologs_motifs_study/clusters/promoter"
common_upstream_folder="/path_to/Transcriptional_regulation_oomycetes/results/python/orthologs_motifs_study/clusters/common_upstream"

# Create the output folder if it doesn't exist
mkdir -p "$common_upstream_folder"

# Loop through all gene files in the gene folder
for gene_file in "$gene_folder"/*.txt; do
    # Get the basename of the gene file (e.g., "1.txt")
    gene_filename=$(basename "$gene_file")
    
    # Check if the corresponding promoter file exists
    promoter_file="$promoter_folder/$gene_filename"
    if [ -f "$promoter_file" ]; then
        # Output file path
        output_file="$common_upstream_folder/$gene_filename"
        
        # Declare an associative array to store promoter lines for fast lookup
        declare -A promoter_hash
        
        # Read the promoter file and populate the hash table
        while read -r promoter_line; do
            # Extract query and subject IDs (without strand information)
            query_subject_ids=$(echo "$promoter_line" | awk -F'\t' '{gsub(/\(.*\)/,"",$1); gsub(/\(.*\)/,"",$2); print $1,$2}')
            promoter_hash["$query_subject_ids"]="$promoter_line"
        done < "$promoter_file"
        
        # Compare gene lines against the promoter hash
        while read -r gene_line; do
            # Extract query and subject IDs from the gene line (without strand information)
            query_subject_ids=$(echo "$gene_line" | awk -F'\t' '{gsub(/\(.*\)/,"",$1); gsub(/\(.*\)/,"",$2); print $1,$2}')
            
            # Check if the IDs exist in the promoter hash
            if [[ -n "${promoter_hash["$query_subject_ids"]}" ]]; then
                # If matched, append the corresponding promoter line to the output file
                echo "${promoter_hash["$query_subject_ids"]}" >> "$output_file"
            fi
        done < "$gene_file"
        
        # Clear the hash table for the next iteration
        unset promoter_hash
    else
        echo "No corresponding promoter file for $gene_filename"
    fi
done

echo "Common sequences saved in the folder: $common_upstream_folder"

