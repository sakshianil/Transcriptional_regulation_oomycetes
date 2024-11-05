#!/bin/bash

# Set paths to directories and executables
clusters_fasta="/path_to/Transcriptional_regulation_oomycetes/results/python/orthologs_motifs_study/clusters/commom_upstream/fasta_output"
jaspar_db="/path_to/motif_databases/JASPAR/JASPAR2020_CORE_vertebrates_non-redundant.meme"
meme_output="/path_to/Transcriptional_regulation_oomycetes/results/meme/orthologs_motif_study"
tomtom="/path_to/meme/bin/tomtom"
# Debugging: Print paths to verify they are correct
echo "Clusters FASTA path: $clusters_fasta"
echo "MEME output directory: $meme_output"

# Ensure output directories exist
mkdir -p "$meme_output"


# Run TOMTOM on MEME output files
for meme_out_dir in "$clusters_fasta"/*_meme_out; do
    if [ -d "$meme_out_dir" ]; then
        meme_txt_file="$meme_out_dir/meme.txt"
        if [ -f "$meme_txt_file" ]; then
            # Extract base name without the _meme_out suffix
            base_name=$(basename "$meme_out_dir" "_meme_out")
            # Construct the path to the corresponding FASTA file
            fasta_file="$clusters_fasta/${base_name}.fasta"
            echo "Checking FASTA file: $fasta_file"  # Debug output
            if [ -f "$fasta_file" ]; then
                tomtom_out_dir="${meme_out_dir}/tomtom_output"
                mkdir -p "$tomtom_out_dir"
                $tomtom --oc "$tomtom_out_dir"  -min-overlap 5 -dist pearson -evalue -thresh 0.5 -no-ssc "$meme_out_dir/meme.txt" "$jaspar_db" 
            else
                echo "Warning: Corresponding FASTA file not found for $meme_out_dir"
            fi
        else
            echo "Warning: meme.txt file not found for $meme_out_dir"
        fi
    fi
done

# Run TOMTOM on STREME output files
for streme_out_dir in "$clusters_fasta"/*_streme_out; do
    if [ -d "$streme_out_dir" ]; then
        streme_txt_file="$streme_out_dir/streme.txt"
        if [ -f "$streme_txt_file" ]; then
            # Extract base name without the _streme_out suffix
            base_name=$(basename "$streme_out_dir" "_streme_out")
            # Construct the path to the corresponding FASTA file
            fasta_file="$clusters_fasta/${base_name}.fasta"
            echo "Checking FASTA file: $fasta_file"  # Debug output
            if [ -f "$fasta_file" ]; then
                tomtom_out_dir="${streme_out_dir}/tomtom_output"
                mkdir -p "$tomtom_out_dir"
                $tomtom --oc "$tomtom_out_dir" -min-overlap 5 -dist pearson -evalue -thresh 0.5 -no-ssc "$streme_out_dir/streme.txt" "$jaspar_db"
            else
                echo "Warning: Corresponding FASTA file not found for $streme_out_dir"
            fi
        else
            echo "Warning: streme.txt file not found for $streme_out_dir"
        fi
    fi
done

END
# Initialize the combined FIMO file with the header
combined_tomtom_file="$meme_output/combined_tomtom.tsv"

# Function to add basename as cluster number column
add_basename_column() {
    local input_file="$1"
    local base_name="$2"
    # Remove comment lines and add basename as cluster number
    awk -v bn="$base_name" 'NR>1 && $0 !~ /^#/ {print bn "\t" $0}' "$input_file"
}

# Get the header from the first non-empty FIMO file
first_tomtom_file=$(find "$clusters_fasta" -type f -path "*/tomtom_output/tomtom.tsv" | head -n 1)
if [ -f "$first_tomtom_file" ]; then
    # Extract header from the first non-commented line of the first file
    grep -v '^#' "$first_tomtom_file" | head -n 1 | awk '{print "basename\t" $0}' > "$combined_tomtom_file"
else
    echo "No TOMTOM files found for header extraction. Exiting."
    exit 1
fi

# Process and merge TOMTOM results from MEME output directories
for meme_out_dir in "$clusters_fasta"/*_meme_out; do
    if [ -d "$meme_out_dir" ]; then
        base_name=$(basename "$meme_out_dir" "_meme_out")
        tomtom_out_dir="$meme_out_dir/tomtom_output"
        if [ -d "$tomtom_out_dir" ]; then
            for tomtom_file in "$tomtom_out_dir"/tomtom.tsv; do
                if [ -f "$tomtom_file" ] && [ -s "$tomtom_file" ]; then  # Check if file is not empty
                    echo "Processing file: $tomtom_file"  # Debug output
                    add_basename_column "$tomtom_file" "$base_name" >> "$combined_tomtom_file"
                else
                    echo "Skipping empty or non-existent file: $tomtom_file"  # Debug output
                fi
            done
        fi
    fi
done


# Process and merge TOMTOM results from STREME output directories
for streme_out_dir in "$clusters_fasta"/*_streme_out; do
    if [ -d "$streme_out_dir" ]; then
        base_name=$(basename "$streme_out_dir" "_streme_out")
        tomtom_out_dir="$streme_out_dir/tomtom_output"
        if [ -d "$tomtom_out_dir" ]; then
            for tomtom_file in "$tomtom_out_dir"/tomtom.tsv; do
                if [ -f "$tomtom_file" ] && [ -s "$tomtom_file" ]; then  # Check if file is not empty
                    echo "Processing file: $tomtom_file"  # Debug output
                    add_basename_column "$tomtom_file" "$base_name" >> "$combined_tomtom_file"
                else
                    echo "Skipping empty or non-existent file: $tomtom_file"  # Debug output
                fi
            done
        fi
    fi
done

echo "Combined FIMO results saved to $combined_tomtom_file"


