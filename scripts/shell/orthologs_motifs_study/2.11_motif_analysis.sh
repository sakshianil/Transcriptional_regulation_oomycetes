#!/bin/bash

# Set paths to directories and executables
bg_model="/path_to/Transcriptional_regulation_oomycetes/results/python/orthologs_motifs_study/bg_promoter_only_aligned.fasta"
bg_model_meme="/path_to/Transcriptional_regulation_oomycetes/results/python/orthologs_motifs_study/bg_promoter_full_400.fasta"
clusters_fasta="/path_to/Transcriptional_regulation_oomycetes/results/python/orthologs_motif_study/clusters/common_upstream/fasta_output"
meme_output="/path_to/Transcriptional_regulation_oomycetes/results/meme/orthologs_motif_study"
meme="/path_to/meme/bin/meme"
perl_script="/path_to/Transcriptional_regulation_oomycetes/scripts/perl"
streme="/path_to/meme-5.3.0/src/streme" 
fimo="/path_to/meme-5.3.0/src/fimo"
gomo="/path_to/meme-5.3.0/src/gomo"
ame="/path_to/meme-5.3.0/src/ame"
# Debugging: Print paths to verify they are correct
echo "Background model path: $bg_model"
echo "Clusters FASTA path: $clusters_fasta"
echo "MEME output directory: $meme_output"
echo "Perl script directory: $perl_script"
echo "STREME executable: $streme"

# Ensure output directories exist
mkdir -p "$meme_output"

# Generate a Markov model for the background sequences
perl "$perl_script/fasta-get-markov.pl" -m 3 < "$bg_model" > "$meme_output/promoter_400.model"

perl "$perl_script/fasta-get-markov.pl" -m 3 < "$bg_model_meme" > "$meme_output/promoter_400_meme.model"

# Check if the Markov model was successfully generated
if [[ -f "$meme_output/promoter_400.model" ]]; then
    echo "Markov model generated at: $meme_output/promoter_400.model"
else
    echo "Error: Markov model was not generated. Check the paths and Perl script."
fi

# Check if the Markov model for meme was successfully generated
if [[ -f "$meme_output/promoter_400_meme.model" ]]; then
    echo "Markov model generated at: $meme_output/promoter_400_meme.model"
else
    echo "Error: Markov model for meme was not generated. Check the paths and Perl script."
fi

$streme -oc $meme_output/streme_400_output_files -p $bg_model --kmer 4


$fimo -oc $meme_output/fimo_streme_100_output_files --norc --parse-genomic-coord $meme_output/streme_100_output_files/streme.txt $bg_model

: << 'END'
for i in "$clusters_fasta"/*.fasta; do
    filtered_file="${i%.fasta}_filtered.fasta"
    awk -v n=5 '/^>/{ if(l>n) print b; b=$0;l=0;next }{l+=length;b=b ORS $0}END{if(l>n) print b }' "$i" > "$filtered_file"
done
END
# Run MEME for motif discovery on filtered FASTA files
for i in "$clusters_fasta"/*.fasta; do
    output_dir="${i%.fasta}_meme_out"
    $meme "$i" -mod zoops -objfun de -test mhg -dna -minw 4 -maxw 16 -nmotifs 10 -evt 0.5 \
    -bfile "$meme_output/promoter_400_meme.model" -markov_order 3 \
    -oc "$output_dir"
done

for i in "$clusters_fasta"/*.fasta; do
    if awk '/^>/{ if(l > 50) exit 0 } {l+=length} END{exit l > 50 ? 0 : 1}' "$i"; then
        $streme -oc "${i%.fasta}_streme_out" -p "$i" --kmer 4 --niter 4000 --seed 10
    else
        echo "No sequences longer than 50 nucleotides in $i, skipping STREME."
    fi
done


# Run FIMO on MEME output files
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
                fimo_out_dir="${meme_out_dir}/fimo_output"
                mkdir -p "$fimo_out_dir"
                $fimo -oc "$fimo_out_dir" --norc --parse-genomic-coord "$meme_txt_file" "$fasta_file"
            else
                echo "Warning: Corresponding FASTA file not found for $meme_out_dir"
            fi
        else
            echo "Warning: meme.txt file not found for $meme_out_dir"
        fi
    fi
done

# Run FIMO on STREME output files
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
                fimo_out_dir="${streme_out_dir}/fimo_output"
                mkdir -p "$fimo_out_dir"
                $fimo -oc "$fimo_out_dir" --norc --parse-genomic-coord "$meme_output/streme_400_output_files/streme.txt" "$fasta_file"
            else
                echo "Warning: Corresponding FASTA file not found for $streme_out_dir"
            fi
        else
            echo "Warning: streme.txt file not found for $streme_out_dir"
        fi
    fi
done


# Initialize the combined FIMO file with the header
combined_fimo_file="$meme_output/combined_fimo.tsv"

# Function to add basename as cluster number column
add_basename_column() {
    local input_file="$1"
    local base_name="$2"
    # Remove comment lines and add basename as cluster number
    awk -v bn="$base_name" 'NR>1 && $0 !~ /^#/ {print bn "\t" $0}' "$input_file"
}

# Get the header from the first non-empty FIMO file
first_fimo_file=$(find "$clusters_fasta" -type f -path "*/fimo_output/fimo.tsv" | head -n 1)
if [ -f "$first_fimo_file" ]; then
    # Extract header from the first non-commented line of the first file
    grep -v '^#' "$first_fimo_file" | head -n 1 | awk '{print "basename\t" $0}' > "$combined_fimo_file"
else
    echo "No FIMO files found for header extraction. Exiting."
    exit 1
fi

# Process and merge FIMO results from MEME output directories
for meme_out_dir in "$clusters_fasta"/*_meme_out; do
    if [ -d "$meme_out_dir" ]; then
        base_name=$(basename "$meme_out_dir" "_meme_out")
        fimo_out_dir="$meme_out_dir/fimo_output"
        if [ -d "$fimo_out_dir" ]; then
            for fimo_file in "$fimo_out_dir"/fimo.tsv; do
                if [ -f "$fimo_file" ] && [ -s "$fimo_file" ]; then  # Check if file is not empty
                    echo "Processing file: $fimo_file"  # Debug output
                    add_basename_column "$fimo_file" "$base_name" >> "$combined_fimo_file"
                else
                    echo "Skipping empty or non-existent file: $fimo_file"  # Debug output
                fi
            done
        fi
    fi
done

# Process and merge FIMO results from STREME output directories
for streme_out_dir in "$clusters_fasta"/*_streme_out; do
    if [ -d "$streme_out_dir" ]; then
        base_name=$(basename "$streme_out_dir" "_streme_out")
        fimo_out_dir="$streme_out_dir/fimo_output"
        if [ -d "$fimo_out_dir" ]; then
            for fimo_file in "$fimo_out_dir"/fimo.tsv; do
                if [ -f "$fimo_file" ] && [ -s "$fimo_file" ]; then  # Check if file is not empty
                    echo "Processing file: $fimo_file"  # Debug output
                    add_basename_column "$fimo_file" "$base_name" >> "$combined_fimo_file"
                else
                    echo "Skipping empty or non-existent file: $fimo_file"  # Debug output
                fi
            done
        fi
    fi
done

echo "Combined FIMO results saved to $combined_fimo_file"


# Run AME on MEME output files
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
                ame_out_dir="${meme_out_dir}/ame_output"
                mkdir -p "$ame_out_dir"
                $ame --oc "$ame_out_dir"  ---scoring avg -method fisher --evalue-report-threshold 0.5 "$fasta_file" "$meme_output/streme_400_output_files/streme.txt" 
            else
                echo "Warning: Corresponding FASTA file not found for $meme_out_dir"
            fi
        else
            echo "Warning: meme.txt file not found for $meme_out_dir"
        fi
    fi
done

# Run AME on STREME output files
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
                ame_out_dir="${streme_out_dir}/ame_output"
                mkdir -p "$ame_out_dir"
                $ame --oc "$ame_out_dir" --scoring avg --method fisher --evalue-report-threshold 0.5 "$fasta_file" "$meme_output/streme_400_output_files/streme.txt"  
            else
                echo "Warning: Corresponding FASTA file not found for $streme_out_dir"
            fi
        else
            echo "Warning: streme.txt file not found for $streme_out_dir"
        fi
    fi
done
