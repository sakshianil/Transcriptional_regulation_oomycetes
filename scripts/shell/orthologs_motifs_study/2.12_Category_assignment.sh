#!/bin/bash

# Define the input and output files
input_file="/path_to/Transcriptional_regulation_oomycetes/results/meme/orthologs_motifs_study/combined_fimo.tsv"  # Change to your actual input file name
output_file="/path_to/Transcriptional_regulation_oomycetes/results/shell/orthologs_motifs_study/categorized_combined_fimo.tsv" # The output file to create

# Create the header for the output file
echo -e "basename\tmotif_alt_id\tcategory\tmotif_id\tsequence_name\tstart\tstop\tstrand\tscore\tp-value\tq-value\tmatched_sequence" > "$output_file"

# Declare an associative array to track seen combinations of basename and motif_alt_id
declare -A seen_combinations

# Read the input file line by line
while IFS=$'\t' read -r basename motif_id motif_alt_id sequence_name start stop strand score p_value q_value matched_sequence
do
    # Skip the header line
    if [[ "$basename" == "basename" && "$motif_alt_id" == "motif_alt_id" ]]; then
        continue
    fi

    # Create a unique key from basename and motif_alt_id
    key="${basename}_${motif_alt_id}"

    # Check if this combination has been seen
    if [[ -n "${seen_combinations[$key]}" ]]; then
        continue  # Skip if this combination has been seen
    fi

    # Mark this combination as seen
    seen_combinations[$key]=1

    # Initialize category variable
    category="E" # Default to category E

    # Check against the regex for each genome
    match_pl_halstedii=$(echo "$sequence_name" | grep -E '^(CEG|CEJ|CEI|ENSRNAG)')
    match_ph_sojae=$(echo "$sequence_name" | grep -E '^(EPrPING|PHYSODRAFT_|ENSRNA|PHYSODRAFT_t)')
    match_ph_infestans=$(echo "$sequence_name" | grep -E '^(PITG_|EPrPING|ENSRNA)')

    # Debugging output to check matches
    echo "Matching for $sequence_name: Pl_halstedii: $match_pl_halstedii, Ph_sojae: $match_ph_sojae, Ph_infestans: $match_ph_infestans"

    # Determine categories based on matches
    if [[ -n "$match_pl_halstedii" && -n "$match_ph_sojae" && -n "$match_ph_infestans" ]]; then
        category="A" # Found in all three
    elif [[ -n "$match_pl_halstedii" && -n "$match_ph_sojae" ]]; then
        category="B" # Found in Pl. halstedii and Ph. sojae
    elif [[ -n "$match_ph_sojae" && -n "$match_ph_infestans" ]]; then
        category="C" # Found in Ph. sojae and Ph. infestans
    elif [[ -n "$match_pl_halstedii" && -n "$match_ph_infestans" ]]; then
        category="D" # Found in Pl. halstedii and Ph. infestans
    fi

    # Append the result to the output file
    echo -e "$basename\t$motif_alt_id\t$category\t$motif_id\t$sequence_name\t$start\t$stop\t$strand\t$score\t$p_value\t$q_value\t$matched_sequence" >> "$output_file"
done < "$input_file"

echo "Processing complete. Output saved to $output_file."

