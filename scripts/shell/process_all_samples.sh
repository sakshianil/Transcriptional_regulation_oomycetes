#!/bin/bash

# Define the raw data directory
RAW_DIR="../../data/raw"

# Loop through each file in the raw data directory
for FILE in $RAW_DIR/*.fq.gz; do
    # Extract the base name of the file (e.g., R5min_2_1 from R5min_2_1.fq.gz)
    BASENAME=$(basename "$FILE" .fq.gz)
    
    # Determine the organism (host or pathogen) based on some naming convention
    if [[ "$BASENAME" == *"host"* ]]; then
        ORGANISM="host"
    else
        ORGANISM="pathogen"
    fi

    # Print the command for user's information (optional)
    echo "Processing $FILE for $ORGANISM"

    # Run the combined_trimming_and_mapping.sh script
    ./trimming_and_mapping.sh "$ORGANISM" "$FILE"
done

