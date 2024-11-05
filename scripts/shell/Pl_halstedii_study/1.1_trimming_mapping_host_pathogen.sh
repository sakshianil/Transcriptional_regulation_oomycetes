#!/bin/bash

# Function to print usage
usage() {
  echo "Usage: $0 -r <raw_data_dir> -o <output_dir>"
  exit 1
}

# Parse command line arguments
while getopts ":r:o:" opt; do
  case ${opt} in
    r )
      RAW_DIR=$OPTARG
      ;;
    o )
      OUTPUT_DIR=$OPTARG
      ;;
    \? )
      usage
      ;;
  esac
done

# Check if both arguments were provided
if [ -z "$RAW_DIR" ] || [ -z "$OUTPUT_DIR" ]; then
  usage
fi

# Define directories for outputs
PROCESSED_DIR="${OUTPUT_DIR}/processed"
RESULTS_DIR="${OUTPUT_DIR}/results"
LOG_FILE="${RESULTS_DIR}/pipeline.log"

# Create necessary directories if they do not exist
mkdir -p "$PROCESSED_DIR"
mkdir -p "$RESULTS_DIR"

# Load configurations(replace path of softwares with local path)
TRM=/path_to/Trimmomatic-0.36/trimmomatic-0.36.jar
ADP=/path_to/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa
STAR_MAPPER=/path_to/STAR-2.5.3a/bin/Linux_x86_64/STAR
F_COUNTS=/path_to/subread-1.5.3-source/bin/featureCounts
STAR_REF_HOST=/path_to/Transcriptional_regulation_oomycetes/data/processed/STAR/STAR_HOST/star_index
GENOME_HOST=/path_to/Transcriptional_regulation_oomycetes/data/support/Pl_halstedii_study/Han.fa
GTF_HOST=/path_to/Transcriptional_regulation_oomycetes/data/support/Pl_halstedii_study/Han.gtf
STAR_REF_PATHOGEN=/path_to/Transcriptional_regulation_oomycetes/data/processed/STAR/star_index
GENOME_PATHOGEN=/path_to/Transcriptional_regulation_oomycetes/data/support/Pl_halstedii_study/Phals.fa
GTF_PATHOGEN=/path_to/Transcriptional_regulation_oomycetes/data/support/Pl_halstedii_study/Phals.gtf

# Log start time
echo "Pipeline started at $(date)" > "$LOG_FILE"

# Function to generate genome index
generate_genome_index() {
  local STAR_REF=$1
  local GENOME_FA=$2
  local GENOME_GTF=$3

  if [ ! -d "$STAR_REF" ]; then
    mkdir -p "$STAR_REF"
   $STAR_MAPPER  --runThreadN 20 --runMode genomeGenerate --genomeDir "$STAR_REF" --genomeFastaFiles "$GENOME_FA" --sjdbGTFfile "$GENOME_GTF" --sjdbOverhang 99
  fi
}

# Generate genome indices if not already present
generate_genome_index "$STAR_REF_HOST" "$GENOME_HOST" "$GTF_HOST"
generate_genome_index "$STAR_REF_PATHOGEN" "$GENOME_PATHOGEN" "$GTF_PATHOGEN"

# Function to process a sample
process_sample() {
  local SAMPLE=$1

  FNAME=$(basename "$SAMPLE" "_1.fq.gz")
  echo "Processing sample: $FNAME" >> "$LOG_FILE"

  # Trim reads
  java -jar $TRM PE -threads 15 -phred33 "$RAW_DIR/${FNAME}_1.fq.gz" "$RAW_DIR/${FNAME}_2.fq.gz" "$PROCESSED_DIR/${FNAME}_paired_1.fq.gz" "$PROCESSED_DIR/${FNAME}_unpaired_1.fq.gz" "$PROCESSED_DIR/${FNAME}_paired_2.fq.gz" "$PROCESSED_DIR/${FNAME}_unpaired_2.fq.gz" ILLUMINACLIP:$ADP:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:5:20 MINLEN:30 2>> "$LOG_FILE"
   
  # Unzip paired-end files

   gunzip "$PROCESSED_DIR/${FNAME}_paired_1.fq.gz" "$PROCESSED_DIR/${FNAME}_paired_2.fq.gz"

   # Remove unpaired files
   rm -f "$PROCESSED_DIR/${FNAME}_unpaired_1.fq.gz" "$PROCESSED_DIR/${FNAME}_unpaired_2.fq.gz"


  # Map to host
   $STAR_MAPPER  --genomeDir $STAR_REF_HOST --outSAMstrandField intronMotif --readFilesIn "$PROCESSED_DIR/${FNAME}_paired_1.fq.gz" "$PROCESSED_DIR/${FNAME}_paired_2.fq.gz" --limitBAMsortRAM 20000000000 --outSAMtype BAM SortedByCoordinate --runThreadN 20 --outFileNamePrefix "$RESULTS_DIR/${FNAME}_host_" 2>> "$LOG_FILE"

  # Map to pathogen
   $STAR_MAPPER  --genomeDir $STAR_REF_PATHOGEN --outSAMstrandField intronMotif --readFilesIn "$PROCESSED_DIR/${FNAME}_paired_1.fq.gz" "$PROCESSED_DIR/${FNAME}_paired_2.fq.gz" --limitBAMsortRAM 20000000000 --outSAMtype BAM SortedByCoordinate --runThreadN 20 --outFileNamePrefix "$RESULTS_DIR/${FNAME}_pathogen_" 2>> "$LOG_FILE"
}

# Process all samples
for SAMPLE in "$RAW_DIR"/*_1.fq.gz; do
  process_sample "$SAMPLE"
done

# Perform feature counting
$F_COUNTS -p -t exon -g gene_id -a $GTF_HOST -o "$RESULTS_DIR/final_counts_host.txt" -T 24 "$RESULTS_DIR"/*host_Aligned.sortedByCoord.out.bam 2>> "$LOG_FILE"
$F_COUNTS -p -t exon -g gene_id -a $GTF_PATHOGEN -o "$RESULTS_DIR/final_counts_pathogen.txt" -T 24 "$RESULTS_DIR"/*pathogen_Aligned.sortedByCoord.out.bam 2>> "$LOG_FILE"

# Log end time
echo "Pipeline completed at $(date)" >> "$LOG_FILE"

