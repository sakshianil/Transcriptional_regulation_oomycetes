#!/bin/bash

# Source the configuration file
if [ -f ./config.sh ]; then
  source ./config.sh
else
  echo "Configuration file config.sh not found!"
  exit 1
fi

# Function to display usage
usage() {
  echo "
Usage: $0 [host|pathogen] *_1.fq.gz
Specify 'host' or 'pathogen' as the first argument.
"
  exit 1
}

# Check for the correct number of arguments
if [ $# -lt 2 ]; then
  usage
fi

# Set variables based on the type of analysis
if [ "$1" == "host" ]; then
  MAPPED=$MAPPED_HOST
  STAR_REF=$STAR_REF_HOST
  GENOME=$GENOME_HOST
  GTF=$GTF_HOST
  shift
elif [ "$1" == "pathogen" ]; then
  MAPPED=$MAPPED_PATHOGEN
  STAR_REF=$STAR_REF_PATHOGEN
  GENOME=$GENOME_PATHOGEN
  GTF=$GTF_PATHOGEN
  shift
else
  usage
fi

# Create directories if they do not exist
[ ! -d "$STAR_REF" ] && mkdir -p $STAR_REF
[ ! -d "$MAPPED" ] && mkdir -p $MAPPED
[ ! -d "$DATA/trimmed" ] && mkdir -p $DATA/trimmed

# Generate genome index if it does not exist
if [ ! -f "$STAR_REF/Genome" ]; then
  STAR --runThreadN 20 --runMode genomeGenerate --genomeDir $STAR_REF --genomeFastaFiles $GENOME --sjdbGTFfile $GTF --sjdbOverhang 99
fi

# Process each input file
for i in $@
do 
  FNAME=$(sed 's/_1.fq.gz//' <<< $(basename $i))

  echo  "## Command used for file: $FNAME ##"

  # Trimmomatic command for trimming
  java -jar $TRM PE -threads 15 -phred33 $DATA/$FNAME\_1.fq.gz $DATA/$FNAME\_2.fq.gz \
  $DATA/trimmed/$FNAME\_paired_1.fq.gz $DATA/trimmed/$FNAME\_unpaired_1.fq.gz \
  $DATA/trimmed/$FNAME\_paired_2.fq.gz $DATA/trimmed/$FNAME\_unpaired_2.fq.gz \
  ILLUMINACLIP:$ADP:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:5:20 MINLEN:30  

  # Unzip paired trimmed files
  gunzip $DATA/trimmed/$FNAME\_paired_1.fq.gz $DATA/trimmed/$FNAME\_paired_2.fq.gz 

  # Remove unpaired trimmed files
  rm -rf $DATA/trimmed/$FNAME\_unpaired_*

  # STAR command for mapping
  STAR --genomeDir $STAR_REF --outSAMstrandField intronMotif \
  --readFilesIn $DATA/trimmed/$FNAME\_paired_1.fq $DATA/trimmed/$FNAME\_paired_2.fq \
  --limitBAMsortRAM 20000000000 --outSAMtype BAM SortedByCoordinate \
  --runThreadN 20 --outFileNamePrefix $MAPPED/$FNAME\_

  echo  "####"

done

