# config.sh
# Configuration file for paths and parameters

# Path to Trimmomatic
export TRM="/path/to/Trimmomatic-0.XX/trimmomatic-0.XX.jar"
export ADP="/path/to/Trimmomatic-0.XX/adapters/TruSeq3-PE-2.fa"

# Data directories
export DATA="/path/to/raw_data"
export MAPPED_HOST="/path/to/STAR/STAR_HOST"
export STAR_REF_HOST="/path/to/STAR/STAR_HOST/star_index"
export GENOME_HOST="/path/to/genome/Han.fa"
export GTF_HOST="/path/to/genome/Han.gtf"

export MAPPED_PATHOGEN="/path/to/STAR"
export STAR_REF_PATHOGEN="/path/to/STAR/star_index"
export GENOME_PATHOGEN="/path/to/genome/Phals.fa"
export GTF_PATHOGEN="/path/to/genome/Phals.gtf"

