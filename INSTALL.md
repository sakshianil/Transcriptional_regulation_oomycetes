
# Installation Guide

## Overview

Welcome to the installation guide for the **Transcriptional Regulation in Oomycetes** project. This guide will help you set up the environment to reproduce the results from our research. Follow the instructions below to install the required software, libraries, and configure the environment correctly.

## System Requirements

- **Operating System**: macOS, Linux, or Windows (using WSL).
- **RAM**: Minimum 8 GB, recommended 16 GB or higher for larger datasets.
- **Storage**: At least 20 GB of free disk space.

## Prerequisites

Before starting, ensure that you have the following installed on your machine:

### Required Software

1. **Git** - for cloning the repository
   ```sh
   # Install Git on Linux/macOS
   sudo apt-get install git
   
   # For macOS using Homebrew
   brew install git
   ```

2. **Python** - version 3.6 or above
   - Install using [Anaconda](https://www.anaconda.com/products/distribution) or:
   ```sh
   sudo apt-get install python3
   ```

3. **R and Perl** - required for certain scripts
   ```sh
   # Install R
   sudo apt-get install r-base

   # Install Perl
   sudo apt-get install perl
   ```

4. **Other Required Software**
   - **Trimmomatic** (version 0.36)
   - **STAR** (version 2.5.3a)
   - **FeatureCounts** (Subread package version 2.0.1)
   - **MEME Suite** (version 5.2.0)
   - **BLAST** (version 2.12.0)
   - **Clustal Omega** (version 1.2.4)

### Python Libraries

Install the following Python packages:

- `pandas`
- `numpy`

You can install these libraries using:

```sh
pip install pandas numpy
```

## Step-by-Step Installation

### Step 1: Clone the Repository

First, clone this repository to your local machine using the command below:

```sh
git clone https://github.com/your-username/Transcriptiona_regulation_oomycetes.git
```

Navigate to the cloned repository:

```sh
cd Transcriptiona_regulation_oomycetes
```

### Step 2: Install Dependencies

#### R Libraries

Some scripts require R packages, primarily from **Bioconductor**. Run the following commands in R to install the required packages:

```R
# Open R
R

# Install Bioconductor dependencies
install.packages("BiocManager")
BiocManager::install(c("DESeq2", "edgeR"))
```

#### Software Installation

- [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
- [STAR](https://github.com/alexdobin/STAR)
- [FeatureCounts](http://subread.sourceforge.net/)
- [MEME Suite](http://meme-suite.org)
- [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
- [Clustal Omega](http://www.clustal.org/omega/)

After downloading and installing these tools, ensure to add each to your system’s PATH variable.

### Step 3: Set Up Environment Variables

Add the installed software to your PATH environment variable. Open your `.bashrc` or `.zshrc` file and add the following lines:

```sh
# Example for STAR and BLAST
export PATH=$PATH:/path/to/STAR
export PATH=$PATH:/path/to/blast/bin
```

Apply changes:

```sh
source ~/.bashrc
```

### Step 4: Directory Structure

Ensure the directory structure is correctly set up:

- **data/**: Contains raw and processed datasets.
  - **raw/**: Raw data from experiments or databases.
  - **processed/**: Processed datasets ready for analysis.
  - **support/**: Supporting files like FASTA, GFF, and GTF.
- **scripts/**: All the computational scripts, organized by language.
- **results/**: Outputs generated from running analyses.
- **docs/**: Project documentation.
- **figures/**: Figures generated during analyses.

All downloaded data files should be stored in `data/raw`, and processed outputs should be in `data/processed`.

### Step 5: Running Scripts

#### Trimming and Mapping

To run the trimming and mapping process, execute the following script:

```sh
bash scripts/shell/Pl_halstedii_study/trimming_mapping_host_pathogen.sh
```

#### Running Python Scripts

Navigate to the `scripts/python` folder and run the Python scripts for data processing:

```sh
python script_name.py
```

#### R Scripts for Statistical Analysis

Execute R scripts:

```sh
Rscript script_name.R
```

Refer to the individual README files within the `scripts/` directories for detailed usage instructions.

## Troubleshooting

- **Missing Dependencies**: Ensure all dependencies are installed and PATH variables are set correctly.
- **Permission Errors**: Add `sudo` before commands if you encounter permission issues on Unix-based systems.
- **Environment Issues**: Use [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) to manage Python environments effectively.

## Contact

For installation issues or general questions, contact [Sakshi Bharti](mailto:sakshi.bharti@senckenberg.de).

