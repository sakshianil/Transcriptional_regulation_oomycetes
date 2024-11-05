# Environment Setup

This document provides detailed information about the environment required for running the scripts and analyses in the **Transcriptional Regulation in Oomycetes** project.

## Software and Tools

Below is a list of software and tools that are essential for this project, along with the recommended versions.

### Operating System
- This project has been developed and tested on the following operating systems:
  - macOS (Catalina or newer)
  - Ubuntu Linux 18.04 or newer

### Programming Languages
- **Python**: Version 3.7 or newer
- **R**: Version 4.0 or newer
- **Perl**: Version 5.30 or newer
- **Bash**: Default shell for Linux and macOS systems

### Required Software
- **Trimmomatic**: Version 0.36 - Used for trimming sequencing reads.
- **STAR**: Version 2.5.3a - Used for aligning RNA-Seq reads to the genome.
- **FeatureCounts (Subread package)**: Version 2.0.1 - Used for quantifying aligned reads.
- **MEME Suite**: Version 5.2.0 - Used for motif analysis.
- **BLAST**: Version 2.12.0 - Used for sequence alignment.
- **Clustal Omega**: Version 1.2.4 - Used for multiple sequence alignment.

### Libraries and Dependencies

#### Python
- **pandas**: Version 1.1.5 or newer - For data manipulation.
- **numpy**: Version 1.19.5 or newer - For numerical operations.
- **Biopython**: Version 1.78 - For parsing biological data.
- **requests**: Version 2.24.0 - For API requests.

To install the Python libraries, you can run:

```bash
pip install pandas==1.1.5 numpy==1.19.5 biopython==1.78 requests==2.24.0
```

#### R
- **Bioconductor Packages**: Various Bioconductor packages are required, including:
  - **edgeR**: For differential gene expression analysis.
  - **ggplot2**: For generating visualizations.

You can install Bioconductor packages with:

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("edgeR", "ggplot2"))
```

### Setting Up the Environment

1. **Clone the Repository**:
   ```bash
   git clone https://github.com/your-username/Transcriptional-regulation-in-oomycetes.git
   cd Transcriptional-regulation-in-oomycetes
   ```

2. **Install Required Software**:
   Ensure you have installed all the necessary software and programming languages according to the versions specified above.

3. **Set Environment Variables**:
   Set the environment variables for BLAST and other tools in your shell profile file (`.bashrc`, `.bash_profile`, etc.):

   ```bash
   export PATH="$PATH:/path/to/ncbi-blast/bin"
   export PATH="$PATH:/path/to/star/bin"
   export PATH="$PATH:/path/to/meme/bin"
   ```

   Modify `/path/to/` according to your actual installation paths.

4. **Install Python and R Libraries**:
   Make sure the Python and R libraries are installed as described above.

### Hardware Requirements

- **Memory**: At least 16 GB of RAM is recommended, especially for processing large RNA-Seq datasets.
- **Storage**: Minimum of 100 GB of free disk space.
- **CPU**: Multi-core processor (quad-core or higher is recommended) to efficiently run STAR, BLAST, and other intensive tasks.

### Additional Notes

- It is recommended to use a virtual environment for Python to prevent conflicts with system libraries. You can create and activate a virtual environment as follows:

  ```bash
  python3 -m venv env
  source env/bin/activate
  ```

  After activating the environment, install the necessary Python packages with `pip`.

## Docker Option (Optional)

For convenience, we provide a Dockerfile that allows you to create a containerized environment with all the tools pre-installed. This is especially useful for maintaining a consistent environment across different systems.

To build the Docker image:

```bash
docker build -t oomycetes_env .
```

To run the container:

```bash
docker run -it -v $(pwd):/workspace oomycetes_env
```

## Troubleshooting

If you encounter any issues while setting up the environment:

- **Check Dependencies**: Verify that all required software and dependencies are correctly installed.
- **Check Paths**: Ensure all necessary paths are correctly added to your `PATH` variable.
- **Memory**: If you experience crashes, consider increasing the memory allocation.

