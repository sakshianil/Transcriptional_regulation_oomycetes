# Testing Guide

This document describes the testing procedures used to verify the correctness, performance, and reliability of the scripts and workflows included in the **Transcriptional Regulation in Oomycetes** project.

## Overview

The purpose of this document is to ensure that all scripts and analyses run as expected, producing accurate and reproducible results. Testing is a crucial step to validate data processing, analysis, and pipeline integrity, and to identify bugs or areas for optimization.

## Types of Testing

### 1. **Unit Testing**
Unit testing ensures that individual functions or sections of code work as intended. It is performed primarily on Python, R, Perl, and shell scripts in the repository.

- **Python Scripts**: The functions defined in the Python scripts are tested with mock datasets to verify their outputs.
- **R Scripts**: Unit tests are conducted using `testthat` for R, checking specific functions like those for data transformations and statistical modeling.

### 2. **Integration Testing**
Integration testing aims to ensure that all components in the pipeline interact with each other correctly.

- **Script Dependencies**: Tests are carried out to ensure data flows smoothly from one stage to another, such as from read trimming (using Trimmomatic) to mapping (using STAR) and quantification (using FeatureCounts).
- **Output File Integrity**: Tests are designed to check that the expected output formats are generated in each processing step and that file paths and parameters align correctly.

### 3. **End-to-End Testing**
End-to-end testing validates the entire pipeline, from raw input data to final output. This ensures that the project can be run smoothly on a local or remote machine, and that the generated results are consistent.

- This type of testing is particularly important for complex workflows like motif analysis, BLAST alignments, and ortholog identification.

### 4. **Performance Testing**
Performance testing is conducted to determine if the pipeline runs efficiently. This includes measuring:

- **Execution Time**: Measuring time taken by each stage of the pipeline and identifying bottlenecks.
- **Memory Usage**: Verifying that memory usage is optimized, especially when processing large datasets.

## Running Tests

### Python Scripts
Tests for Python functions are defined using `pytest`. To run the tests:

1. Navigate to the `scripts/python` directory.
2. Run the command:

   ```bash
   pytest
   ```

3. Any failed tests will be listed with detailed information.

### R Scripts
To run tests for R functions:

1. Make sure you have the `testthat` library installed:

   ```R
   install.packages("testthat")
   ```

2. Run tests by executing:

   ```R
   library(testthat)
   test_dir("scripts/R/tests")
   ```

### Shell Scripts
Shell scripts are tested manually in a safe environment. We recommend creating a test folder where intermediate outputs can be generated for comparison to expected results.

- Run the trimming, alignment, and counting shell scripts with a subset of data to verify the correct processing of reads.
- Test script paths and dependencies to ensure all required tools are available and accessible.

## Testing Datasets

For testing purposes, we have included sample datasets in the `data/support/testing` directory. These datasets are smaller and representative of the full datasets to ensure:

1. Quick and effective testing of functions.
2. Representative testing to identify possible issues.

## Testing Procedures

1. **Cloning the Repository**: Clone the repository and navigate to the project root.
   ```bash
   git clone https://github.com/your-username/Transcriptional-regulation-in-oomycetes.git
   cd Transcriptional-regulation-in-oomycetes
   ```

2. **Environment Setup**: Ensure the software dependencies are installed as per `ENVIRONMENT.md`. Activate the virtual environment, if used.

3. **Run Test Scripts**:
   - **Python**: Execute `pytest` to validate Python script functionality.
   - **R**: Use `testthat` to ensure R functions behave as expected.
   - **Shell**: Execute test runs of `shell` scripts using the small sample dataset provided.

4. **Validate Outputs**: Compare outputs with pre-defined expected results to verify consistency and accuracy.

## Troubleshooting

- **Environment Issues**: Ensure all dependencies are installed and set up as per `ENVIRONMENT.md`.
- **Path Issues**: Make sure all input/output paths are correct and the files exist in the expected locations.
- **Outdated Libraries**: Update software and libraries if issues persist with compatibility.

## Known Issues

- **Memory Limitations**: Some of the full dataset analyses require a significant amount of RAM. Testing should be performed on subsets of the data to prevent system crashes.
- **Execution Time**: End-to-end tests may take several hours. It is recommended to run them in a well-resourced environment or on high-performance computing clusters.

## Continuous Integration

To automate testing:
- Consider using GitHub Actions to implement Continuous Integration (CI) workflows.
- Define a `.github/workflows/test.yml` file to automatically run the tests each time a change is made to the repository.

## Contact

For further questions related to testing, please contact [Sakshi Bharti](mailto:sakshi.bharti@senckenberg.de).

