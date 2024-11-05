
# Changelog

All notable changes to this project will be documented in this file. This project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]
- Planned additions to the project, including new data analysis methods and support for additional species.

## [1.3.0] - 2024-10-29
### Added
- Introduced `INSTALL.md` to guide users through setting up the software and dependencies.
- Added effector protein motif information from Sharma et al. (2015) for the `orthologs_motifs_study`.
- Included a link to the JASPAR 2020 database in the `orthologs_TF_study` directory.
- Created `LICENSE` file to provide project licensing details.
- New sections for README files in `/data/processed/` and `/data/support/` folders.

### Changed
- Modified directory structure documentation to include links to pipeline diagrams in `figures/`.
- Updated `README.md` files across the project to reflect additional datasets and clearer descriptions.

### Fixed
- Corrected broken links in `/data/support/README.md`.
- Resolved incorrect reference in the `orthologs_TF_study` documentation pointing to outdated resources.

## [1.2.0] - 2024-10-20
### Added
- Introduced `CODE_OF_CONDUCT.md` and `CONTRIBUTING.md` to facilitate community participation.
- Added links to effector sequences and blast results in `orthologs_motifs_study`.
- `00_Overview.md` in the `/docs` folder was expanded to provide comprehensive context on transcriptional regulation in oomycetes.

### Changed
- Refined the `04_Results.md` to provide a better summary of the research results, including adding images and tables.
- Improved formatting in README files to ensure species names are italicized properly.

### Fixed
- Addressed incorrect directory path issues in `/data/processed/README.md`.
- Fixed typos and broken markdown formatting in various README files.

## [1.1.0] - 2024-10-15
### Added
- Initial creation of `00_Overview.md`, `01_Methodology.md`, and other documentation files in `/docs`.
- `scripts/shell/Pl_halstedii_study/trimming_mapping_host_pathogen.sh` script added for preprocessing and analysis.
- Added support for additional R Bioconductor packages required for some analyses.

### Changed
- Reorganized `scripts/` folder to separate shell, Python, and R scripts into different subdirectories for clarity.
- Updated `README.md` to include comprehensive funding and institutional acknowledgments.

### Fixed
- Fixed incorrect labels in figures used in `/figures/pipeline_1.png`.
- Resolved file path issues for support data in `/data/support/README.md`.

## [1.0.0] - 2024-10-01
### Added
- Initial release of the project repository for **Transcriptional Regulation in Oomycetes**.
- Added core scripts for data processing, analysis, and visualization.
- Documentation in `README.md` to describe repository structure and usage guidelines.

