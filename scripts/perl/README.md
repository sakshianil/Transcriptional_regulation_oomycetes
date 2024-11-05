# Perl Script for Markov Model Estimation

## Overview

This script estimates a Markov model from a FASTA file of DNA or protein sequences. The Markov model can be used for background correction in motif analysis. The script processes sequences by removing ambiguous characters and combining reverse complements for DNA by default.

Please note that for your specific analysis, **you must set `$USE_RC = 0`** to avoid combining reverse complements, as the Python script used earlier has already reverse-complemented sequences on the negative strand.

## Requirements
- Perl (version 5.10 or higher)
- Input should be in FASTA format

## Usage
To run the script, use the following command:

```sh
perl fasta-get-markov.pl -m <order> <input.fasta >output.model
```

### Arguments
- **`-m <order>`**: Specify the order of the Markov model to use. Default is `0`.
- **`-p`**: Use the protein alphabet instead of DNA. By default, the script uses the DNA alphabet.

### Example Usage

```sh
perl fasta-get-markov.pl -m 3 <bg_promoter.fasta >bg_promoter.model
```

This command will estimate a 3rd-order Markov model for the sequences in `bg_promoter.fasta` and save the output to `bg_promoter.model`.

## Important Note
The script has a default parameter **`$USE_RC = 1`**, which means that the Markov model combines both strands of DNA by default.

For your analysis, **set `$USE_RC = 0`** before running the script, because the Python script previously used has already taken care of reverse complements for negative strands. This adjustment prevents redundant handling of reverse complement sequences.

You can do this by editing the following line in the script:

```perl
$USE_RC = 0;
```

## Script Details
- **Author**: Timothy L. Bailey
- **Creation Date**: July 1, 2002
- **Version History**: The script has undergone various updates, initially created in 2002.

## Input and Output
- **Input**: A FASTA file of DNA or protein sequences.
- **Output**: The estimated Markov model, written to standard output.

## License
This script is provided under the copyright of Timothy L. Bailey, 2002. Please refer to the license terms for any usage or modification requirements.

## Contact
For any questions related to using this script, please reach out to [Sakshi Bharti](mailto:sakshi.bharti@senckenberg.de).


