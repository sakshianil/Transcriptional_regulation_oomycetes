import os
import time
from Bio import Entrez, SeqIO
from Bio.Seq import Seq

# Set the email address for NCBI Entrez
Entrez.email = "sakshi.bharti@senckenberg.de"

/Users/sbharti/Documents/GitHub/Transcriptional_regulation_oomycetes/results/python/orthologs_motifs_study/clusters/gene 
# Directory paths
input_dir = "/path_to/Transcriptional_regulation_oomycetes/results/python/orthologs_motifs_study/clusters/gene"
output_dir = "/path_to/Transcriptional_regulation_oomycetes/results/python/orthologs_motifs_study/clusters/promoter"

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

def extract_accession_from_fasta(fasta_file_path):
    """
    Extract the accession number from the FASTA file header.
    Assumes that the header line starts with '>' and the accession follows directly.
    """
    with open(fasta_file_path, "r") as fasta_file:
        for line in fasta_file:
            if line.startswith(">"):
                # Extract the accession number (text after '>' until the first space or dot)
                accession = line.split()[0][1:]  # Skip the '>'
                return accession
    return None

def fetch_upstream_region(accession, upstream_length=400):
    try:
        # Use efetch to get the nucleotide record and location
        handle = Entrez.efetch(db="nuccore", id=accession, rettype="gb", retmode="text")
        nucleotide_record = SeqIO.read(handle, "genbank")
        handle.close()

        # Check the strand and fetch upstream sequence
        for feature in nucleotide_record.features:
            if feature.type == "gene" and 'locus_tag' in feature.qualifiers:
                # Assuming the gene is the first feature with a locus tag
                start = int(feature.location.start)
                strand = feature.location.strand

                # Calculate upstream region start and stop positions
                if strand == 1:  # Positive strand
                    upstream_start = max(1, start - upstream_length)
                    upstream_end = start - 1
                elif strand == -1:  # Negative strand
                    upstream_start = start + 1
                    upstream_end = start + upstream_length

                # Fetch the upstream region
                handle = Entrez.efetch(db="nuccore", id=accession, seq_start=upstream_start, seq_stop=upstream_end, rettype="fasta", retmode="text")
                upstream_sequence = handle.read()
                handle.close()

                # If on the negative strand, take the reverse complement
                if strand == -1:
                    seq_record = SeqIO.read(Entrez.efetch(db="nuccore", id=accession, seq_start=upstream_start, seq_stop=upstream_end, rettype="fasta", retmode="text"), "fasta")
                    seq_record.seq = seq_record.seq.reverse_complement()
                    upstream_sequence = seq_record.format("fasta")

                return upstream_sequence

        print(f"No gene feature found in nucleotide record for {accession}.")
        return None

    except Exception as e:
        print(f"Error fetching upstream region for nucleotide {accession}: {e}")
        return None

# Loop over all fasta files in the input directory
for filename in os.listdir(input_dir):
    if filename.endswith("_gene.fasta"):
        fasta_file_path = os.path.join(input_dir, filename)

        # Extract the accession number from the FASTA header
        accession = extract_accession_from_fasta(fasta_file_path)
        if not accession:
            print(f"Could not extract a valid accession from {filename}")
            continue

        print(f"Fetching upstream region for accession {accession}")

        # Fetch upstream sequence
        upstream_sequence = fetch_upstream_region(accession)

        # Save the upstream sequence to a file
        if upstream_sequence:
            # Constructing the filename with protein_id and gene_id
            gene_id = accession  # Adjust if needed to extract protein_id
            output_file = os.path.join(output_dir, f"{gene_id}_upstream.fasta")
            try:
                with open(output_file, "w") as f:
                    f.write(upstream_sequence)
                print(f"Saved upstream region for {accession} to {output_file}")
            except Exception as e:
                print(f"Error writing file {output_file}: {e}")

        # Introduce a delay to avoid hitting NCBI rate limits
        time.sleep(1)

