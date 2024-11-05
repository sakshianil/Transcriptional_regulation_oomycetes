import requests
import pandas as pd

# Input and output file paths
input_file_path = "/path_to/Transcriptional_regulation_oomycetes/results/shell/orthologs_TF_study/combined_with_profile_summaries.tsv"  # #Adjust path to your input file for the output .fasta file
fasta_output_file = "/path_to/Transcriptional_regulation_oomycetes/results/python/orthologs_TF_study/output_query_sequences.fasta"  # Path for the output .fasta file

# Load the input data
data = pd.read_csv(input_file_path, sep='\t')

# Function to fetch FASTA sequence from UniProt using UniProt ID
def fetch_fasta_sequence(uniprot_id):
    try:
        url = f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta"
        response = requests.get(url)

        # Check if the request was successful
        if response.status_code == 200:
            return response.text  # Return the FASTA sequence
        else:
            print(f"Error fetching sequence for UniProt ID {uniprot_id}, Status Code: {response.status_code}")
            return None
    except Exception as e:
        print(f"Error fetching sequence for UniProt ID {uniprot_id}: {e}")
        return None

# Open the output file for writing
with open(fasta_output_file, "w") as fasta_file:
    for index, row in data.iterrows():
        uniprot_id = row['Uniprot ID']

        # Skip if no Uniprot ID is available or it is 'N/A'
        if pd.isna(uniprot_id) or uniprot_id == 'N/A':
            continue

        # Fetch the FASTA sequence for the given Uniprot ID
        fasta_sequence = fetch_fasta_sequence(uniprot_id)

        # Write the FASTA sequence to the output file if available
        if fasta_sequence:
            fasta_file.write(fasta_sequence)

print(f"FASTA sequences saved to {fasta_output_file}.")

