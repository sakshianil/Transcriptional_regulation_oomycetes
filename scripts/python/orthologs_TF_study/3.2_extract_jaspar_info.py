import requests
import pandas as pd

# Input and output file paths
input_file_path = "/path_to/Transcriptional_regulation_oomycetes/results/shell/orthologs_TF_study/combined_tomtom_with_urls.tsv"
output_file_path = "/path_to/Transcriptional_regulation_oomycetes/results/python/orthologs_TF_study/combined_with_profile_summaries.tsv"

# Load the input data
data = pd.read_csv(input_file_path, sep='\t')

# JASPAR API URL
jaspar_base_url = "http://jaspar.genereg.net/api/v1/matrix/"

# Function to extract required information from the JASPAR API using matrix ID
def extract_info_from_api(matrix_id):
    if pd.isna(matrix_id):  # Check if the Matrix ID is NaN
        return {key: 'N/A' for key in ['Name', 'MatrixId', 'Class', 'Family', 'Collection', 'Taxon', 'Species', 'Validation ID', 'Uniprot ID', 'Source']}

    try:
        # Send a request to the JASPAR API with the matrix ID
        url = f"{jaspar_base_url}{matrix_id}/"
        response = requests.get(url, timeout=10)

        # Check if request was successful
        if response.status_code != 200:
            print(f"Error: Unable to fetch data from {url}, Status Code: {response.status_code}")
            return {key: 'Error' for key in ['Name', 'MatrixId', 'Class', 'Family', 'Collection', 'Taxon', 'Species', 'Validation ID', 'Uniprot ID', 'Source']}
        
        # Parse the JSON response
        data = response.json()

        # Extract relevant fields
        extracted_info = {
            'Name': data.get('name', 'N/A'),
            'MatrixId': data.get('matrix_id', 'N/A'),
            'Class': data.get('class', 'N/A'),
            'Family': data.get('family', 'N/A'),
            'Collection': data.get('collection', 'N/A'),
            'Taxon': data.get('taxon', 'N/A'),
            'Species': data.get('species', [{'name': 'N/A'}])[0].get('name', 'N/A'),  # List of species, extract first one
            'Validation ID': data.get('validation', 'N/A'),
            'Uniprot ID': data.get('uniprot_ids', ['N/A'])[0],  # List of Uniprot IDs, extract first one
            'Source': data.get('source', 'N/A'),
        }

        return extracted_info
    except Exception as e:
        print(f"Error fetching data for {matrix_id}: {e}")
        return {key: 'Error' for key in ['Name', 'MatrixId', 'Class', 'Family', 'Collection', 'Taxon', 'Species', 'Validation ID', 'Uniprot ID', 'Source']}

# Apply the extraction function to each row using the matrix ID
results = data['Target_ID'].apply(lambda matrix_id: extract_info_from_api(matrix_id))

# Convert results to a DataFrame and concatenate with the original data
results_df = pd.DataFrame(results.tolist())
output_data = pd.concat([data, results_df], axis=1)

# Save the output to a new TSV file
output_data.to_csv(output_file_path, sep='\t', index=False)
print(f"Output saved to {output_file_path}.")

