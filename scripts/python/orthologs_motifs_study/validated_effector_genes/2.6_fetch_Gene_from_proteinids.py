import time
from Bio import Entrez
from urllib.error import HTTPError

Entrez.email = "write_here_your_email_id"

def fetch_gene_sequence(protein_id, retries=3):
    for _ in range(retries):
        try:
            # Search for the protein ID in the protein database
            search_handle = Entrez.esearch(db="protein", term=protein_id, retmax=1)
            search_result = Entrez.read(search_handle)
            search_handle.close()

            if search_result["IdList"]:
                # Get the protein ID
                protein_ncbi_id = search_result["IdList"][0]
                
                # Link the protein to the corresponding gene in the nucleotide database
                link_handle = Entrez.elink(dbfrom="protein", db="nucleotide", id=protein_ncbi_id, linkname="protein_nuccore")
                link_result = Entrez.read(link_handle)
                link_handle.close()

                # Check if there is a link to a gene
                if link_result and link_result[0]["LinkSetDb"]:
                    gene_id = link_result[0]["LinkSetDb"][0]["Link"][0]["Id"]

                    # Fetch the gene sequence in FASTA format
                    fetch_handle = Entrez.efetch(db="nucleotide", id=gene_id, rettype="fasta", retmode="text")
                    sequence = fetch_handle.read()
                    fetch_handle.close()
                    return sequence
            return None
        except HTTPError as e:
            if e.code == 429:
                print(f"Hit rate limit for {protein_id}. Retrying in 10 seconds...")
                time.sleep(10)
            else:
                raise
    print(f"Failed to retrieve sequence for {protein_id} after {retries} attempts.")
    return None

#create a new file of protein Ids named as "protein_ids_only.txt" by fetching second column from the file "path_to/Transcriptional_regulation_oomycetes/data/support/orthologs_motifs_study/validated_effectors/effectors_blast_results.txt" 
#Read protein IDs from file and fetch gene sequences
with open("/path_to/Transcriptional_regulation_oomycetes/data/support/orthologs_motifs_study/validated_effectors/protein_ids_only.txt", "r") as file:
    for line in file:
        protein_id = line.strip()  # Read each line as a protein ID
        if protein_id:  # Check if the line is not empty
            sequence = fetch_gene_sequence(protein_id)
            if sequence:
                with open(f"/path_to/Transcriptional_regulation_oomycetes/results/python/orthologs_motifs_study/clusters/gene/{protein_id}_gene.fasta", "w") as out_file:
                    out_file.write(sequence)
            else:
                print(f"Gene sequence not found for protein ID: {protein_id}")
            
            # Introduce a delay to avoid overwhelming NCBI servers
            time.sleep(1)
        else:
            print("Skipping empty line.")

