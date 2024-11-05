import re
from collections import Counter
from collections import defaultdict
import pandas as pd
import string
import sys
import math
import csv
from pandas import read_csv
import pprint
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio import motifs

# Dictionaries to store genome, annotation, and sequence data
phals_dict = defaultdict(dict)
seq_dict = defaultdict(dict)
phal_gtf_dict = defaultdict(dict)


#Main function to execute the script
def main():
    """
    Main function to handle processing of genome and annotation data,
    and generate upstream promoter sequences.
    """
    # Define input,output file paths
    annotation_gtf = '/path_to/Transcriptional_regulation_oomycetes/data/support/Pl_halstedii_study/Phals.gtf'
    genome_fasta = '/path_to/Transcriptional_regulation_oomycetes/data/support/Pl_halstedii_study/Phals.fa'
    output_file = "/path_to/Transcriptional_regulation_oomycetes/results/python/Pl_halstedii_study/bg_promoter_100.fasta"

    # Load genome and GTF data
    print("Processing genome and annotation files...")
    genome_handler(genome_fasta)
    gtf_handler(annotation_gtf)

    # Generate upstream promoter sequences and save to output file
    print("Generating upstream promoter sequences...")
    with open(output_file, 'w') as output:
        for gene_id in phal_gtf_dict:
            upstream_generator(gene_id, output)
    print("Promoter sequences saved to", output_file)



#Function to generate upstream promoter sequences

def upstream_generator(gene_id, output):
      """
      Generates upstream promoter sequences for a given gene ID.
    
      Parameters:
      gene_id (str): The gene identifier.
      output (file object): File object to write the promoter sequences.
      """
        try:
                seq = seq_dict[phal_gtf_dict[gene_id]['scaffold_id']]

        except:
                print("Error with gene_id: %s" % (gene_id))


        strand = phal_gtf_dict[gene_id]['strand']

        promoter_len = 100

        if (strand == '+'):

                if (int(phal_gtf_dict[gene_id]['start']) >= promoter_len):
                        extract = seq [int(phal_gtf_dict[gene_id]['start'])-promoter_len:int(phal_gtf_dict[gene_id]['start'])-1]

                        seq_from = int(phal_gtf_dict[gene_id]['start'])-promoter_len-1
                        seq_to = int(phal_gtf_dict[gene_id]['start'])-1

                elif (int(phal_gtf_dict[gene_id]['start']) < promoter_len):
                        extract = seq[0:int(phal_gtf_dict[gene_id]['start'])-1]
                        seq_from = 0
                        seq_to = int(phal_gtf_dict[gene_id]['start'])-1
                output.write(">%s (%s:%i-%i strand:%s pattern:%s)\n%s" % (gene_id, phal_gtf_dict[gene_id]['scaffold_id'], seq_from, seq_to, strand, phals_dict[clusters], extract))

        elif (strand == '-'):


                if ((len(seq)-int(phal_gtf_dict[gene_id]['end'])) >= promoter_len):

                        extract = seq[int(phal_gtf_dict[gene_id]['end'])+1:int(phal_gtf_dict[gene_id]['end'])+promoter_len]
                        seq_from = int(phal_gtf_dict[gene_id]['end'])+1
                        seq_to = int(phal_gtf_dict[gene_id]['end'])+promoter_len

                elif ((len(seq)-int(phal_gtf_dict[gene_id]['end'])) < promoter_len):
                        extract = seq[int(phal_gtf_dict[gene_id]['end'])+1:len(seq)]
                        seq_from = int(phal_gtf_dict[gene_id]['end'])+1
                        seq_to = len(seq)

                extract = reverse_complement(extract)
                output.write(">%s (%s:%i-%i strand:%s)\n%s" % (gene_id, phal_gtf_dict[gene_id]['scaffold_id'], seq_from, seq_to, strand, phals_dict[clusters], extract))


        return None
#Function to get reverse complement of sequence on negative strand

def reverse_complement(seq):
    return Seq(seq).reverse_complement()

#Function to handle genome FASTA file and populate sequence dictionary

def genome_handler(genome_fasta):

    fil = genome_fasta
    cur_scaf = ''
    cur_seq = []
    for line in open(fil):
        if line.startswith(">") and cur_scaf == '':
            cur_scaf = line.split(' ')[0]

        elif line.startswith(">") and cur_scaf != '':
            seq_dict[cur_scaf.replace(">", "")] = ''.join(cur_seq)
            cur_scaf = line.split(' ')[0]
            cur_seq = []
        else:
            cur_seq.append(line.rstrip())
    seq_dict[cur_scaf.replace(">", "")] = ''.join(cur_seq)

    return None

#Function to handle GTF annotation file and populate annotation dictionary

def gtf_handler(phal_gtf):

        file = open(phal_gtf, "r")
        for line in file:
            la = line.strip().split("\t")
            try:
                if la[2] == "gene":
                    gene_id = la[8]
                    phal_gtf_dict[gene_id]['scaffold_id'] = la[0]
                    phal_gtf_dict[gene_id]['start'] = la[3]
                    phal_gtf_dict[gene_id]['end'] = la[4]
                    phal_gtf_dict[gene_id]['strand'] = la[6]
            except:
                continue
        return None

#Entry point for the script
if __name__ == '__main__':
    main()
