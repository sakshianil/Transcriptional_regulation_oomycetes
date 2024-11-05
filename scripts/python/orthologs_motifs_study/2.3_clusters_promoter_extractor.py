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

pattern = "N"
phals_dict = defaultdict(dict)
clusters = defaultdict(dict)
seq_dict = defaultdict(dict)
phal_gtf_dict = defaultdict(dict)

#Main function to execute the script

def main():

    clusters_degreport = '/path_to/Transcriptional_regulation_oomycetes/results/R/gene_clusters_minc15_reduceT.txt'
    phal_gtf = '/path_to/Transcriptional_regulation_oomycetes/data/support/Phals.gtf'
    genome_fasta = '/path_to/Transcriptional_regulation_oomycetes/data/support/Phals.fa'
    
    read_cluster(clusters_degreport)
    genome_handler(genome_fasta)
    gtf_handler(phal_gtf)


    pp = pprint.PrettyPrinter(indent=4)

    for query_pattern in clusters:
        if len(clusters[query_pattern].split(",")) >= 1:

            ptrn = query_pattern.replace(',', '')
            output = open(("%s.fasta") % (ptrn), 'w')

            for genes in clusters[query_pattern].split(","):
                upstream_generator(genes, query_pattern, output)

            output.close()



#Function to read cluster file geenerated from R script

def read_cluster(clust_file):

    pp = pprint.PrettyPrinter(indent=1)

    clust = open(clust_file, "r")

    for line in clust:
        if line.startswith("genes"):
            continue
        line = line.strip().split('\t')
        id = line[1]
        value = line[0]
        try:
            clusters[id] = clusters[id] + "," + value

        except:
            clusters[id] = value


#Function to generate upstream promoter sequences (400 nucleotides upstream)

def upstream_generator(gene_id, cluster_id, output):

        query_pattern = cluster_id

        try:
                seq = seq_dict[phal_gtf_dict[gene_id]['scaffold_id']]

        except:
                print("Error with gene_id: %s" % (gene_id))


        ptrn = query_pattern.replace(',', '')
        strand = phal_gtf_dict[gene_id]['strand']

        promoter_len = 400

        if (strand == '+'):

                if (int(phal_gtf_dict[gene_id]['start']) >= promoter_len):
                        extract = seq [int(phal_gtf_dict[gene_id]['start'])-promoter_len-1:int(phal_gtf_dict[gene_id]['start'])-1]
                        seq_from = int(phal_gtf_dict[gene_id]['start'])-promoter_len-1
                        seq_to = int(phal_gtf_dict[gene_id]['start'])-1
                elif (int(phal_gtf_dict[gene_id]['start']) < promoter_len):
                        extract = seq[0:int(phal_gtf_dict[gene_id]['start'])-1]
                        seq_from = 0
                        seq_to = int(phal_gtf_dict[gene_id]['start'])-1
                output.write(">%s(%s:%i-%i:strand:%s)\n%s\n" % (gene_id, phal_gtf_dict[gene_id]['scaffold_id'], seq_from, seq_to, strand, extract))

        elif (strand == '-'):

                if ((len(seq)-int(phal_gtf_dict[gene_id]['end'])) >= promoter_len):

                        extract = seq[int(phal_gtf_dict[gene_id]['end'])+1:int(phal_gtf_dict[gene_id]['end'])+promoter_len+1]
                        seq_from = int(phal_gtf_dict[gene_id]['end'])+1
                        seq_to = int(phal_gtf_dict[gene_id]['end'])+promoter_len+1
                elif ((len(seq)-int(phal_gtf_dict[gene_id]['end'])) < promoter_len):
                        extract = seq[int(phal_gtf_dict[gene_id]['end'])+1:len(seq)]
                        seq_from = int(phal_gtf_dict[gene_id]['end'])+1
                        seq_to = len(seq)

                extract = reverse_complement(extract)
                output.write(">%s(%s:%i-%i:strand:%s)\n%s\n" % (gene_id, phal_gtf_dict[gene_id]['scaffold_id'], seq_from, seq_to, strand, extract))


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
                    gene_id = re.findall(r'(?:CEG|CEJ|CEI|ENSRNAG)\d+', la[8])[0]
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
