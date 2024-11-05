import csv
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
# global variables
pattern = "N"
phals_dict = defaultdict(dict) 
clusters = defaultdict(dict) 
seq_dict = defaultdict(dict) 
phal_gtf_dict = defaultdict(dict) # dict storing gff information


#Main function to execute the script
def main():
        phals = read_csv('/path_to/Transcriptional_regulation_oomycetes/results/R/gene_clusters_minc15_reduceT.txt', sep='\t')
        phal_gtf = '/path_to/Transcriptional_regulation_oomycetes/data/support/Phals.gtf'
        genome_fasta = '/path_to/Transcriptional_regulation_oomycetes/data/support/Phals.fa'
        pattern_generator(phals)
        genome_handler(genome_fasta)
        gtf_handler(phal_gtf)

        pp = pprint.PrettyPrinter(indent=4)

        for query_pattern in clusters:
                if len(clusters[query_pattern]['genes']) >= 1:

                        ptrn = query_pattern.replace(',', '')
                        output = open(("%s.fasta") % (ptrn), 'w')

                        for genes in clusters[query_pattern]['gene_ids']:
                                gene_generator(genes, query_pattern, output)

                        output.close()

#Function to generate gene sequences
def gene_generator(gene_id, query_pattern, output):

        try:
                seq = seq_dict[phal_gtf_dict[gene_id]['scaffold_id']]

        except:
                print("Error with gene_id: %s" % (gene_id))


        ptrn = query_pattern.replace(',', '')
        strand = phal_gtf_dict[gene_id]['strand']


        if (strand == '+'):
                extract = seq[int(phal_gtf_dict[gene_id]['start'])-1:int(phal_gtf_dict[gene_id]['end'])]
                seq_from = int(phal_gtf_dict[gene_id]['start'])
                seq_to = int(phal_gtf_dict[gene_id]['end'])
                output.write(">%s (%s:%i-%i:strand:%s) \n%s\n" % (gene_id, phal_gtf_dict[gene_id]['scaffold_id'], seq_from, seq_to, strand, extract))

        elif (strand == '-'):
                extract = seq[int(phal_gtf_dict[gene_id]['start'])-1:int(phal_gtf_dict[gene_id]['end'])]
                seq_from = int(phal_gtf_dict[gene_id]['start'])
                seq_to = int(phal_gtf_dict[gene_id]['end'])
                extract = reverse_complement(extract)
                output.write(">%s (%s:%i-%i:strand:%s) \n%s\n" % (gene_id, phal_gtf_dict[gene_id]['scaffold_id'], seq_from, seq_to, strand, extract))

        return None
#Function to get reverse complement of sequence on negative strand

def reverse_complement(seq):
    return Seq(seq).reverse_complement()

#Function to get gene ids
def pattern_generator(phals):

        for i in range(len(phals['genes'])):
                phals_dict[i]['pattern'] = ",".join(phals.iloc[i, 1:].apply(str).values)
                phals_dict[i]['N_count'] = phals_dict[i]['pattern'].count(pattern)
                phals_dict[i]['genes'] = phals['genes'][i]
        
        for index in phals_dict.keys():
                if phals_dict[index]['N_count'] == 0:
                        try:
                                clusters[phals_dict[index]['pattern']]['genes'] += "," + phals_dict[index]['genes']

                        except:

                                clusters[phals_dict[index]['pattern']]['genes'] = phals_dict[index]['genes']
        key_list = clusters.keys()

        for ch in key_list:
                p = re.findall(r'(?:CEG|CEJ|CEI|ENSRNAG)\d+', clusters[ch]['genes'])

                clusters[ch]['gene_ids'] = p

        return None

#Function to handle genome FASTA file and populate sequence dictionary

def genome_handler(genome_fasta):

        fil = genome_fasta
        cur_scaf = ''
        cur_seq = []
        for line in open(fil):
                if line.startswith(">") and cur_scaf == '':
                        cur_scaf = line.split(' ')[0]

                elif line.startswith(">") and cur_scaf != '':
                        seq_dict[cur_scaf.replace(">","")] = ''.join(cur_seq)
                        cur_scaf = line.split(' ')[0]
                        cur_seq = []
                else:
                        cur_seq.append(line.rstrip())
        seq_dict[cur_scaf.replace(">","")] = ''.join(cur_seq)

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
