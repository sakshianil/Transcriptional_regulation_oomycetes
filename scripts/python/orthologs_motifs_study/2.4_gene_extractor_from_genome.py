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

#Dictionaries to store genome, annotation, and sequence data
# global variables
pattern = "N"
p_dict = defaultdict(dict) 
clusters = defaultdict(dict) 
seq_dict = defaultdict(dict) 
p_gtf_dict = defaultdict(dict) # dict storing gtf information

#Replace the path of gtf and fasta files of all genomes (Ph. sojae, Ph. infestans, Pl. halstedii) in order to gene sequences
def main():
       #Replace the path of gtf and fasta files of specific genome in order to retreive gene sequences
        p_gtf = '/path_to/Transcriptional_regulation_oomycetes/data/support/orthologs_motifs_study/Psoj.gtf'
        genome_fasta = '/path_to/Transcriptional_regulation_oomycetes/data/support/orthologs_motifs_study/Psoj.fa'

        genome_handler(genome_fasta)
        gtf_handler(p_gtf)

        pp = pprint.PrettyPrinter(indent=4)
        #replace output file name (optional)
        output = open(("Psoj_gene.fasta"),'w')

        for genes in p_gtf_dict:

                gene_generator(genes, output)
        output.close()


def gene_generator(gene_id, output):


        try:
                seq = seq_dict[p_gtf_dict[gene_id]['scaffold_id']]

        except:
                print("Error with gene_id: %s" % (gene_id))


        strand = p_gtf_dict[gene_id]['strand']


        if (strand == '+'):
                extract = seq[int(p_gtf_dict[gene_id]['start'])-1:int(p_gtf_dict[gene_id]['end'])]
                seq_from = int(p_gtf_dict[gene_id]['start'])
                seq_to = int(p_gtf_dict[gene_id]['end'])
                output.write("\n>%s (%s:%i-%i strand:%s)\n%s" % (gene_id, p_gtf_dict[gene_id]['scaffold_id'], seq_from, seq_to, strand, extract))



        elif (strand == '-'):
                extract = seq[int(p_gtf_dict[gene_id]['start'])-1:int(p_gtf_dict[gene_id]['end'])]
                seq_from = int(p_gtf_dict[gene_id]['start'])
                seq_to = int(p_gtf_dict[gene_id]['end'])
                extract = reverse_complement(extract)
                output.write("\n>%s (%s:%i-%i strand:%s)\n%s" % (gene_id, p_gtf_dict[gene_id]['scaffold_id'], seq_from, seq_to, strand, extract))

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
                        seq_dict[cur_scaf.replace(">","")] = ''.join(cur_seq)
                        cur_scaf = line.split(' ')[0]
                        cur_seq = []
                else:
                        cur_seq.append(line.rstrip())
        seq_dict[cur_scaf.replace(">","")] = ''.join(cur_seq)

        return None

#Function to handle GTF annotation file and populate annotation dictionary

def gtf_handler(p_gtf):
        file = open(p_gtf, "r")
        for line in file:
                la = line.strip().split("\t")
                try:
                        if la[2] == "gene":
                                gene_id = la[8]
                                p_gtf_dict[gene_id]['scaffold_id'] = la[0]
                                p_gtf_dict[gene_id]['start'] = la[3]
                                p_gtf_dict[gene_id]['end'] = la[4]
                                p_gtf_dict[gene_id]['strand'] = la[6]
                except:
                        continue
        return None

#Entry point for the script
if __name__ == '__main__':
    main()
