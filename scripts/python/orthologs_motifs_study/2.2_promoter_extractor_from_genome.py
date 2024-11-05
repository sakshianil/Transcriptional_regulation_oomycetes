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
p_dict = defaultdict(dict)
seq_dict = defaultdict(dict)
p_gtf_dict = defaultdict(dict)


#Main function to execute the script
def main():

    
#Replace the path of gtf and fasta files of all genomes (Ph. sojae, Ph. infestans, Pl. halstedii) in order to retreive upstream 400 nucleotides 
    annotation_gtf = '/path_to/Transcriptional_regulation_oomycetes/data/support/orthologs_motifs_study/Psoj.gtf' 
    genome_fasta = '/path_to/Transcriptional_regulation_oomycetes/data/support/orthologs_motifs_study/Psoj.fa'

    genome_handler(genome_fasta)
    gtf_handler(annotation_gtf)


    pp = pprint.PrettyPrinter(indent=4)

    output = open(("Psoj_promoter.fasta"), 'w')

    for genes in p_gtf_dict:

        upstream_generator(genes, output)
    output.close()


#Function to generate upstream promoter sequences

def upstream_generator(gene_id, output):



        try:
                seq = seq_dict[p_gtf_dict[gene_id]['scaffold_id']]

        except:
                print("Error with gene_id: %s" % (gene_id))


        strand = p_gtf_dict[gene_id]['strand']

        promoter_len = 400

        if (strand == '+'):

                if (int(p_gtf_dict[gene_id]['start']) >= promoter_len):
                        extract = seq [int(p_gtf_dict[gene_id]['start'])-promoter_len:int(p_gtf_dict[gene_id]['start'])-1]

                        seq_from = int(p_gtf_dict[gene_id]['start'])-promoter_len-1
                        seq_to = int(p_gtf_dict[gene_id]['start'])-1

                elif (int(p_gtf_dict[gene_id]['start']) < promoter_len):
                        extract = seq[0:int(p_gtf_dict[gene_id]['start'])-1]
                        seq_from = 0
                        seq_to = int(p_gtf_dict[gene_id]['start'])-1
                output.write(">%s (%s:%i-%i strand:%s pattern:%s)\n%s" % (gene_id, p_gtf_dict[gene_id]['scaffold_id'], seq_from, seq_to, strand, p_dict[clusters], extract))

        elif (strand == '-'):


                if ((len(seq)-int(p_gtf_dict[gene_id]['end'])) >= promoter_len):

                        extract = seq[int(p_gtf_dict[gene_id]['end'])+1:int(p_gtf_dict[gene_id]['end'])+promoter_len]
                        seq_from = int(p_gtf_dict[gene_id]['end'])+1
                        seq_to = int(p_gtf_dict[gene_id]['end'])+promoter_len

                elif ((len(seq)-int(p_gtf_dict[gene_id]['end'])) < promoter_len):
                        extract = seq[int(p_gtf_dict[gene_id]['end'])+1:len(seq)]
                        seq_from = int(p_gtf_dict[gene_id]['end'])+1
                        seq_to = len(seq)

                extract = reverse_complement(extract)
                output.write(">%s (%s:%i-%i strand:%s)\n%s" % (gene_id, p_gtf_dict[gene_id]['scaffold_id'], seq_from, seq_to, strand, p_dict[clusters], extract))


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
