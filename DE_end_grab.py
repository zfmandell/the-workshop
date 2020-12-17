#!/usr/bin/env python

"""
This script will take the output text file created by DESeq2, corresponding to the 3' end, Log2FC, and padj
will sort through file and filter based on user input criteria of minimum Log2FC and padj (default = abs(Log2FC) > 4 & padj < 0.001)
will also take fasta file output from local_peak.py
final return will be be two files:
all upstream fasta sequences corresponding to +log2FC above thresholds
all upstream fasta sequences corresponding to -log2FC above thresholds
"""


import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

def read_fasta(genome_fasta):
    #reads in fasta, returns dict of sequences
    fasta_sequences,fasta_dict =SeqIO.parse(open(genome_fasta),'fasta'),{}
    for fasta in fasta_sequences:
        fasta_dict[fasta.id] = str(fasta.seq)
    return fasta_dict

def read_DE(DESeq2,FC,padj):
    #reads in DESeq2 textfile, returns dict of all ends that pass threshold parameters with associated log2FC
    final_ends = {}
    with open(DESeq2,"r") as inp:
        firstline = inp.readline()
        for line in inp:
            if abs(float(line.strip().split()[3])) > FC and float(line.strip().split()[6]) < padj:
                final_ends[line.strip().split()[1]] = float(line.strip().split()[3])
    return final_ends

def main():
    parser = argparse.ArgumentParser(description='will create .csv file of all delta values at all called peaks ')
    parser.add_argument('DESeq2',type=str,help='<.txt> output text file of DESeq2')
    parser.add_argument('fasta',type=str,help='<.fasta> Experimental fasta file to pull sequences from')
    parser.add_argument('outfile',type=str,default=None,help=' base name of output files')
    parser.add_argument('-Log2FC',type=float,default=4.0,help='<.txt> Log2FC threshold')
    parser.add_argument('-padj',type=float,default=0.001,help='<.txt> padj threshold')
    args = parser.parse_args()

    seqs = read_fasta(args.fasta)
    ends = read_DE(args.DESeq2,args.Log2FC,args.padj)

    for key,value in ends.iteritems():
        if value > 0.0:
            for k,v in seqs.iteritems():
                if str(key) == str(k.split("_")[0].split(":")[1].split(".")[0]):
                    with open("pos_"+args.outfile,"a+") as outp:
                        outp.write(">"+k+"\n")
                        outp.write(seqs[k]+"\n")

        if value < 0.0:
            for k,v in seqs.iteritems():
                if key == k.split("_")[0].split(":")[1].split(".")[0]:
                    with open("neg_"+args.outfile,"a+") as outp:
                        outp.write(">"+k+"\n")
                        outp.write(seqs[k]+"\n")



if __name__ == '__main__':
    main()
