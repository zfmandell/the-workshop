#!/usr/bin/env python

"""
This script will take the output text file created by DESeq2, corresponding to the 3' end, Log2FC, and padj
will sort through file and filter based on user input criteria of minimum Log2FC and padj (default = abs(Log2FC) > 4 & padj < 0.001)
will also take genbank gff3 file corresponding to BSub Reference Genome
final return will be a file corresponding to the genomic coordinates of the DE stop, and the percentage into the feature the DE stop is
"""


from __future__ import division
import argparse
from collections import Counter
import operator

def read_gff3(gff3):
    #reads in fasta, returns dict of sequences, dict of strand info for gene
    return_dict = {}
    strand_dict = {}
    with open(gff3,'r') as inp:
        for line in inp:
            if line[0] != '#' and len(line.strip().split()) != 0:
                if str(line.strip().split()[2]) == 'gene':
                    return_dict[(line.strip().split()[8].split(";")[2].split("=")[1]).replace("-",".")] = [int(line.strip().split()[3]),int(line.strip().split()[4])]
                    strand_dict[(line.strip().split()[8].split(";")[2].split("=")[1]).replace("-",".")] = line.strip().split()[6]
    return sorted(return_dict.items(), key=operator.itemgetter(1))


def read_DE(DESeq2,FC,padj):
    #reads in DESeq2 textfile, returns dict of all ends that pass threshold parameters with associated log2FC
    final_ends = {}
    end_mean = {}
    with open(DESeq2,"r") as inp:
        firstline = inp.readline()
        for line in inp:
            if abs(float(line.strip().split()[3])) > FC and float(line.strip().split()[6]) < padj:
                final_ends[int(line.strip().split()[1])] = float(line.strip().split()[3])
                end_mean[int(line.strip().split()[1])] = float(line.strip().split()[2])
    return final_ends

def gene_compile(ends,seqs):
    #returns two lists, for each stop that passes threshold criteria, what feauture it is in, and coords of that feature
    pos_genes = {}
    neg_genes = {}

    for key,value in ends.iteritems():
        if value > 0.0:
            for item in seqs:
                if item[1][0] <= key <= item[1][1]:
                    pos_genes[key] = [item[0],item[1]]
                elif seqs.index(item) < len(seqs) - 1:
                    if item[1][1] < key < seqs[seqs.index(item)+1][1][0]:
                        pos_genes[key] = ["-".join([item[0],seqs[seqs.index(item)+1][0]]),[item[1][1],seqs[seqs.index(item)+1][1][0]]]
                elif seqs.index(item) == len(seqs) - 1:
                    if item[1][1] < key:
                        pos_genes[key] = ["-".join([item[0],seqs[0][0]]),[item[1][1],seqs[0][1][0]]]
        if value < 0.0:
            for item in seqs:
                if item[1][0] <= key <= item[1][1]:
                    neg_genes[key] = [item[0],item[1]]
                elif seqs.index(item) < len(seqs) - 1:
                    if item[1][1] < key < seqs[seqs.index(item)+1][1][0]:
                        neg_genes[key] = ["-".join([item[0],seqs[seqs.index(item)+1][0]]),[item[1][1],seqs[seqs.index(item)+1][1][0]]]
                elif seqs.index(item) == len(seqs) - 1:
                    if item[1][1] < key:
                        neg_genes[key] = ["-".join([item[0],seqs[0][0]]),[item[1][1],seqs[0][1][0]]]

    return pos_genes,neg_genes

def perc_calc(pos,neg,strand):
    #takes two dicts, corresponding to output of gene_compile, returns the percent into feature for each DE end in each list
    pos_perc = {}
    neg_perc = {}

    if strand.lower() == 'plus':
        for key,value in pos.iteritems():
            pos_perc[key] = [value[0],(float(key)-value[1][0])/(value[1][1]-value[1][0])*100]
        for key,value in neg.iteritems():
            neg_perc[key] = [value[0],(float(key)-value[1][0])/(value[1][1]-value[1][0])*100]
    elif strand.lower() == 'minus':
        for key,value in pos.iteritems():
            pos_perc[key] = [value[0],100-((float(key)-value[1][0])/(value[1][1]-value[1][0])*100)]
        for key,value in neg.iteritems():
            neg_perc[key] = [value[0],100-((float(key)-value[1][0])/(value[1][1]-value[1][0])*100)]

    return pos_perc,neg_perc

def writer(pos,neg,outname):
    #outputs to a file, DE End, Feature,perc into feature

    with open("pos_"+outname,"w") as outp:
        outp.write("DE_End,feature,%_into_feature\n")
        for key,value in pos.iteritems():
            if len(value[0].split("-")) == 1:
                outp.write(str(key)+',CDS,'+str(value[1])+"\n")
            else:
                outp.write(str(key)+',inter,'+str(value[1])+"\n")

    with open("neg_"+outname,"w") as outp:
        outp.write("DE_End,feature,%_into_feature\n")
        for key,value in neg.iteritems():
            if len(value[0].split("-")) == 1:
                outp.write(str(key)+',CDS,'+str(value[1])+"\n")
            else:
                outp.write(str(key)+',inter,'+str(value[1])+"\n")

def main():
    parser = argparse.ArgumentParser(description='will create .csv file of all delta values at all called peaks ')
    parser.add_argument('DESeq2',type=str,help='<.txt> output text file of DESeq2')
    parser.add_argument('genbank',type=str,help='<.gff3> reference genome from genbank in gff3 format')
    parser.add_argument('outfile',type=str,default=None,help=' base name of output files')
    parser.add_argument('strand',type=str,default=None,help='strand information <plus/minus>')
    parser.add_argument('-Log2FC',type=float,default=4.0,help='<.txt> Log2FC threshold')
    parser.add_argument('-padj',type=float,default=0.001,help='<.txt> padj threshold')
    args = parser.parse_args()

    seqs = read_gff3(args.genbank)
    ends = read_DE(args.DESeq2,args.Log2FC,args.padj)
    pos,neg = gene_compile(ends,seqs)
    pos_perc,neg_perc = perc_calc(pos,neg,args.strand)
    writer(pos_perc,neg_perc,args.outfile)

if __name__ == '__main__':
    main()
