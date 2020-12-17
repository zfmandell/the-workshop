#!/usr/bin/env python

"""
This script will take the output text file created by DESeq2, corresponding to the 3' end, Log2FC, and padj
will sort through file and filter based on user input criteria of minimum Log2FC and padj (default = abs(Log2FC) > 4 & padj < 0.001)
will also take genbank gff3 file corresponding to BSub Reference Genome
final return will be be two files:
text file containing all genes that contain differential 3' ends corresponding to +log2FC above thresholds, with metainformation
"""

from __future__ import division
import argparse
from collections import Counter
import operator

def read_gff3(gff3,strand):
    #reads in fasta, returns dict of sequences, dict of strand info for gene
    return_dict = {}
    strand_dict = {}
    with open(gff3,'r') as inp:
        for line in inp:
            if line[0] != '#' and len(line.strip().split()) != 0:
                if str(line.strip().split()[2]) == 'gene':
                    return_dict[(line.strip().split()[8].split(";")[2].split("=")[1]).replace("-",".")] = [int(line.strip().split()[3]),int(line.strip().split()[4])]
                    strand_dict[(line.strip().split()[8].split(";")[2].split("=")[1]).replace("-",".")] = line.strip().split()[6]
    if strand.lower() == 'plus':
        return sorted(return_dict.items(), key=operator.itemgetter(1)),strand_dict,return_dict
    elif strand.lower() == 'minus':
        return sorted(return_dict.items(), key=operator.itemgetter(1)),strand_dict,return_dict

def read_DE(DESeq2,FC,padj):
    #reads in DESeq2 textfile, returns dict of all ends that pass threshold parameters with associated log2FC
    final_ends = {}
    end_mean = {}
    with open(DESeq2,"r") as inp:
        firstline = inp.readline()
        for line in inp:
            if abs(float(line.strip().split()[3])) > FC and float(line.strip().split()[6]) < padj:
                final_ends[int(line.strip().split()[1])] = float(line.strip().split()[3])
    return final_ends

def read_rpkm(rpkm):
    #reads in rpkm txt file, return dict feature : rpkm
    # for intergenic features, will be sorted
    return_dict = {}
    with open(rpkm,'r') as inp:
        for line in inp:
            if line.strip().split("\t")[0].split("-") == 1:
                return_dict[line.strip().split("\t")[0]] = float(line.strip().split("\t")[1])
            else:
                return_dict["-".join(sorted(line.strip().split("\t")[0].split("-")))] = float(line.strip().split("\t")[1])
    return return_dict

def gene_compile(ends,seqs):
    #returns two lists, for each stop that passes threshold criteria, what gene it is in, pos if change is pos, neg if change is neg
    pos_genes = []
    neg_genes = []
    feature_dict = {}

    for key,value in ends.iteritems():
        if value > 0.0:
            for item in seqs:
                if item[1][0] <= key <= item[1][1]:
                    pos_genes.append(",".join([str(item[0]),'CDS']))
                elif seqs.index(item) < len(seqs) - 1:
                    if item[1][1] < key < seqs[seqs.index(item)+1][1][0]:
                        pos_genes.append(",".join(["-".join([item[0],seqs[seqs.index(item)+1][0]]),'inter']))
                elif seqs.index(item) == len(seqs) - 1:
                    if item[1][1] < key:
                        pos_genes.append(",".join(["-".join([item[0],seqs[0][0]]),'inter']))
        if value < 0.0:
            for item in seqs:
                if item[1][0] <= key <= item[1][1]:
                    neg_genes.append(",".join([str(item[0]),'CDS']))
                elif seqs.index(item) < len(seqs) - 1:
                    if item[1][1] < key < seqs[seqs.index(item)+1][1][0]:
                        neg_genes.append(",".join(["-".join([item[0],seqs[seqs.index(item)+1][0]]),'inter']))
                elif seqs.index(item) == len(seqs) - 1:
                    if item[1][1] < key:
                        neg_genes.append(",".join(["-".join([item[0],seqs[0][0]]),'inter']))

    for item in pos_genes:
        if item.split(",")[0] not in feature_dict.keys():
            feature_dict[item.split(",")[0]] = item.split(",")[1]

    for item in neg_genes:
        if item.split(",")[0] not in feature_dict.keys():
            feature_dict[item.split(",")[0]] = item.split(",")[1]

    return pos_genes,neg_genes,feature_dict

def writer(pos_genes,neg_genes,outfile,strands,strand,feature_dict,rpkm_dict,raw_dict):
    #writes two files, one for pos changes, one for neg changes, has gene, plus number of differential ends within gene
    #count of differential 3' ends are all normalized by rpkm as determined by coverage files, for each, mean of two biological conditions

    pos = sorted(Counter(pos_genes).items(), key=operator.itemgetter(1),reverse=True)
    neg = sorted(Counter(neg_genes).items(), key=operator.itemgetter(1),reverse=True)

    raw_pos = {}
    raw_neg = {}
    norm_count_pos = {}
    norm_count_neg = {}
    if strand.lower() == 'plus':
        for item in pos:
            features = item[0].split(",")
            if str(features[1]) == 'CDS':
                norm_count_pos[features[0]] = float(item[1])/rpkm_dict[features[0]]
                raw_pos[features[0]] = item[1]
            else:
                norm_count_pos[features[0]] = float(item[1])/rpkm_dict["-".join(sorted(features[0].split("-")))]
                raw_pos[features[0]] = item[1]
        for item in neg:
            features = item[0].split(",")
            if str(features[1]) == 'CDS':
                norm_count_neg[features[0]] = float(item[1])/rpkm_dict[features[0]]
                raw_neg[features[0]] = item[1]
            else:
                norm_count_neg[features[0]] = float(item[1])/rpkm_dict["-".join(sorted(features[0].split("-")))]
                raw_neg[features[0]] = item[1]
    if strand.lower() == 'minus':
        for item in pos:
            features = item[0].split(",")
            if str(features[1]) == 'CDS':
                norm_count_pos[features[0]] = float(item[1])/rpkm_dict[features[0]]
                raw_pos[features[0]] = item[1]
            else:
                norm_count_pos[features[0]] = float(item[1])/rpkm_dict["-".join(sorted(features[0].split("-")))]
                raw_pos[features[0]] = item[1]
        for item in neg:
            features = item[0].split(",")
            if str(features[1]) == 'CDS':
                norm_count_neg[features[0]] = float(item[1])/rpkm_dict[features[0]]
                raw_neg[features[0]] = item[1]
            else:
                norm_count_neg[features[0]] = float(item[1])/rpkm_dict["-".join(sorted(features[0].split("-")))]
                raw_neg[features[0]] = item[1]

    pos_norm = sorted(norm_count_pos.items(), key=operator.itemgetter(1),reverse=True)
    neg_norm = sorted(norm_count_neg.items(), key=operator.itemgetter(1),reverse=True)

    if strand.lower() == 'plus':
        with open("pos_"+outfile,"a+") as outp:
            outp.write("gene"+","+"feature"+","+"strand_of_feature"+","+"strand_of_3'end"+","+"normalized_count"+","+"raw_count"+",length\n")
            for item in pos_norm:
                if len(item[0].split("-")) == 1:
                    length = str(raw_dict[item[0]][1]-raw_dict[item[0]][0])
                    outp.write(str(item[0])+","+str(feature_dict[item[0]])+","+str(strands[item[0]])+","+"+"+","+str(item[1])+","+str(raw_pos[item[0]])+","+length+"\n")
                else:
                    genes = item[0].split("-")
                    length = str(raw_dict[genes[1]][0]-raw_dict[genes[0]][1])
                    outp.write(str(item[0])+","+str(feature_dict[item[0]])+","+":".join([strands[genes[0]],strands[genes[1]]])+","+"+"+","+str(item[1])+","+str(raw_pos[item[0]])+","+length+"\n")
        with open("neg_"+outfile,"a+") as outp:
            outp.write("gene"+","+"feature"+","+"strand_of_feature"+","+"strand_of_3'end"+","+"normalized_count"+","+"raw_count"+",length\n")
            for item in neg_norm:
                if len(item[0].split("-")) == 1:
                    length = str(raw_dict[item[0]][1]-raw_dict[item[0]][0])
                    outp.write(str(item[0])+","+str(feature_dict[item[0]])+","+str(strands[item[0]])+","+"+"+","+str(item[1])+","+str(raw_neg[item[0]])+","+length+"\n")
                else:
                    genes = item[0].split("-")
                    length = str(raw_dict[genes[1]][0]-raw_dict[genes[0]][1])
                    outp.write(str(item[0])+","+str(feature_dict[item[0]])+","+":".join([strands[genes[0]],strands[genes[1]]])+","+"+"+","+str(item[1])+","+str(raw_neg[item[0]])+","+length+"\n")

    elif strand.lower() == 'minus':
        with open("pos_"+outfile,"a+") as outp:
            outp.write("gene"+","+"feature"+","+"strand_of_feature"+","+"strand_of_3'end"+","+"normalized_count"+","+"raw_count"+",length\n")
            for item in pos_norm:
                if len(item[0].split("-")) == 1:
                    length = str(raw_dict[item[0]][1]-raw_dict[item[0]][0])
                    outp.write(str(item[0])+","+str(feature_dict[item[0]])+","+str(strands[item[0]])+","+"-"+","+str(item[1])+","+str(raw_pos[item[0]])+","+length+"\n")
                else:
                    genes = item[0].split("-")
                    length = str(raw_dict[genes[1]][0]-raw_dict[genes[0]][1])
                    outp.write(str(item[0])+","+str(feature_dict[item[0]])+","+":".join([strands[genes[1]],strands[genes[0]]])+","+"-"+","+str(item[1])+","+str(raw_pos[item[0]])+","+length+"\n")
        with open("neg_"+outfile,"a+") as outp:
            outp.write("gene"+","+"feature"+","+"strand_of_feature"+","+"strand_of_3'end"+","+"normalized_count"+","+"raw_count"+",length\n")
            for item in neg_norm:
                if len(item[0].split("-")) == 1:
                    length = str(raw_dict[item[0]][1]-raw_dict[item[0]][0])
                    outp.write(str(item[0])+","+str(feature_dict[item[0]])+","+str(strands[item[0]])+","+"-"+","+str(item[1])+","+str(raw_neg[item[0]])+","+length+"\n")
                else:
                    genes = item[0].split("-")
                    length = str(raw_dict[genes[1]][0]-raw_dict[genes[0]][1])
                    outp.write(str(item[0])+","+str(feature_dict[item[0]])+","+":".join([strands[genes[1]],strands[genes[0]]])+","+"-"+","+str(item[1])+","+str(raw_neg[item[0]])+","+length+"\n")

def main():
    parser = argparse.ArgumentParser(description='will create .csv file of all delta values at all called peaks ')
    parser.add_argument('DESeq2',type=str,help='<.txt> output text file of DESeq2')
    parser.add_argument('genbank',type=str,help='<.gff3> reference genome from genbank in gff3 format')
    parser.add_argument('rpkm',type=str,help='<.txt> rpkm file, output of rpkm_calc.py')
    parser.add_argument('outfile',type=str,default=None,help=' base name of output files')
    parser.add_argument('strand',type=str,default='plus',help='strand information <plus/minus> default = plus')
    parser.add_argument('-Log2FC',type=float,default=4.0,help='<.txt> Log2FC threshold')
    parser.add_argument('-padj',type=float,default=0.001,help='<.txt> padj threshold')
    args = parser.parse_args()

    seqs,strands,raw_dict = read_gff3(args.genbank,args.strand)
    ends = read_DE(args.DESeq2,args.Log2FC,args.padj)
    rpkm = read_rpkm(args.rpkm)
    pos,neg,features = gene_compile(ends,seqs)
    writer(pos,neg,args.outfile,strands,args.strand,features,rpkm,raw_dict)

if __name__ == '__main__':
    main()
