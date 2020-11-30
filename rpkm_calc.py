#!/usr/bin/env python

"""
This script will take a coverage file from an rnaseq run (.cov) NEEDS TO BE ALL, NOT OLIGO ONLY
Will also take a genbank .gff3 file, and will generate a RPKM value for each CDS sequence, and intergenic region

RPKM = 10^9 * (N/L) * (1/C)

N = # reads of feature
L = lengh of feature
C = total # of reads

output file will be tab delimited text file eg) feature RPKM
"""

from __future__ import division
import argparse
import numpy as np
import operator

def read_gff3(gff3):
    #reads in fasta, returns dict of sequences, dict of strand info for gene
    return_dict = {}
    with open(gff3,'r') as inp:
        for line in inp:
            if line[0] != '#' and len(line.strip().split()) != 0:
                if str(line.strip().split()[2]) == 'gene':
                    return_dict[(line.strip().split()[8].split(";")[2].split("=")[1]).replace("-",".")] = [int(line.strip().split()[3]),int(line.strip().split()[4])]
    return sorted(return_dict.items(), key=operator.itemgetter(1))

def read_cov(cov):
    #reads in a coverage file, loads all coverage values into a single list, index 1
    return_list = []
    with open(cov,'r') as inp:
        for line in inp:
            return_list.append(float(line.strip().split("\t")[2]))
    return return_list,sum(return_list)

def feature_place(gff3,cov):
    #will place each gene+intergenic feature into coverage context, return dict key = feature (gene if CDS or gene-gene if intra), value = length,coverage sum
    feature_dict = {}
    for item in gff3:
        length = item[1][1]-item[1][0]
        coverage = sum(cov[item[1][0]:item[1][1]+1])
        feature_dict[item[0]] = [length,coverage]
    for item in gff3:
        if gff3.index(item) < len(gff3)-1:
            length = gff3[gff3.index(item)+1][1][0]-gff3[gff3.index(item)][1][1]
            if length < 0:
                length = length*-1
                coverage = sum(cov[gff3[gff3.index(item)+1][1][0]+1:gff3[gff3.index(item)][1][1]])
            else:
                coverage = sum(cov[gff3[gff3.index(item)][1][1]:gff3[gff3.index(item)+1][1][0]+1])
            feature_dict["-".join([gff3[gff3.index(item)][0],gff3[gff3.index(item)+1][0]])] = [length,coverage]
    return feature_dict

def rpkm_calc(feature_dict,total_cov):
    #will calculate rpkm of each feature, return dict of feature : rpkm
    rpkm_dict = {}
    for key,value in feature_dict.iteritems():
        try:
            rpkm_dict[key] = 1000000000*(value[1]/value[0])*(1/total_cov)
        except ZeroDivisionError:
            rpkm_dict[key] = 0
    return rpkm_dict


def writer(rpkm_dict,outname):
    #output text file
    with open(outname,"w") as outp:
        for key,value in rpkm_dict.iteritems():
            outp.write(key+"\t"+str(value)+"\n")

def main():
    parser = argparse.ArgumentParser(description='will create .csv file of all delta values at all called peaks ')
    parser.add_argument('cov',type=str,help='<.cov> coverage file for RNA-Seq')
    parser.add_argument('genbank',type=str,help='<.gff3> reference genome from genbank in gff3 format')
    parser.add_argument('outfile',type=str,default=None,help=' base name of output files')
    args = parser.parse_args()

    gff3 = read_gff3(args.genbank)
    cov,total_cov = read_cov(args.cov)
    placed = feature_place(gff3,cov)
    rpkm = rpkm_calc(placed,total_cov)
    writer(rpkm,args.outfile)

if __name__ == '__main__':
    main()
