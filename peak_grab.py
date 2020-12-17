#!/usr/bin/env python

"""
This script will take the fasta file created by local_peak.py, and all of the .bedgraph slope replicate files in a directory
Output will be the delta values for each replicate at each peak in .csv format
Intended to create count matrix for DESeq2 -> differential expression of 3' ends

Can only hand two biological conditions, ex) WT and KO, WT and Treated

"""

import argparse
import glob
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import linecache


def delta_grab(fastas):
    #accepts fasta file, will output list of position of all called 3' ends, sorted in ascending order
    deltas = []
    for item in fastas:
        fasta_sequences = SeqIO.parse(open(item),'fasta')
        for seq in fasta_sequences:
            deltas.append(int(str(seq.id).split("_")[0].split(":")[1].split(".")[0]))
    return sorted(set(deltas),reverse=False)

def writer(delta_list,fyle_list,output):
    #applies delta list to coverage files, returns CSV corresponding to the delta value at called 3' end in all replicate coverage files
    with open(output,"w") as outp:
        outp.write("pos,")
        for item in fyle_list[:-1]:
            outp.write(str(item)+",")
        outp.write(str(fyle_list[-1])+"\n")
        for item in delta_list:
            outp.write(str(item)+",")
            for cov in fyle_list[:-1]:
                line = linecache.getline(cov,item-1)
                outp.write(str(line.strip().split("\t")[3])+",")
            line = linecache.getline(fyle_list[-1],item-1)
            outp.write(str(line.strip().split("\t")[3])+"\n")

def main():
    parser = argparse.ArgumentParser(description='will create .csv file of all delta values at all called peaks ')
    parser.add_argument('fasta_WT',type=str,help='<.fasta> WT fasta file to pull sequences from ')
    parser.add_argument('fasta_EXP',type=str,help='<.fasta> Experimental fasta file to pull sequences from')
    parser.add_argument('outfile',type=str,default=None,help='name of output file')
    args = parser.parse_args()

    deltas = delta_grab([args.fasta_WT,args.fasta_EXP])

    fyle_list = []
    for fyle in sorted(glob.glob('*.bedgraph')):
        fyle_list.append(fyle)
    writer(deltas,fyle_list,args.outfile)

if __name__ == '__main__':
    main()
