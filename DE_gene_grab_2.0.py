#!/usr/bin/env python

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC


def read_fasta(genome_fasta):
    #reads in fasta
    fasta_sequences,fasta_dict =SeqIO.parse(open(genome_fasta),'fasta'),{}
    for fasta in fasta_sequences:
        fasta_dict['fasta'] = str(fasta.seq)
    return fasta_dict

def reverse_complement(seq):
    alt_map = {'ins':'0'}
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    for k,v in alt_map.iteritems():
        seq = seq.replace(k,v)
    bases = list(seq)
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    for k,v in alt_map.iteritems():
        bases = bases.replace(v,k)
    return bases

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
    parser.add_argument('strand',type=str,default='plus',help='strand information <plus/minus> default = plus')
    parser.add_argument('-Log2FC',type=float,default=4.0,help='<.txt> Log2FC threshold')
    parser.add_argument('-padj',type=float,default=0.001,help='<.txt> padj threshold')
    args = parser.parse_args()

    fasta = read_fasta(args.fasta)
    ends = read_DE(args.DESeq2,args.Log2FC,args.padj)

    if args.strand.lower() == 'plus':
        for key,value in ends.iteritems():
            if value > 0.0:
                with open("pos_"+args.outfile,'a+') as outp:
                    outp.write(">position:"+str(key)+"_DE:"+str(value)+"\n")
                    outp.write(str(fasta['fasta'][(int(key)-50):int(key)+1])+"\n")
            if value < 0.0:
                with open("neg_"+args.outfile,'a+') as outp:
                    outp.write(">position:"+str(key)+"_DE:"+str(value)+"\n")
                    outp.write(str(fasta['fasta'][(int(key)-50):int(key)+1])+"\n")

    if args.strand.lower() == 'minus':
        for key,value in ends.iteritems():
            if value > 0.0:
                with open("pos_"+args.outfile,'a+') as outp:
                    outp.write(">position:"+str(key)+"_DE:"+str(value)+"\n")
                    outp.write(reverse_complement(str(fasta['fasta'][(int(key)):int(key)+51]))+"\n")

            if value < 0.0:
                with open("neg_"+args.outfile,'a+') as outp:
                    outp.write(">position:"+str(key)+"_DE:"+str(value)+"\n")
                    outp.write(reverse_complement(str(fasta['fasta'][(int(key)):int(key)+51]))+"\n")

if __name__ == '__main__':
    main()
