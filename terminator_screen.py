#!/usr/bin/env python

"takes a list of WT terminators, filters creates a fasta file from the two provided fasta files composed of all sequences that are not terminators"


import sys
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC


def term_build(terms):
    term = []
    with open(terms,"r") as inp:
        firstline = inp.readline()
        secondline = inp.readline()
        for line in inp:
            if line.strip().split(",")[2] == 'Terminator':
                term.append(float(line.strip().split(",")[0]))
    return list(set(term))

def read_fasta(genome_fasta):
    fasta_sequences,fasta_dict =SeqIO.parse(open(genome_fasta),'fasta'),{}
    for fasta in fasta_sequences:
        fasta_dict[fasta.id] = str(fasta.seq)
    return fasta_dict

def fasta_prune(terms,fasta):
    shared = []
    for item in fasta.keys():
        for term in terms:
            if (float(term)-3.0 <= float(item.split("_")[0].split(":")[1]) <= float(term)+3.0):
                shared.append(item)
    return list(set(shared))

def main():
    parser = argparse.ArgumentParser(description='will create .fasta file corresponding to all non terminator 3 prime ends')
    parser.add_argument('terminators',type=str,default=None,help='<.txt> terminator list')
    parser.add_argument('fasta_fwd',type=str,default=None,help='<.txt> fwd fasta file')
    parser.add_argument('fasta_rev',type=str,default=None,help=' rev fasta file')
    parser.add_argument('outfile',type=str,default=None,help='output file name in fasta format')
    args = parser.parse_args()

    terms = term_build(args.terminators)

    fasta_fwd = read_fasta(args.fasta_fwd)
    fasta_rev = read_fasta(args.fasta_rev)

    term_fwd = fasta_prune(terms,fasta_fwd)
    term_rev = fasta_prune(terms,fasta_rev)

    with open(args.outfile,'w') as outp:
        for key,value in fasta_fwd.iteritems():
            if key not in term_fwd:
                outp.write('>'+str(key)+'\n')
                outp.write(str(value)+"\n")

        for key,value in fasta_rev.iteritems():
            if key not in term_rev:
                outp.write('>'+str(key)+'\n')
                outp.write(str(value)+"\n")

if __name__ == '__main__':
    main()
