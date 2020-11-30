#!/usr/bin/env python

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC


def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('fwd',type=str,help='<.fasta> fwd fasta file of upstream sequences')
    parser.add_argument('rev',type=str,help='<.fasta> rev fasta file of upstream sequences')
    parser.add_argument('term_efficiency',type=str,help='<.csv> term efficiency csv for this strain')
    parser.add_argument('outfile',type=str,help='name of outfile')
    args = parser.parse_args()

    peaks = []
    with open(args.term_efficiency,'r') as inp:
        firstline = inp.readline()
        for line in inp:
            peaks.append(int(line.strip().split(",")[0]))

    records_fwd = list(SeqIO.parse(args.fwd, "fasta"))
    records_rev = list(SeqIO.parse(args.rev, "fasta"))

    with open(args.outfile,'w') as outp:
        for item in records_fwd:
            if int(item.id.split("_")[0].split(":")[1]) in peaks:
                outp.write('>'+str(item.id)+"\n"+str(item.seq)+"\n")
        for item in records_rev:
            if int(item.id.split("_")[0].split(":")[1]) in peaks:
                outp.write('>'+str(item.id)+"\n"+str(item.seq)+"\n")

if __name__ == '__main__':
    main()
