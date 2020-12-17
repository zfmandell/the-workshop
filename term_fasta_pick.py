#!/usr/bin/env python

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

def main():
    parser = argparse.ArgumentParser(description='calculates the delta value at each position [ x upstream - x downstream].')
    parser.add_argument('ends',type=str,help='<.fasta> fwd fasta file of upstream sequences')
    parser.add_argument('fasta',type=str,help='<.fasta> rev fasta file of upstream sequences')
    parser.add_argument('outfile',type=str,help='name of outfile')
    args = parser.parse_args()

    peaks = []
    with open(args.ends,'r') as inp:
        for line in inp:
            peaks.append(int(line.strip().split(",")[0]))

    records = list(SeqIO.parse(args.fasta, "fasta"))

    to_keep = []
    for item in records:
        if int(item.id.split("_")[0].split(":")[1]) in peaks:
            to_keep.append(item)

    with open(args.outfile, "w") as output_handle:
        SeqIO.write(to_keep, output_handle, "fasta")

if __name__ == '__main__':
    main()
