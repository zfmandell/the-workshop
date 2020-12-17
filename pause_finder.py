#!/usr/bin/env python

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

def build_end_dict(fyle):
    return_dyct = {}
    with open(fyle,'r') as inp:
        firstline = inp.readline()
        for line in inp:
            return_dyct[int(line.strip().split(",")[0])] = [line.strip().split(",")[1],line.strip().split(",")[3]]
    return return_dyct

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

def complement(seq):
    alt_map = {'ins':'0'}
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    for k,v in alt_map.iteritems():
        seq = seq.replace(k,v)
    bases = list(seq)
    bases = [complement.get(base,base) for base in bases]
    bases = ''.join(bases)
    for k,v in alt_map.iteritems():
        bases = bases.replace(v,k)
    return bases

def U_find(seq,U):
    return seq.replace('T','U').find(U)

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('ends',type=str,help='')
    parser.add_argument('genome',type=str,help='')
    parser.add_argument('outfile',type=str,help='')
    args = parser.parse_args()

    peaks,fasta_fyle = build_end_dict(args.ends),read_fasta(args.genome)


    with open(args.outfile,"w") as outp:
        for key,value in peaks.iteritems():
            U = value[1]
            if value[0] == '+':
                base_seq = fasta_fyle['fasta'][key-700:key+700]
                outp.write(">position:"+str(key)+"_strand:"+str(value[0])+"\n")
                if U_find(base_seq,U) == -1:
                    print 'no'
                    print '+'
                    print key
                else:
                    end = U_find(base_seq,U)+len(U)+501
                    start = U_find(base_seq,U)+len(U)
                    outp.write(base_seq[start:end])
                    outp.write("\n")
            else:
                base_seq = reverse_complement(fasta_fyle['fasta'][key-700:key+700])
                outp.write(">position:"+str(key)+"_strand:"+str(value[0])+"\n")
                if U_find(base_seq,U) == -1:
                    print 'no'
                    print '-'
                    print key
                else:
                    end = U_find(base_seq,U)+len(U)+501
                    start = U_find(base_seq,U)+len(U)
                    outp.write(base_seq[start:end])
                    outp.write("\n")


if __name__ == '__main__':
    main()
