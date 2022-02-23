#!/usr/bin/env python3

import sys
from Bio import SeqIO

def read_in_fasta(afasta):
    '''Reads in a fasta file to a dictionary'''
    fasta_dict = {}
    fasta_sequences = SeqIO.parse(open(afasta),'fasta')
    for fasta in fasta_sequences:
        fasta_dict[fasta.id] = str(fasta.seq).upper()
    return fasta_dict

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
    bases = list(seq)
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    return bases

fastas = read_in_fasta(sys.argv[1])
with open(sys.argv[2],'w') as outp:
    for key,value in fastas.items():
        outp.write('>'+key+'\n'+reverse_complement(value)+'\n')