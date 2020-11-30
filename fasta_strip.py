#!/usr/bin/env python

import sys
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

fasta = read_fasta(sys.argv[1])

ends = {}
with open(sys.argv[2],'r') as inp:
    firstline = inp.readline()
    for line in inp:
        ends[int(line.strip().split(',')[0])] = line.strip().split(',')[1]

with open(sys.argv[3],'w') as outp:
    for key,value in ends.iteritems():
        if value == '+':
            outp.write('>'+str(key)+'_'+str(value)+'\n'+fasta['fasta'][key-20:key+21]+'\n')
        else:
            outp.write('>'+str(key)+'_'+str(value)+'\n'+reverse_complement(fasta['fasta'][key-20:key+21])+'\n')
